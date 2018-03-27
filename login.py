# -*- coding: utf-8 -*-
from __future__ import division
import bempp.api
import bempp.api as bem
import numpy as np
import scipy.sparse.linalg
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator
from bempp.api.grid import grid_from_element_data
from scipy.sparse.linalg.interface import LinearOperator as _LinearOperator
from bempp.api.utils.logging import timeit as _timeit
from bempp.api.assembly.discrete_boundary_operator import DiscreteBoundaryOperator, SparseDiscreteBoundaryOperator

def get_h(grid):
    #Routine that returns the Mesh Quality Information
    elements1 = list(grid.leaf_view.entity_iterator(1))
    vol = []
    for el1 in elements1:
        vol.append(el1.geometry.volume)
    vol = np.array(vol)
    return [vol.min(), vol.max(), vol.mean()]

def translate(grid, a=0, b=0, c=0):
    #Routine that translates a grid by a vector (a,b,c)
    vertices = grid.leaf_view.vertices.T
    vertices[:, 0] += a
    vertices[:, 1] += b
    vertices[:, 2] += c
    elements = grid.leaf_view.elements
    return grid_from_element_data(vertices.T,elements)


class _it_counter(object):
    def __init__(self, store_residuals):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []

    def __call__(self, x):
        self._count += 1
        if self._store_residuals:
            self._residuals.append(np.linalg.norm(x))
            print ('iteration -', self._count,"|| residual -", np.linalg.norm(x))


    @property
    def count(self):
        return self._count

    @property
    def residuals(self):
        return self._residuals

"""
Modification of linalg/iterative_solvers/gmres function: Allows to plug numpy matrices, LinearOperator or weak_forms() as an input
"""

def gmres(A, b, tol=1E-5, restart=None, maxiter=None, use_strong_form=False, return_residuals=False):
    """Interface to the scipy.sparse.linalg.gmres function.

    This function behaves like the scipy.sparse.linalg.gmres function. But
    instead of a linear operator and a vector b it can take a boundary operator
    and a grid function. In this case, the result is returned as a grid function in the
    correct space.
    
    """
    import bempp.api
    import time

    if not isinstance(A, BoundaryOperator) and use_strong_form==True:
        raise ValueError("Strong form only with BoundaryOperator")
    if isinstance(A, BoundaryOperator) and not isinstance(b, GridFunction):
        raise ValueError("Instance Error")
    if not isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        raise ValueError("Instance Error")

    # Assemble weak form before the logging messages

    if isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        if use_strong_form:
            if not A.range.is_compatible(b.space):
                raise ValueError(
                    "The range of A and the space of A must have the same number of unknowns if the strong form is used.")
            A_op = A.strong_form()
            b_vec = b.coefficients
        else:
            A_op = A.weak_form()
            b_vec = b.projections(A.dual_to_range)
    else:
        A_op = A
        b_vec = b

    callback = _it_counter(return_residuals)

    
    start_time = time.time()
    x, info = scipy.sparse.linalg.gmres(A_op, b_vec,
                                        tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    
    if isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        res_fun = GridFunction(A.domain, coefficients=x.ravel())

    else:
        res_fun = x
    if return_residuals:
        return res_fun, info, callback.residuals
    else:
        return res_fun, info


"""
Modification of assembly.assembler,assemble_singular_part function: Added a distance parameter to consider more interactions. Highly sequential assembly, the code could be improved
"""

def assemble_singular_part(operator, distance=1):
    """Assemble the singular part of an integral operator."""

    from scipy.sparse import csc_matrix
    from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator



    test_space = operator.dual_to_range
    trial_space = operator.domain

    test_dof_count = operator.dual_to_range.global_dof_count
    trial_dof_count = operator.domain.global_dof_count 

    # If the test and trial grid are different return a zero matrix

    #if operator.domain.grid != operator.dual_to_range.grid:
        #return SparseDiscreteBoundaryOperator(csc_matrix((test_dof_count, trial_dof_count)))

    grid = operator.domain.grid

    # Now get adjacent element pairs

    vertex_to_element_matrix = grid.leaf_view.vertex_to_element_matrix
    
    element_to_element_matrix = vertex_to_element_matrix.transpose().dot(vertex_to_element_matrix)

    element_to_element_matrix = np.power(element_to_element_matrix, distance)
    #element_to_element_matrix2 = np.dot(element_to_element_matrix,element_to_element_matrix)
    nonzero_pairs = element_to_element_matrix.nonzero()
    index_pairs = zip(nonzero_pairs[0], nonzero_pairs[1])

    # Now get all pairs of basis functions who partially overlap via adjacent elements

    all_test_trial_function_pairs = []

    for pair in index_pairs:

        test_element = grid.leaf_view.element_from_index(pair[0])
        trial_element = grid.leaf_view.element_from_index(pair[1])

        global_test_dofs = test_space.get_global_dofs(test_element)
        global_trial_dofs = trial_space.get_global_dofs(trial_element)

        for test_dof_index in global_test_dofs:
            if test_dof_index > -1:
                for trial_dof_index in global_trial_dofs:
                    if trial_dof_index > -1:
                        all_test_trial_function_pairs.append((test_dof_index, trial_dof_index))
        
    # Remove duplicates 
    all_test_trial_function_pairs = list(set(all_test_trial_function_pairs))

    # Now get all integraton element pairs associated 

    all_integration_element_pairs = []

    for function_pair in all_test_trial_function_pairs:

        local_dofs = test_space.global_to_local_dofs(function_pair)[0]
        test_local_dofs = local_dofs[0]
        trial_local_dofs = local_dofs[1]

        for test_dof in test_local_dofs:
            for trial_dof in trial_local_dofs:
                all_integration_element_pairs.append(
                        (test_dof.entity_index, trial_dof.entity_index))

     # Remove duplicates
    all_integration_element_pairs = list(set(all_integration_element_pairs))

    # Now compute all local dof interactions

    assembler = operator.local_assembler

    weak_forms = assembler.evaluate_local_weak_forms(all_integration_element_pairs)

    # Create a dictionary

    weak_form_lookup = dict(zip(all_integration_element_pairs, weak_forms))

    # Now need to create the sparse matrix

    data = np.zeros(len(all_test_trial_function_pairs), dtype=assembler.dtype)
    row_indices = np.zeros(len(all_test_trial_function_pairs), dtype=np.int)
    col_indices = np.zeros(len(all_test_trial_function_pairs), dtype=np.int)

    #import ipdb; ipdb.set_trace()

    for index, (test_index, trial_index) in enumerate(all_test_trial_function_pairs):

        row_indices[index] = test_index
        col_indices[index] = trial_index

        # Get local dofs and dof_weights
        local_dofs, weights = test_space.global_to_local_dofs([test_index, trial_index])
        test_dofs, trial_dofs = local_dofs
        test_weights, trial_weights = weights

        for i, test_dof in enumerate(test_dofs):
            for j, trial_dof in enumerate(trial_dofs):
                element_pair = (test_dof.entity_index, trial_dof.entity_index)
                data[index] += (weak_form_lookup[element_pair][test_dof.dof_index, trial_dof.dof_index] *
                        np.conj(test_weights[i]) * np.conj(trial_weights[j]))

    return SparseDiscreteBoundaryOperator(csc_matrix((data, (row_indices, col_indices)), 
            shape=(test_dof_count, trial_dof_count)))


class InverseSparseDiscreteBoundaryOperator(DiscreteBoundaryOperator):
    """
    Apply the (pseudo-)inverse of a sparse operator.

    This class uses a Sparse LU-Decomposition (in the case of a square matrix)
    or a sparse normal equation to provide the application of an inverse to
    a sparse operator.

    This class derives from
    :class:`scipy.sparse.linalg.interface.LinearOperator`
    and thereby implements the SciPy LinearOperator protocol.

    Parameters
    ----------
    operator : bempp.api.SparseDiscreteBoundaryOperator
        Sparse operator to be inverted.

    """


    def __init__(self, operator, mu=1e-4):

        self._solver = _Solver(operator, mu=mu)
        self._adjoint_op = None
        self._operator = operator
        super(InverseSparseDiscreteBoundaryOperator, self).__init__(
            self._solver.dtype, self._solver.shape)

    @_timeit("Inverse Sparse Operator matvec ")
    def _matvec(self, vec):
        """Implemententation of matvec."""

        return self._solver.solve(vec)

    def _rmatvec(self, vec):
        """Implemententation of rmatvec."""
        #pylint: disable=protected-access
        if self._adjoint_op is None:
            self._adjoint_op = self.adjoint()
        return self._adjoint_op * vec

    def _transpose(self):
        return InverseSparseDiscreteBoundaryOperator(
            self._operator.transpose())

    def _adjoint(self):
        return InverseSparseDiscreteBoundaryOperator(
            self._operator.adjoint())

    def elementary_operators(self):
        """Return the elementary operators that make up this operator."""
        return self._operator.elementary_operators()



class _Solver(object):  # pylint: disable=too-few-public-methods
    """Actual solve of a sparse linear system."""

    #pylint: disable=too-many-locals
    def __init__(self, operator, mu):

        from scipy.sparse import csc_matrix


        if isinstance(operator, SparseDiscreteBoundaryOperator):
            mat = operator.sparse_operator
        elif isinstance(operator, csc_matrix):
            mat = operator
        else:
            raise ValueError("op must be either of type " +
                             "SparseDiscreteBoundaryOperator or of type " +
                             "csc_matrix. Actual type: " +
                             str(type(operator)))

        from scipy.sparse.linalg import splu, spilu
        self._solve_fun = None
        self._shape = (mat.shape[1], mat.shape[0])
        self._dtype = mat.dtype
        
        import time
        import bempp.api

        SolverInterface = spilu
        actual_mat = mat

        start_time = time.time()
        
        # Square matrix case
        solver = SolverInterface(actual_mat, drop_tol=mu)
        self._solve_fun = solver.solve
        
        end_time = time.time()


    def solve(self, vec):
        """Solve with right-hand side vec."""

        if self._dtype == 'float64' and _np.iscomplexobj(vec):
            return (self.solve(_np.real(vec)) +
                    1j * self.solve(_np.imag(vec)))

        result = self._solve_fun(vec.squeeze())

        if vec.ndim > 1:
            return result.reshape(self.shape[0], 1)
        else:
            return result

    @property
    def shape(self):
        """Return the shape of the inverse operator."""
        return self._shape

    @property
    def dtype(self):
        """Return the dtype."""
        return self._dtype


def assemble_diagonal(operator):
    """Assemble the singular part of an integral operator."""

    from scipy.sparse import csc_matrix
    from bempp.api.assembly.discrete_boundary_operator import SparseDiscreteBoundaryOperator



    test_space = operator.dual_to_range
    trial_space = operator.domain

    test_dof_count = operator.dual_to_range.global_dof_count
    trial_dof_count = operator.domain.global_dof_count 

    # If the test and trial grid are different return a zero matrix

    #if operator.domain.grid != operator.dual_to_range.grid:
        #return SparseDiscreteBoundaryOperator(csc_matrix((test_dof_count, trial_dof_count)))

    grid = operator.domain.grid

    # Now get adjacent element pairs

    vertex_to_element_matrix = grid.leaf_view.vertex_to_element_matrix
    
    element_to_element_matrix = vertex_to_element_matrix.transpose().dot(vertex_to_element_matrix)

    #element_to_element_matrix2 = np.dot(element_to_element_matrix,element_to_element_matrix)
    nonzero_pairs = element_to_element_matrix.nonzero()
    index_pairs = zip(nonzero_pairs[0], nonzero_pairs[1])

    # Now get all pairs of basis functions who partially overlap via adjacent elements

    all_test_trial_function_pairs = []

    for pair in index_pairs:

        test_element = grid.leaf_view.element_from_index(pair[0])
        trial_element = grid.leaf_view.element_from_index(pair[1])

        global_test_dofs = test_space.get_global_dofs(test_element)
        global_trial_dofs = trial_space.get_global_dofs(trial_element)

        #for test_dof_index in global_test_dofs:
            #if test_dof_index > -1:
                #for trial_dof_index in global_trial_dofs:
                    #if trial_dof_index > -1:
                        #all_test_trial_function_pairs.append((test_dof_index, trial_dof_index))
        for test_dof_index in global_test_dofs:
            #for trial_dof_index in global_trial_dofs:
            if test_dof_index > -1:
                all_test_trial_function_pairs.append((test_dof_index, test_dof_index))
    # Remove duplicates 
    all_test_trial_function_pairs = list(set(all_test_trial_function_pairs))

    # Now get all integraton element pairs associated 

    all_integration_element_pairs = []

    for function_pair in all_test_trial_function_pairs:

        local_dofs = test_space.global_to_local_dofs(function_pair)[0]
        test_local_dofs = local_dofs[0]
        trial_local_dofs = local_dofs[1]

        for test_dof in test_local_dofs:
            for trial_dof in trial_local_dofs:
                all_integration_element_pairs.append(
                        (test_dof.entity_index, trial_dof.entity_index))

     # Remove duplicates
    all_integration_element_pairs = list(set(all_integration_element_pairs))

    # Now compute all local dof interactions

    assembler = operator.local_assembler

    weak_forms = assembler.evaluate_local_weak_forms(all_integration_element_pairs)

    # Create a dictionary

    weak_form_lookup = dict(zip(all_integration_element_pairs, weak_forms))

    # Now need to create the sparse matrix

    data = np.zeros(len(all_test_trial_function_pairs), dtype=assembler.dtype)
    row_indices = np.zeros(len(all_test_trial_function_pairs), dtype=np.int)
    col_indices = np.zeros(len(all_test_trial_function_pairs), dtype=np.int)

    #import ipdb; ipdb.set_trace()

    for index, (test_index, trial_index) in enumerate(all_test_trial_function_pairs):

        row_indices[index] = test_index
        col_indices[index] = trial_index

        # Get local dofs and dof_weights
        local_dofs, weights = test_space.global_to_local_dofs([test_index, trial_index])
        test_dofs, trial_dofs = local_dofs
        test_weights, trial_weights = weights

        for i, test_dof in enumerate(test_dofs):
            for j, trial_dof in enumerate(trial_dofs):
                element_pair = (test_dof.entity_index, trial_dof.entity_index)
                data[index] += (weak_form_lookup[element_pair][test_dof.dof_index, trial_dof.dof_index] *
                        np.conj(test_weights[i]) * np.conj(trial_weights[j]))

    return SparseDiscreteBoundaryOperator(csc_matrix((data, (row_indices, col_indices)), 
            shape=(test_dof_count, trial_dof_count)))
