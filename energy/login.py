# -*- coding: utf-8 -*-
from __future__ import division
import bempp.api as bem
import numpy as np
import scipy.sparse.linalg
from bempp.api.assembly import GridFunction
from bempp.api.assembly import BoundaryOperator
from bempp.api.grid import grid_from_element_data

def get_h(grid):
    # routine that gives the exact value of h defined as the maximum length of the triangles
    # of a given mesh
    elements1 = list(grid.leaf_view.entity_iterator(1))
    vol = []
    for el1 in elements1:
        vol.append(el1.geometry.volume)
    vol = np.array(vol)
    return [vol.min(), vol.max(), vol.mean()]

class _it_counter(object):

    def __init__(self, store_residuals):
        self._count = 0
        self._store_residuals = store_residuals
        self._residuals = []

    def __call__(self, x):
        self._count += 1
        if self._store_residuals:
            self._residuals.append(np.linalg.norm(x))
            print('iteration -', self._count,"|| residual -", np.linalg.norm(x))


    @property
    def count(self):
        return self._count

    @property
    def residuals(self):
        return self._residuals

    
def gmres(A, b, tol=1E-5, restart=None, H=None, solver = 'scipy', maxiter=None, use_strong_form=False, return_residuals=False):
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

    #bempp.api.log("Starting GMRES iteration")

    start_time = time.time()
    if solver == 'scipy':
        x, info = scipy.sparse.linalg.gmres(A_op, b_vec,
                                            tol=tol, restart=restart, maxiter=maxiter, callback=callback)
    elif solver == 'pyamg':
        from gmres import gmres
        x, info = gmres(A_op, b_vec, H=H,
                                            tol=tol, restrt=restart, maxiter=maxiter, callback=callback)
    end_time = time.time()
    #bempp.api.log("GMRES finished in {0} iterations and took {1:.2E} sec.".format(
    #    callback.count, end_time - start_time))

    if isinstance(A, BoundaryOperator) and isinstance(b, GridFunction):
        res_fun = GridFunction(A.domain, coefficients=x.ravel())

    else:
        res_fun = x
    if return_residuals:
        return res_fun, info, callback.residuals
    else:
        return res_fun, info

def relative_error(self, fun, element=None):
    """Compute the relative L^2 error compared to a given analytic function."""
    from bempp.api.integration import gauss_triangle_points_and_weights
    import numpy as np
    p_normal = np.array([[0]])
    global_diff = 0
    fun_l2_norm = 0
    accuracy_order = self.parameters.quadrature.far.single_order
    points, weights = gauss_triangle_points_and_weights(accuracy_order)
    npoints = points.shape[1]
    element_list = [element] if element is not None else list(self.grid.leaf_view.entity_iterator(0))
    for element in element_list:
        integration_elements = element.geometry.integration_elements(points)
        normal = element.geometry.normals(p_normal).reshape(3)
        global_dofs = element.geometry.local2global(points)
        fun_vals = np.zeros((self.component_count, npoints), dtype=self.dtype)
        for j in range(npoints):
            fun_vals[:, j] = fun(global_dofs[:, j], normal)
        diff = np.sum(np.abs(self.evaluate(element, points) - fun_vals)**2, axis=0)
        global_diff += np.sum(diff * integration_elements * weights)
        abs_fun_squared = np.sum(np.abs(fun_vals)**2, axis=0)
        fun_l2_norm += np.sum(abs_fun_squared * integration_elements * weights)
    return np.sqrt(global_diff/fun_l2_norm)



def translate(grid, a=0, b=0, c=0):
    vertices = grid.leaf_view.vertices.T
    vertices[:, 0] += a
    vertices[:, 1] += b
    vertices[:, 2] += c
    elements = grid.leaf_view.elements
    return grid_from_element_data(vertices.T,elements)


def plot_amplitude(ref):
    grid = ref.grid
    space = ref.space 
    P1 = bem.function_space(grid, "P", 1)
    N_P1 = int(P1.global_dof_count)
    elements0 = list(grid.leaf_view.entity_iterator(0))

    c = np.zeros(N_P1)
    for el0 in elements0:
        val1,val2,val3 = ref.evaluate(el0, np.array([[0,0],[1,0],[0,1]]).T).T
        val1 = np.abs(np.vdot(val1,val1))
        val2 = np.abs(np.vdot(val2,val2))
        val3 = np.abs(np.vdot(val3,val3))
        dof = P1.get_global_dofs(el0)

        c[dof[0]] = val1
        c[dof[1]] = val2
        c[dof[2]] = val3

    u_P1 = bem.GridFunction(P1, coefficients=c)
    return u_P1



def calderon_electric_field(grid, wave_number, parameters=None):
    """Return a pair (E^2, E) of the squared EFIE operator E^2 and E itself"""
    import bempp.api

    class EfieSquared(bempp.api.assembly.BoundaryOperator):
        """Implementation of the squared electric field operator."""

        def __init__(self, grid, wave_number, parameters):
            from bempp.api.assembly import InverseSparseDiscreteBoundaryOperator
            from bempp.api.space import project_operator

            bc_space = bempp.api.function_space(grid, "BC", 0)
            rbc_space = bempp.api.function_space(grid, "RBC", 0)
            rwg_space = bempp.api.function_space(grid, "B-RWG", 0)
            snc_space = bempp.api.function_space(grid, "B-SNC", 0)
            rwg_bary_space = bempp.api.function_space(
                grid.barycentric_grid(), "RWG", 0)
            snc_bary_space = bempp.api.function_space(
                grid.barycentric_grid(), "SNC", 0)
            super(EfieSquared, self).__init__(
                rwg_space, rwg_space, rbc_space, label="EFIE_SQUARED")

            self._efie_fine = electric_field(
                rwg_bary_space, rwg_bary_space, snc_bary_space, wave_number,
                parameters=parameters)
            self._efie = project_operator(
                self._efie_fine, domain=rwg_space,
                range_=rwg_space, dual_to_range=snc_space)
            self._efie2 = project_operator(
                self._efie_fine, domain=bc_space,
                range_=rwg_space, dual_to_range=rbc_space)
            self._ident = bempp.api.operators.boundary.sparse.identity(
                bc_space, rwg_space, snc_space)
            self._inv_ident = InverseSparseDiscreteBoundaryOperator(
                self._ident.weak_form())

        def _weak_form_impl(self):

            efie_weak = self._efie.weak_form()
            efie2_weak = self._efie2.weak_form()

            return efie2_weak * self._inv_ident * efie_weak, efie2_weak

    operator = EfieSquared(grid, wave_number, parameters)
    #pylint: disable=protected-access
    return operator, operator._efie2