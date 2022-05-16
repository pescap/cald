# -*- coding: utf-8 -*-
# copied from cald_energy/min_comparison_uniform-energy.py

from __future__ import division
import bempp.api as bem
import numpy as np
from login import get_h, gmres, relative_error, translate, calderon_electric_field,assemble_singular_part
from numpy.linalg import norm
import time
#from geometry import destroyer,reentrant_cube_exterior

########################################################################
########################################################################

####################################
# PARAMETERS 
####################################
# Penser a decaler
from bempp.core.common.global_parameters import global_parameters

#quad_mu= [4,3,2,6]
quad_mu= [1,1,1,2]
#quad_mu= [1,1,1,2]
quad_nu= [4,3,2,6]
quad_dense= [10,10,10,12]
quad_eps= [5,4,3,7]

distance = 5
mu_nf = 0

mu = 0.1
nu = 0.001
eps = 0.0001

# wavenumber
kappa = 10
# precision (points per wavelength)
precision = 10
h = 2.0 * np.pi / (precision * kappa)

bem.global_parameters.assembly.boundary_operator_assembly_type = "dense"

# To change
#grid = bem.import_grid("cube_k_5.msh")
#name = "grid_k_10_p_10.msh"
name = "cube_uniform.msh"
grid = bem.import_grid(name)



from bempp.api.space import project_operator

# Specify the incident field and its tangential component

def incident_field(x):
    return np.array([0.*x[0], 0.*x[1], -np.exp(1j * kappa * x[0])])
    #return np.array([np.exp(1j * kappa * x[2]), 0. * x[2], 0. * x[2]])

def tangential_trace(x, n, domain_index, result):
    result[:] = np.cross(incident_field(x), n, axis=0)

################################################
# Parameters nu

parameters_nu = global_parameters()
if nu == -1:
	parameters_nu.assembly.boundary_operator_assembly_type='dense'
else:
	parameters_nu.assembly.boundary_operator_assembly_type='hmat'
parameters_nu.hmat.eps = nu
parameters_nu.quadrature.double_singular = quad_nu[3]
parameters_nu.quadrature.far.double_order = quad_nu[2]
parameters_nu.quadrature.medium.double_order = quad_nu[1]
parameters_nu.quadrature.near.double_order = quad_nu[0]

################################################
# Parameters mu

parameters_mu = global_parameters()
if mu == -1:
	parameters_mu.assembly.boundary_operator_assembly_type='dense'
else:
	parameters_mu.assembly.boundary_operator_assembly_type='hmat'
parameters_mu.hmat.eps = mu
parameters_mu.quadrature.double_singular = quad_mu[3]
parameters_mu.quadrature.far.double_order = quad_mu[2]
parameters_mu.quadrature.medium.double_order = quad_mu[1]
parameters_mu.quadrature.near.double_order = quad_mu[0]
########################################################################
########################################################################



parameters_dense = global_parameters()
parameters_dense.assembly.boundary_operator_assembly_type='dense'
parameters_dense.hmat.eps = -1
parameters_dense.quadrature.double_singular = quad_dense[3]
parameters_dense.quadrature.far.double_order = quad_dense[2]
parameters_dense.quadrature.medium.double_order = quad_dense[1]
parameters_dense.quadrature.near.double_order = quad_dense[0]



parameters_eps = global_parameters()
if eps == -1:
	parameters_eps.assembly.boundary_operator_assembly_type='dense'
else:
	parameters_eps.assembly.boundary_operator_assembly_type='hmat'
parameters_eps.hmat.eps = eps
parameters_eps.quadrature.double_singular = quad_eps[3]
parameters_eps.quadrature.far.double_order = quad_eps[2]
parameters_eps.quadrature.medium.double_order = quad_eps[1]
parameters_eps.quadrature.near.double_order = quad_eps[0]





########################################################################
########################################################################

bc_space = bem.function_space(grid, "BC", 0)
rbc_space = bem.function_space(grid, "RBC", 0)
rwg_space = bem.function_space(grid, "RWG", 0)
snc_space = bem.function_space(grid, "SNC", 0)

rwg_bary_space = bem.function_space(grid.barycentric_grid(), "RWG", 0)
snc_bary_space = bem.function_space(grid.barycentric_grid(), "SNC", 0)

b_rwg_space = bem.function_space(grid, "B-RWG", 0)
b_snc_space = bem.function_space(grid, "B-SNC", 0)

N = int(rwg_space.global_dof_count)
#print N ,'NDOF'

efie_none_ref = bem.operators.boundary.maxwell.electric_field(rwg_space,rwg_space, snc_space, kappa, parameters=parameters_dense)
# define the operators for standard formulation
efie_none = bem.operators.boundary.maxwell.electric_field(rwg_space,rwg_space, snc_space, kappa, parameters = parameters_nu)
identity_none = bem.operators.boundary.sparse.identity(rwg_space,rwg_space, snc_space)
efie_bary = bem.operators.boundary.maxwell.electric_field(rwg_bary_space,rwg_bary_space, snc_bary_space, kappa, parameters = parameters_eps)


# define the operators for preconditioned formulation
efie_bc = bem.operators.boundary.maxwell.electric_field(bc_space,b_rwg_space, rbc_space, kappa, parameters=parameters_mu)
identity_hyp = bem.operators.boundary.sparse.identity(bc_space, b_rwg_space, b_snc_space)
identity_hyp_t = bem.operators.boundary.sparse.identity(b_rwg_space, b_rwg_space, rbc_space)
inv_identity_hyp = bem.assembly.InverseSparseDiscreteBoundaryOperator(identity_hyp.weak_form())
inv_identity_hyp_t = bem.assembly.InverseSparseDiscreteBoundaryOperator(identity_hyp_t.weak_form())

trace_fun_none = bem.GridFunction(rwg_space, fun=tangential_trace)
trace_fun_b_none = bem.GridFunction(b_rwg_space, fun=tangential_trace)
rhs_none = (identity_none * trace_fun_none).projections()
# Get the weak form for the standard formulation

t0 = time.time()
efie_none_wf = efie_none.weak_form()
ta_none = time.time()-t0

print('assembly of Cald prec bi')
t0 = time.time()
efie_bc_sf = efie_bc.strong_form()
tc_bi = time.time()-t0

efie_hyp_sf = efie_bc_sf * inv_identity_hyp * efie_none_wf
rhs_hyp_paul =  efie_bc_sf * inv_identity_hyp * rhs_none


tolerance = 1e-5
restart = 200
maxiter = 100000

print("solver of bi scipy")
t0 = time.time()
c_paul,  info_paul, res_paul  = gmres(efie_hyp_sf, rhs_hyp_paul, tol=tolerance, restart=restart, return_residuals=True, maxiter=maxiter)
ts_bi = time.time()-t0

