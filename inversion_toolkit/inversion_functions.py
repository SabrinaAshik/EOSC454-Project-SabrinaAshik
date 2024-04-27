import numpy as np
from discretize import TensorMesh
from SimPEG import (
    data_misfit,
    inverse_problem,
    regularization,
    optimization,
    directives,
    inversion,
)

def define_inverse_problem_1DLayered(simulation_L2, data_object, regularization_mesh, reference_resistivity_model):
    # Define the data misfit. Here the data misfit is the L2 norm of the weighted
    # residual between the observed data and the data predicted for a given model.
    # The weighting is defined by the reciprocal of the uncertainties.
    dmis_L2 = data_misfit.L2DataMisfit(simulation=simulation_L2, data=data_object)
    reg_L2 = regularization.WeightedLeastSquares(
    regularization_mesh,
    length_scale_x=10.0,
    reference_model=reference_resistivity_model,
    reference_model_in_smooth=False,
    )
    opt_L2 = optimization.InexactGaussNewton(
    maxIter=100, maxIterLS=20, maxIterCG=20, tolCG=1e-3
    )
    inv_prob_L2 = inverse_problem.BaseInvProblem(dmis_L2, reg_L2, opt_L2)
    
    update_jacobi = directives.UpdatePreconditioner(update_every_iteration=True)
    starting_beta = directives.BetaEstimate_ByEig(beta0_ratio=1e1)
    beta_schedule = directives.BetaSchedule(coolingFactor=2.0, coolingRate=3)
    target_misfit = directives.TargetMisfit(chifact=1.0)

    directives_list_L2 = [update_jacobi, 
                        starting_beta, 
                        beta_schedule, 
                        target_misfit]
    
    inv_L2 = inversion.BaseInversion(inv_prob_L2, directives_list_L2)

    return inv_L2

def define_inverse_problem_1DLayered_irls(simulation_irls, data_object, regularization_mesh, reference_resistivity_model):
    dmis_irls = data_misfit.L2DataMisfit(simulation=simulation_irls, data=data_object)
    reg_irls = regularization.Sparse(
    regularization_mesh,
    alpha_s=1,
    alpha_x=1,
    reference_model=reference_resistivity_model,
    reference_model_in_smooth=True,
    norms=[1.0, 1.0],
    )
    opt_irls = optimization.InexactGaussNewton(
    maxIter=100, maxIterLS=20, maxIterCG=30, tolCG=1e-3
    )
    inv_prob_irls = inverse_problem.BaseInvProblem(dmis_irls, reg_irls, opt_irls)
    
    starting_beta_irls = directives.BetaEstimate_ByEig(beta0_ratio=1e-1)
    update_jacobi_irls = directives.UpdatePreconditioner(update_every_iteration=True)
    update_irls = directives.Update_IRLS(
        coolingFactor=2,
        coolingRate=3,
        f_min_change=1e-4,
        max_irls_iterations=30,
        chifact_start=1.0,
    )

    directives_list_irls = [update_irls, starting_beta_irls, update_jacobi_irls]
    
    inv_irls = inversion.BaseInversion(inv_prob_irls, directives_list_irls)

    return inv_irls

def get_regularization_mesh_and_models(estimated_resistivity):
    depth_min = 50  # top layer thickness
    depth_max = 5000.0  # depth to lowest layer
    geometric_factor = 1.2  # rate of thickness increase
    # Increase subsequent layer thicknesses by the geometric factors until
    # it reaches the maximum layer depth.
    layer_thicknesses = [depth_min]
    while np.sum(layer_thicknesses) < depth_max:
        layer_thicknesses.append(geometric_factor * layer_thicknesses[-1])

    n_layers = len(layer_thicknesses) + 1  # Number of layers

    # Define 1D cell widths
    h = np.r_[layer_thicknesses, layer_thicknesses[-1]]
    h = np.flipud(h)

    # Create regularization mesh
    regularization_mesh = TensorMesh([h], "N")

    starting_model = np.log((estimated_resistivity) * np.ones(n_layers))
    
    return regularization_mesh, layer_thicknesses, starting_model

def get_regularization_mesh_and_models_irls(estimated_resistivity):
    depth_min = 50  # top layer thickness
    depth_max = 5000.0  # depth to lowest layer
    geometric_factor = 1.5  # rate of thickness increase
    # Increase subsequent layer thicknesses by the geometric factors until
    # it reaches the maximum layer depth.
    layer_thicknesses_irls = [depth_min]
    while np.sum(layer_thicknesses_irls) < depth_max:
        layer_thicknesses_irls.append(geometric_factor * layer_thicknesses_irls[-1])

    n_layers = len(layer_thicknesses_irls) + 1  # Number of layers

    # Define 1D cell widths
    h = np.r_[layer_thicknesses_irls, layer_thicknesses_irls[-1]]
    h = np.flipud(h)

    # Create regularization mesh
    regularization_mesh_irls = TensorMesh([h], "N")

    starting_model = np.log((estimated_resistivity) * np.ones(n_layers))
    
    return regularization_mesh_irls, layer_thicknesses_irls, starting_model

def get_dobs_from_1d_data(imported_1d_data, location_index):
    dobs_real = imported_1d_data[:, 1::2].T
    dobs_imag = imported_1d_data[:, 2::2].T
    bz_real = dobs_real[location_index]
    bz_imag = dobs_imag[location_index]

    dobs = np.empty(bz_real.size + bz_imag.size, dtype=bz_real.dtype)
    dobs[0::2] = bz_real
    dobs[1::2] = bz_imag
    return dobs