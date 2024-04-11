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
