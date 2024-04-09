import numpy as np
from discretize import TensorMesh
    
import SimPEG.electromagnetics.frequency_domain as fdem
from SimPEG import (
    maps,
    data,
    data_misfit,
    inverse_problem,
    regularization,
    optimization,
    directives,
    inversion,
)

def recover_model(location_index, x_position):
    # path to the directory containing our data
    dir_path = "./3d_simulation_data/"
    data_filename = dir_path + "em3dfm_1d_data.txt"

    dobs = np.loadtxt(str(data_filename), skiprows=1)
    dobs_real = dobs[:, 1::2].T
    dobs_imag = dobs[:, 2::2].T
    bz_real = dobs_real[location_index]
    bz_imag = dobs_imag[location_index]

    dobs_recovered = np.empty(bz_real.size + bz_imag.size, dtype=bz_real.dtype)
    dobs_recovered[0::2] = bz_real
    dobs_recovered[1::2] = bz_imag

    # Define Survey
    frequencies = dobs[:, 0]

    source_location = np.array([x_position, 0.0, 30.0])
    receiver_location = np.array([x_position + 10, 0.0, 30.0])

    # Each unique location and frequency defines a new transmitter
    receiver_list = []
    receiver_list.append(
        fdem.receivers.PointMagneticFluxDensitySecondary(
            receiver_location,
            orientation='z',
            component="real",
        )
    )
    receiver_list.append(
        fdem.receivers.PointMagneticFluxDensitySecondary(
            receiver_location,
            orientation='z',
            component="imag",
        )
    )

    source_list = []
    for freq in frequencies:
        source_list.append(
            fdem.sources.MagDipole(
                receiver_list=receiver_list,
                frequency=freq,
                location=source_location,
                orientation="z",
                moment=1.0,
            )
        )
    survey = fdem.survey.Survey(source_list)

    uncertainties = 0.01 * np.abs(dobs_recovered) * np.ones(np.shape(dobs_recovered))
    data_object = data.Data(survey, dobs=dobs_recovered, noise_floor=uncertainties)

    # estimated host conductivity (S/m)
    estimated_resistivity = 100

    depth_min = 10  # top layer thickness
    depth_max = 6000.0  # depth to lowest layer
    geometric_factor = 1.1  # rate of thickness increase

    # Increase subsequent layer thicknesses by the geometric factors until
    # it reaches the maximum layer depth.
    layer_thicknesses = [depth_min]
    while np.sum(layer_thicknesses) < depth_max:
        layer_thicknesses.append(geometric_factor * layer_thicknesses[-1])

    n_layers = len(layer_thicknesses) + 1  # Number of layers

    log_resistivity_map = maps.ExpMap(nP=n_layers)

    # Starting model is log-conductivity values (S/m)
    starting_resistivity_model = np.log((estimated_resistivity) * np.ones(n_layers))

    # Reference model is also log-resistivity values (S/m)
    reference_resistivity_model = starting_resistivity_model.copy()

    simulation_L2 = fdem.Simulation1DLayered(
        survey=survey, thicknesses=layer_thicknesses, rhoMap=log_resistivity_map
    )

    dmis_L2 = data_misfit.L2DataMisfit(simulation=simulation_L2, data=data_object)

    # Define 1D cell widths
    h = np.r_[layer_thicknesses, layer_thicknesses[-1]]
    h = np.flipud(h)

    # Create regularization mesh
    regularization_mesh = TensorMesh([h], "N")

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
    starting_beta = directives.BetaEstimate_ByEig(beta0_ratio=5)
    beta_schedule = directives.BetaSchedule(coolingFactor=2.0, coolingRate=3)
    target_misfit = directives.TargetMisfit(chifact=1.0)

    directives_list_L2 = [update_jacobi, 
                        starting_beta, 
                        beta_schedule, 
                        target_misfit]

    # Here we combine the inverse problem and the set of directives
    inv_L2 = inversion.BaseInversion(inv_prob_L2, directives_list_L2)

    # Run the inversion
    recovered_model_L2 = inv_L2.run(starting_resistivity_model)

    return recovered_model_L2, layer_thicknesses