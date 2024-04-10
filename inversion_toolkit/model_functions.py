import numpy as np
import SimPEG.electromagnetics.frequency_domain as fdem
from SimPEG.utils import mkvc
from SimPEG import (
    maps,
)


def get_1d_model(x_position, m_true, inversion_mesh_2d):
    """
    Function for constructring 1D model at a specified x-position from a 2D model

    Parameters
    ----------

    x_position: int, float
        specified x-position on the 2D model

    m_true: np.ndarray
        2D model

    inversion_mesh_2d: type(inversion_mesh_2d)
        2D mesh that the model which represents the model

    Returns
    -------
    layer_thicknesses: np.ndarray 
        array containing thicknesses of each layer
    
    log_resistivity_map: type(log_resistivity_map)

    log_resistivity_model: np.ndarray
        array containing log(resistivity) values at the specified x-position

    z_values_model: np.ndarray
        array containing the model values at every z-position at the specified x-position
    """
    # Find index of the x-position and the indices of grid cell centers that have the specified x-position
    x_index = np.where(inversion_mesh_2d.gridCC[:, 0]>= x_position)[0][0]
    x_indices = np.where(np.isclose(inversion_mesh_2d.gridCC[:, 0], inversion_mesh_2d.gridCC[:, 0][x_index]))

    # Get the z-values corresponding to the specified x-position
    z_values_mesh = inversion_mesh_2d.gridCC[x_indices, 1][0]
    z_values_mesh = z_values_mesh[::-1] + (inversion_mesh_2d.h[1][::-1])/2
    z_values_mesh = np.round(z_values_mesh)

    z_values_model = m_true[x_indices][::-1]

    unique_values, indices = np.unique(z_values_model, return_index=True)

    # Define 1D model from the 3D model
    log_resistivity_model = unique_values[np.argsort(indices)]

    layers = []
    layer_indices = np.sort(indices)
    for i in range(len(layer_indices)-1):
        idx = layer_indices[i+1]
        thickness = abs(z_values_mesh[idx])
        if i > 0:
            thickness = thickness - sum(layers[:i])
        layers.append(thickness)

    layer_thicknesses = np.array(layers)

    log_resistivity_map = maps.ExpMap(nP=len(layer_thicknesses)+1)

    return layer_thicknesses, log_resistivity_map, log_resistivity_model, z_values_model

def generate_survey(x_position, frequencies, moment, sep=10):

    # Defining transmitter locations
    xtx, ytx, ztx = np.meshgrid(x_position, [0], [30])
    source_locations = np.c_[mkvc(xtx), mkvc(ytx), mkvc(ztx)]
    ntx = np.size(xtx)

    # Define receiver locations
    xrx, yrx, zrx = np.meshgrid(x_position + sep, [0], [30])
    receiver_locations = np.c_[mkvc(xrx), mkvc(yrx), mkvc(zrx)]

    source_list = []  # Create empty list to store sources

    # Each unique location and frequency defines a new transmitter
    for ii in range(ntx):
        # Define receivers of different type at each location
        bzr_receiver = fdem.receivers.PointMagneticFluxDensitySecondary(
            receiver_locations[ii, :], "z", "real"
        )

        bzi_receiver = fdem.receivers.PointMagneticFluxDensitySecondary(
            receiver_locations[ii, :], "z", "imag"
        )
        receivers_list = [bzr_receiver, bzi_receiver]
        for jj in range(len(frequencies)):
            # Must define the transmitter properties and associated receivers
            source_list.append(
                fdem.sources.MagDipole(
                    receivers_list,
                    frequencies[jj],
                    source_locations[ii, :],
                    orientation="z",
                    momemt = moment
                )
            )

    # create the survey and problem objects for running the forward simulation
    survey = fdem.Survey(source_list)

    return survey