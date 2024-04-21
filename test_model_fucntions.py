import numpy as np
import discretize
from inversion_toolkit.model_functions import get_1d_model, generate_survey

m_true = np.load("./models/2d-model-array.npy")
inversion_mesh_2d = discretize.load_mesh("./meshes/inversion_mesh_2d.json")

def test_left_edge_model():
    correct_layers = [500]
    correct_rhos = [1000, 100]

    layers, _, log_rho, _ = get_1d_model(0, m_true=m_true, inversion_mesh_2d=inversion_mesh_2d)

    print("correct layers: ", correct_layers)
    print("function returned layers: ", layers)

    assert np.allclose(layers, correct_layers)

    print("correct log(rho): ", np.log(correct_rhos))
    print("function returned log(rho) ", log_rho)

    assert np.allclose(log_rho, np.log(correct_rhos))

def test_right_edge_model():
    correct_layers = [1000, 1500]
    correct_rhos = [1000, 10, 100]

    layers, _, log_rho, _ = get_1d_model(10000, m_true=m_true, inversion_mesh_2d=inversion_mesh_2d)

    print("correct layers: ", correct_layers)
    print("function returned layers: ", layers)

    assert np.allclose(layers, correct_layers)

    print("correct log(rho): ", np.log(correct_rhos))
    print("function returned log(rho) ", log_rho)

    assert np.allclose(log_rho, np.log(correct_rhos))

if __name__ == "__main__":
    print("Testing the left edge model...")
    test_left_edge_model()
    print("Testing the right edge model...")
    test_right_edge_model()