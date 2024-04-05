# EOSC 454 Final Project: Frequncy Domain Electromagnetic Data Inversion of Multi-Layered Model

## Description

This project focuses on conducting Frequency Domain Electromagnetic (FDEM) Data Inversion of the following  multi-Layered synthetic model with the resistivity values listed for Model A:

<img width="245" alt="image" src="https://github.com/SabrinaAshik/EOSC454-Project-SabrinaAshik/assets/70414346/aa39796a-020f-413b-8bba-9faf54f8a7c0">

(De Groot‐Hedlin, C., & Constable, S. (2004). Inversion of magnetotelluric data for 2D structure with sharp resistivity contrasts. Geophysics, 69(1), 78–86. https://doi.org/10.1190/1.1649377)

FDEM is a geophysical method widely used in hydrogeology, mineral exploration, and environmental studies. The inversion process aims to estimate the subsurface electrical conductivity distribution by fitting observed electromagnetic data to a forward model.

## Objectives
### Starting in 1D: What does a sounding look like on the far left vs far right end of the model?
- The notebooks in the 1d_files folder explores the 1D forward simulation and inversion for the soundings on the far left and far right.
### Moving on to 3D: Run a 3D forward simulation to generate data, and then invert each sounding independently in 1D.
- The notebooks in the 3d_files folder currently contains the 3D forward simulation file. The files in this folder is still work in progress.

## Installation

To run the project, follow these installation instructions:

1. Clone the repository:

    ```
    git clone https://github.com/SabrinaAshik/EOSC454-Project-SabrinaAshik.git
    ```

2. Navigate to the project directory:

    ```
    cd EOSC454-Project-SabrinaAshik
    ```

3. Create a conda environment using the provided environment.yml file:

    ```
    conda env create -f environment.yml
    ```

4. Activate the conda environment:

    ```
    conda activate eosc454
    ```

5. Run the project scripts or notebooks.

## License

This project is licensed under the terms of the MIT License. See the [LICENSE](LICENSE) file for details.

