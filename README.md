# Optimization of VLSI chip area using Mixed Integer Linear Programming

NOTE: This repository uses the MOSEK library in C++ to solve the problem. Click [here](https://github.com/ShuvoNewaz/MILP-VLSI-Floorplanning-Python/) to view the repository in Python.

This repository does not yet have an argument parser like it's Python counterpart, but it has all the other features. This repository is for the Linux operating system and is based on C++11. The visualization uses the [matplotlib](https://matplotlib.org/) library with Python. To use this repository, please follow these steps:

- Clone this repository or download as a zip.
- This project makes use of an LP-solver named MOSEK. The tool can be downloaded and the license can be obtained from [MOSEK's website](https://www.mosek.com/resources/getting-started/). Once registered, follow the instructions regarding the directory setup for MOSEK in the email. In particular, pay attention to the directory where the license file is kept.
- The setup for the MOSEK libraries and header files are dependent on the operating system. For instance, g++ may not be used with Windows to run MOSEK. The [run_template.sh](run_template.sh) outlines the template paths for the required header files and libraries.
- The input arguments are very similar to the Python version. Edit the "main.cpp" file to run as required.

The models are created on the basis of the work by [Sutanthavibul et al](https://dl.acm.org/doi/abs/10.1145/123186.123255).