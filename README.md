# Exporting data

## Overview
ION is a project that includes various C++ files and scripts for simulating and analyzing materials using the LAMMPS molecular dynamics software. The project provides functionality to export simulation data and manage simulations.

## Features
- **C++ Codebase**: Core functionality written in C++.
- **LAMMPS Export**: Tools to export data for use in LAMMPS simulations.
- **Simulation Management**: Scripts and configurations to manage simulation parameters.

## Usage
- Configure your simulations using `lammps_simulation.table`.
- Use the provided C++ executables to export and manage simulation data.
- Run LAMMPS simulations with the exported data.

## Repository Structure
- **.idea/**: Project files for IDE configurations.
- **CMakeLists.txt**: CMake build configuration.
- **QERFC_Li_Li/**: Directory containing additional resources or scripts.
- **export_lammps.cpp**: C++ source file for exporting LAMMPS data.
- **export_lammps.h**: Header file for exporting LAMMPS data.
- **lammps_simulation.table**: Table file for LAMMPS simulation configuration.
- **main.cpp**: Main C++ source file for the project.
