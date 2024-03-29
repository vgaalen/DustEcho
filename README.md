# Dust Reprocessing Echo Analiser Module for Tidal disruption events

This radiative transfer simulation is designed to model Dust Echos of Tidal Disruption Events

This work builds upon SKIRT developed by Astronomical Observatory, Ghent University; see https://skirt.ugent.be/ or https://github.com/SKIRT/SKIRT9.

## Installation Instructions

This software has only been tested on Unix based systems
Installation on another OS is possible, but not supported by these instructions

 - A recent C++ compiler with full support for the C++14 standard:
        On macOS, use Apple Clang v10 (included in Xcode v10) or later.
        On other Unix systems, use GNU g++ v5.4.0 or later or Intel icc v19.0.2 or later.
 - CMake v3.2.2 or later.
 - git
 - python3
 - the following python-packages (a venv is included in the repo)
    - astropy
    - scipy
    - matplotlib
    - datetime
    - tqdm
    - reportlab
    - lxml
    - jupyter notebook / jupyterlab

#### Get the source code

Clone the github repository, build SKIRT's C++ code, and retrieve SKIRT's resources (these contain files too large for github and are therefore fetched seperately)
```
git clone https://github.com/vgaalen/DustEcho [YOURDIR]
cd [YOURDIR]/SKIRT/git
chmod +rx makeSKIRT.sh
chmod +rx downloadResources.sh
./makeSKIRT.sh
./downloadResources.sh
```

The downloadResources.sh executable prompts you to download different resource packages, only SKIRT9_Resources_Core is required.

#### Trying out the simulation

This repository includes two notebooks that show the capabilities of this simulation.
DustEcho.ipynb shows how to generate a simulated lightcurve for ASASSN-15lh, while DustEchoCompare.ipynb does this using multiple values for the minimum grain size.
