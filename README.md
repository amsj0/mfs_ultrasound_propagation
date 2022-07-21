# Meshless boundary collocation techniques applied to Ultrasound Field Modelling

This code solves for ultrasound field modelling (Helmholtz PDE) using among others the idea behind the Distributed Point Source Method (DPSM). 
Solution for incident wave on homogeneous domain is sought by using Huygens principle based on real sources.
The boundary behaviour is sought by solving the Boundary Value Problem (BVP) using the Method of Fundamental Solution (MFS) based on virtual sources.

## Infinite Plane interface

Current development implements truncated solution to the infinite plane interface between two fluids.



## Octave code

Requires GNU Octave to run.

### Running

Setup **inputfile.txt** with model parameter ranges and run 

````
octave --no-gui run_from_inputfile.m
````

## Python code

Update the requirements with

````
pip install -U -r requirements.txt
````

### Running

For running the MFS code: 

````
python main.py inputfile.yaml
````

or simply

````
python mfsolution.py inputfile.yaml
````

For running the MFS and analysis code: 

````
python main.py inputfile.yaml output_path
````

For running the analysis code: 

````
python analyse.py inputfile.yaml /path/to/output
````

For running the MFS, analysis and tseries code : 

````
python main.py inputfile.yaml /path/to/output config.yaml
````

For running the tseries code: 

````
python tseries.py /path/to/output config.yaml
````
