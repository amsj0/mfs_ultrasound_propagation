# Meshless boundary collocation techniques applied to Ultrasound Field Modelling

This code solves for ultrasound field modelling (Helmholtz PDE) using among others the idea behind the Distributed Point Source Method (DPSM). 
Solution for incident wave on homogeneous domain is sought by using Huygens principle based on real sources.
The boundary behaviour is sought by solving the Boundary Value Proble (BVP) using the Method of Fundamental Solution (MFS) based on virtual sources.

## Infinite Plane interface

Current development implements truncated solution to the infinite plane interface between two fluids.

## Requirements

Requires GNU Octave to run.

### Running

Setup **inputfile.txt** with model parameter ranges and run 

````
octave --no-gui run_from_inputfile.m
````
