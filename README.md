# SU(N) scalar EMT measurement code

- License: GPLv2
## Author
- Joseph K. L. Lee joseph.lee@ed.ac.uk

## Introduction
This code computes the scalar Energy-Momentum Tensor two point functions with Wilson Flow using the Hadrons library (https://github.com/aportelli/Hadrons.git). 


## Install
Download, compile, and install the Hadrons library from https://github.com/aportelli/Hadrons.git .

This measurement code can be built using:
``` bash
./bootstrap.sh
mkdir build; cd build
../configure --with-hadrons=<dir>
make
```
`<dir>` is the installation prefix of Hadrons.

## Run
The measurement code can be run using:
```
scalar-emt-meas <param-xml>
```
`<param-xml>` is the .xml file containing the measurement run parameters. An example is found in `param-example.xml`.