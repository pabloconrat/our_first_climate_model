#!/bin/bash -l

module load netcdf/4.7.0-gcc

g++ -o main main.cpp ./repwvl_V1.0_cpp/repwvl_thermal.cpp -lnetcdf_c++4 -L/opt/local/lib