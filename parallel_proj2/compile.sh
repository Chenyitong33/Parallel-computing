#!/usr/bin/env bash
g++ -o nbody_omp -fopenmp nbody_omp.cpp
mpic++ -fopenmp -o nbody nbody.cpp
