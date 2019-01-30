[![Build Status](https://travis-ci.org/Ecotrust/pyFIRS.svg?branch=master)](https://travis-ci.org/Ecotrust/pyFIRS)
# pyFIRS
pyFIRS: A Python-based Forest Information Retrieval System

This toolkit is intended to integrate several building blocks for processing point clouds from airborne laser scanning, remotely-sensed imagery, and ground-based forest inventory measurements to generate forest type maps including automated delineation of stands and imputation of stand- and plot-level attributes.  

It is under active development.

This toolkit includes a series of functions to generate raster and vector layers useful for forest management planning. It supports processing of raw point cloud data in LAS/LAZ format into geospatial data layers of forest canopy cover, height, etc. Routines for forest type classification (in terms of dominant species, size class, and stocking level) and the generation of tree lists and plot-level attributes to enable integration with growth-and-yield models are under development.

The lidar processing components in this toolkit currently include thin wrappers around executables available in the [FUSION](http://forsys.cfr.washington.edu/fusion/fusionlatest.html) and [LAStools](https://rapidlasso.com/lastools/) software packages, which are executed using the Python subprocess module. FUSION and LAStools are designed for use on Windows. 
