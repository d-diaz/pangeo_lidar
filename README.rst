=============================
pangeo_lidar
=============================

Processing lidar point clouds with a Dask distributed pipeline using Pangeo.

The lidar processing components in this toolkit are thin wrappers around
command line tools from the FUSION_ and LAStools_ software packages executed
using the Python subprocess module.

FUSION and LAStools are designed for Windows, so this Pangeo-Binder is set up
using a Dockerfile which installs these software packages as well as wine to
enable running these Windows binaries on Linux.

Try the notebooks on pangeo.binder.io_ : |Binder|

See http://pangeo.io for more information.

Features
--------

* TODO

.. _FUSION: http://forsys.cfr.washington.edu/fusion/fusionlatest.html

.. _LAStools: https://rapidlasso.com/lastools/

.. _pangeo.binder.io: http://binder.pangeo.io/

.. |Binder| image:: http://binder.pangeo.io/badge.svg
    :target: http://binder.pangeo.io/v2/gh/d-diaz/pangeo_lidar/master
