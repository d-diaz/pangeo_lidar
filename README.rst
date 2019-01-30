=============================
pangeo_lidar
=============================

Processing lidar point clouds with a dask pipeline using Pangeo.

The lidar processing components in this toolkit currently include thin wrappers
around executables available in the FUSION_ and LAStools_ software packages,
which are executed using the Python subprocess module. FUSION and LAStools are
designed for use on Windows, so this Pangeo-Binder has been set up using a
Dockerfile which will install these software packages as well as wine on the
Linux server running these notebooks.

Try the notebooks on pangeo.binder.io_ : |Binder|

See http://pangeo.io for more information.

Features
--------

* TODO

.. _FUSION: http://forsys.cfr.washington.edu/fusion/fusionlatest.html

.. _LAStools: (https://rapidlasso.com/lastools/

.. _pangeo.binder.io: http://binder.pangeo.io/

.. |Binder| image:: http://binder.pangeo.io/badge.svg
    :target: http://binder.pangeo.io/v2/gh/d-diaz/pangeo_lidar/master
