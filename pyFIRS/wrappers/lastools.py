import os
import subprocess
import platform
import shutil
from pyFIRS.utils import listlike, PipelineError
import urllib.request
import geopandas as gpd
import numpy as np


# helper function for formatting command line arguments
def format_lastools_kws(**kwargs):
    '''Formats keyword arguments for LAStools command line usage.'''
    kws = []
    for key, value in kwargs.items():
        if isinstance(value, bool):
            kws.append('-{}'.format(key))
        elif listlike(value):
            kws.append('-{}'.format(key))
            for arg in value:
                kws.append(str(arg))
        else:
            kws.append('-{}'.format(key))
            kws.append(str(value))
    return kws


class LAStools_base(object):
    "A class for executing LAStools functions as methods"

    def __init__(self, src='C:\\lastools\\bin'):
        "Initialize with a path to the LAStools executables"
        self.src = src
        self.system = platform.system()

        # retrieve the documentation for each LAStool from the web
        tools = [
            self.lasview.__func__, self.lasinfo.__func__,
            self.lasground.__func__, self.lasclassify.__func__,
            self.las2dem.__func__, self.las2iso.__func__,
            self.lascolor.__func__, self.lasgrid.__func__,
            self.lasoverlap.__func__, self.lasoverage.__func__,
            self.lasboundary.__func__, self.lasclip.__func__,
            self.lasheight.__func__, self.lastrack.__func__,
            self.lascanopy.__func__, self.lasthin.__func__,
            self.lassort.__func__, self.lasduplicate.__func__,
            self.lascontrol.__func__, self.lastile.__func__,
            self.lassplit.__func__, self.txt2las.__func__,
            self.las2dem.__func__, self.las2iso.__func__,
            self.las2las.__func__, self.las2shp.__func__,
            self.las2tin.__func__, self.lasvoxel.__func__,
            self.lasreturn.__func__, self.laszip.__func__,
            self.lasindex.__func__, self.lasvalidate.__func__
        ]

        for tool in tools:
            name = tool.__name__
            URL = "http://www.cs.unc.edu/~isenburg/laszip/download/{}_README.txt".format(
                name)
            try:
                docstring = urllib.request.urlopen(URL).read().decode(
                    'utf-8', 'ignore')
            except:
                print('Error retrieving docstring for {}'.format(name))
                docstring = URL
            setattr(tool, '__doc__', docstring)

    def run(self, cmd, **kwargs):
        """Executes a LAStools command line tool.

        Formats kwargs provided in Pythonic format into the format expected
        by LAStools command line tools and executes the command line tool using
        the subprocess module.

        Parameters
        ----------
        cmd: string
            name of LAStools command line tool
        input: string, path to file(s)
            path to files to process with command line tool

        Returns
        -------
        CompletedProcess, a class from the subprocess module that includes
        attributes such as args, stdout, stderr, and returncode.
        """
        # check to see if echo was requested
        if 'echo' in kwargs:
            echo = kwargs['echo']
            del kwargs['echo']
        else:
            echo = False

        # check to see if output directory exists, if not, make it
        if 'odir' in kwargs:
            path = kwargs['odir']
            # makedirs will create whole directory tree recursively if needed
            os.makedirs(path, exist_ok=True)

        if 'wine_prefix' in kwargs:
            wine_prefix = kwargs['wine_prefix']
            del kwargs['wine_prefix']
        else:
            wine_prefix = None

        # format the kwargs
        kws = format_lastools_kws(**kwargs)

        # format the command to include the path to executables
        cmd = os.path.join(self.src, cmd)

        if self.system == 'Linux':
            # if we're on a linux system, execute the commands using WINE
            if wine_prefix:  # if we're using specific WINE server
                proc = subprocess.run(
                    'WINEPREFIX={} wine {}.exe {}'.format(
                        wine_prefix, cmd, ' '.join(kws)),
                    stderr=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    shell=True)
            else:  # no wine_prefix defined
                try:
                    proc = subprocess.run(['wine', cmd + '.exe', *kws],
                                          stderr=subprocess.PIPE,
                                          stdout=subprocess.PIPE)
                except OSError:  # we're probably running Windows Subsystem for Linux
                    # or don't have wine installed
                    proc = subprocess.run([cmd + '.exe', *kws],
                                          stderr=subprocess.PIPE,
                                          stdout=subprocess.PIPE)

        else:  # we're not on a Linux machine, use windows executable directly
            proc = subprocess.run([cmd + '.exe', *kws],
                                  stderr=subprocess.PIPE,
                                  stdout=subprocess.PIPE)
        if echo:
            print(proc.stdout.decode())
            print(proc.stderr.decode())

        if proc.returncode != 0:
            cmd_name = os.path.basename(cmd)
            error_msg = proc.stderr.decode().split('\r')[0]
            raise PipelineError(
                '''{} failed on "{}" with the following error message
                {}'''.format(cmd_name, kwargs['i'], error_msg))

        return proc

    def lasview(self, **kwargs):
        cmd = 'lasview'
        return self.run(cmd, **kwargs)

    def lasinfo(self, **kwargs):
        cmd = 'lasinfo'
        return self.run(cmd, **kwargs)

    def lasground(self, **kwargs):
        cmd = 'lasground'
        return self.run(cmd, **kwargs)

    def lasnoise(self, **kwargs):
        cmd = 'lasnoise'
        return self.run(cmd, **kwargs)

    def lasclassify(self, **kwargs):
        cmd = 'lasclassify'
        return self.run(cmd, **kwargs)

    def las2dem(self, **kwargs):
        cmd = 'las2dem'
        return self.run(cmd, **kwargs)

    def las2iso(self, **kwargs):
        cmd = 'las2iso'
        return self.run(cmd, **kwargs)

    def lascolor(self, **kwargs):
        cmd = 'lascolor'
        return self.run(cmd, **kwargs)

    def lasgrid(self, **kwargs):
        cmd = 'lasgrid'
        return self.run(cmd, **kwargs)

    def lasoverlap(self, **kwargs):
        cmd = 'lasoverlap'
        return self.run(cmd, **kwargs)

    def lasoverage(self, **kwargs):
        cmd = 'lasoverage'
        return self.run(cmd, **kwargs)

    def lasboundary(self, **kwargs):
        cmd = 'lasboundary'
        return self.run(cmd, **kwargs)

    def lasclip(self, **kwargs):
        cmd = 'lasclip'
        return self.run(cmd, **kwargs)

    def lasheight(self, **kwargs):
        cmd = 'lasheight'
        return self.run(cmd, **kwargs)

    def lastrack(self, **kwargs):
        cmd = 'lastrack'
        return self.run(cmd, **kwargs)

    def lascanopy(self, **kwargs):
        cmd = 'lascanopy'
        return self.run(cmd, **kwargs)

    def lasthin(self, **kwargs):
        cmd = 'lasthin'
        return self.run(cmd, **kwargs)

    def lassort(self, **kwargs):
        cmd = 'lassort'
        return self.run(cmd, **kwargs)

    def lasduplicate(self, **kwargs):
        cmd = 'lasduplicate'
        return self.run(cmd, **kwargs)

    def lascontrol(self, **kwargs):
        cmd = 'lascontrol'
        return self.run(cmd, **kwargs)

    def lastile(self, **kwargs):
        cmd = 'lastile'
        return self.run(cmd, **kwargs)

    def lassplit(self, **kwargs):
        cmd = 'lassplit'
        return self.run(cmd, **kwargs)

    def txt2las(self, **kwargs):
        cmd = 'txt2las'
        return self.run(cmd, **kwargs)

    def blast2dem(self, **kwargs):
        cmd = 'blast2dem'
        return self.run(cmd, **kwargs)

    def blast2iso(self, **kwargs):
        cmd = 'blast2iso'
        return self.run(cmd, **kwargs)

    def las2las(self, **kwargs):
        cmd = 'las2las'
        return self.run(cmd, **kwargs)

    def las2shp(self, **kwargs):
        cmd = 'las2shp'
        return self.run(cmd, **kwargs)

    def las2tin(self, **kwargs):
        cmd = 'las2shp'
        return self.run(cmd, **kwargs)

    def lasvoxel(self, **kwargs):
        cmd = 'lasvoxel'
        return self.run(cmd, **kwargs)

    def lasreturn(self, **kwargs):
        cmd = 'lasreturn'
        return self.run(cmd, **kwargs)

    def laszip(self, **kwargs):
        cmd = 'laszip'
        return self.run(cmd, **kwargs)

    def lasindex(self, **kwargs):
        cmd = 'lasindex'
        return self.run(cmd, **kwargs)

    def lasvalidate(self, **kwargs):
        cmd = 'lasvalidate'
        return self.run(cmd, **kwargs)


# Pythonic wrappers for LAStools command line tools
def get_bounds(lasinfo):
    '''Parses the minimum and maximum X, Y, and Z values from LASinfo output.

    This is a helper function used by the pitfree function.

    Parameters
    ----------
    lasinfo: string
        result produced by executing the lasinfo command line tool on a lidar
        data file

    Returns
    -------
    bounds: tuple
        a 6-tuple containing (xmin, ymin, zmin, xmax, ymax, zmax)
    '''
    min_start = lasinfo.index('min x y z:')
    min_stop = lasinfo.index('\r\n', min_start)
    min_line = lasinfo[min_start:min_stop]
    min_vals = min_line.split(':')[1].strip().split(' ')
    mins = (float(min_vals[0]), float(min_vals[1]), float(min_vals[2]))

    max_start = lasinfo.index('max x y z:', min_stop)
    max_stop = lasinfo.index('\r\n', max_start)
    max_line = lasinfo[max_start:max_stop]
    max_vals = max_line.split(':')[1].strip().split(' ')
    maxs = (float(max_vals[0]), float(max_vals[1]), float(max_vals[2]))
    return mins + maxs


class useLAStools(LAStools_base):
    """A class which inherits the command-line tools from LAStools_base and
    provides additional methods for processing lidar data that chain together
    these lower-level commands and which may integrate some Python processing.
    """

    def pitfree(self,
                lasfile,
                outdir,
                units,
                xy_res=None,
                z_res=None,
                splat_radius=None,
                max_TIN_edge=None,
                blast=False,
                cleanup=True,
                echo=False,
                wine_prefix=None):
        '''Creates a pit-free Canopy Height Model from a lidar point cloud.

        This function chains together several LAStools command line tools to
        produce a pit-free Canopy Height Model (CHM) from a raw lidar point
        cloud. A working subdirectory is created in the same folder where the
        input lidar data file is located to hold intermediate files from the
        process, which will be deleted (by default) upon completion of the CHM.
        A GeoTiff of the CHM is output to the same folder as the inputfile named
        {inputfile}_chm_pitfree.tif

        This method was first described by:

            Khosravipour, A. et al. (2014) "Generating Pit-free Canopy Height
            Models from Airborne Lidar."" Photogrammetric Engineering & Remote
            Sensing 80(9): 863–872.

        The method was elaborated by co-author/LAStools developer Martin
        Isenburg in a blog post, "Rasterizing Perfect Canopy Height Models from
        LiDAR". The default values used here are drawn from this blog post:
        https://rapidlasso.com/2014/11/04/rasterizing-perfect-canopy-height-models-from-lidar/

        The method employed here varies from the blog post by using blast2dem to
        allow out-of-core processing of larger point clouds to avoid running out
        of memory when calculating intermediate DEMs.

        Parameters
        ----------
        lasfile: string, path to file (required)
            Lidar point cloud input file to process.
        outdir: string, path to directory (required)
            Output directory where pit free CHM will be saved
        units: string (required)
            'm' for meters or 'ft' for feet
        xy_res: numeric (optional)
            Size of grid cells for Canopy Height Model, in same units as lidar
            data. Used in the `step` argument of las2dem. Default is 0.33333 if
            units are in meters or 1.0 if units are in feet.
        z_res: numeric (optional)
            Height of vertical slices used to build CHM layers. Will always use
            layers from 0-2m and 2-5m, then will stack on layers z_res thick.
        splat_radius: numeric (optional)
            Distance by which lidar points will be replicated in each of eight
            directions, which used in the `subcircle` argument of lasthin.
            Default is 0.1 if units are in meters or 0.3 if units are in feet.
        max_TIN_edge: numeric (optional)
            Maximum length of edges for points to remain connected in TIN
            created by las2dem. Used in the `kill` argument of las2dem. Always
            considered by las2dem in units of meters, and conversion to feet
            is handled by las2dem if the LAS header indicates units are in ft.
            Default is 1.0 meters.
        blast: boolean (optional)
            Whether or not to use BLAST commands from LAStools to handle larger
            files. If True, will employ blast2dem and rather than las2dem.
            Defaults to False (i.e., to use las2dem).
        cleanup: boolean (optional)
            Whether or not to remove the temporary working directory and
            intermediate files produced. Defaults to True.
        echo: boolean (optional)
            If true, will echo to stdout and stderr the calls of all the
            LAStools.
        wine_prefix: integer or string (optional)
            If provided when run on a Linux OS, identifies a specific WINE
            server to use for executing the command. Defaults to None.
        '''
        path_to_file = os.path.abspath(lasfile)
        path, fname = os.path.split(path_to_file)
        basename = fname.split('.')[0]

        # make a temporary working directory
        tmpdir = os.path.join(outdir, 'work_{}'.format(basename))
        os.makedirs(tmpdir, exist_ok=True)

        # run lasheight to normalize point cloud
        odir = os.path.join(tmpdir, 'normalized')
        proc_height = self.lasheight(
            i=path_to_file,
            odir=odir,
            olaz=True,
            replace_z=True,
            keep_class=(1, 2, 5),
            drop_below=-0.1,  # drop points below the ground
            echo=echo,
            wine_prefix=wine_prefix)

        # get the minimum and maximum normalized heights
        # we'll use these later for creating layered canopy height models
        infile = os.path.join(tmpdir, 'normalized', '*.laz')
        info_proc = self.lasinfo(i=infile, wine_prefix=wine_prefix)
        lasinfo = info_proc.stderr.decode()
        _, _, zmin, _, _, zmax = get_bounds(lasinfo)

        # check to see if we need to use defaults
        if units.lower() in ('m', 'meter', 'meters'):
            if not xy_res:
                xy_res = 0.33333
            if not z_res:
                z_res = 5.0
            if not splat_radius:
                splat_radius = 0.1
            if not max_TIN_edge:
                max_TIN_edge = 1.0
            hts = [0.0, 2.0] + [x for x in np.arange(5.0, zmax, z_res)]
        elif units.lower() in ('f', 'ft', 'feet'):
            if not xy_res:
                xy_res = 1.0
            if not z_res:
                z_res = 15.0
            if not splat_radius:
                splat_radius = 0.3
            if not max_TIN_edge:
                max_TIN_edge = 1.0  # blast2dem converts from meters to feet
                # so we use the same value for meters or feet
            hts = [0.0, 6.56168
                   ] + [x for x in np.arange(16.4042, zmax, z_res).tolist()]
        else:
            raise ValueError('{} is not recognized units'.format(units))

        # create DEM of ground for minimum value of pitfree CHM
        infile = os.path.join(tmpdir, 'normalized', '*.laz')  # *_height.laz
        odir = os.path.join(tmpdir, 'chm_layers')
        odix = '_chm_ground'
        if blast:
            proc_dem1 = self.blast2dem(
                i=infile,
                odir=odir,
                odix=odix,
                obil=True,
                drop_z_above=0.1,
                step=xy_res,  # resolution of ground model
                use_tile_bb=True,  # trim the tile buffers
                echo=echo,
                wine_prefix=wine_prefix)
        else:
            proc_dem1 = self.las2dem(
                i=infile,
                odir=odir,
                odix=odix,
                obil=True,
                drop_z_above=0.1,
                step=xy_res,  # resolution of ground model
                use_tile_bb=True,  # trim the tile buffers
                echo=echo,
                wine_prefix=wine_prefix)
        # "splat" and thin the lidar point cloud to get highest points using a
        # finer resolution than our final CHM will be
        infile = os.path.join(tmpdir, 'normalized', '*.laz')  # *_height.laz
        odir = os.path.join(tmpdir, 'splatted')
        proc_thin = self.lasthin(
            i=infile,
            odir=odir,
            olaz=True,
            highest=True,
            subcircle=splat_radius,
            step=xy_res / 2.0,
            echo=echo,
            wine_prefix=wine_prefix)

        # using the "splatted" lidar point cloud, generate CHM layers above
        # ground, above 2m, and then in 5m increments up to zmax...
        # las2dem first makes a TIN and then rasterizes to grid
        infile = os.path.join(tmpdir, 'splatted', '*.laz')

        # loop through the layers
        dem2_procs = []
        for i, ht in enumerate(hts):
            odix = '_chm_{:02d}_{:03d}'.format(i, int(ht))
            odir = os.path.join(tmpdir, 'chm_layers')
            if blast:
                proc_dem2 = self.blast2dem(
                    i=infile,
                    odir=odir,
                    odix=odix,
                    obil=True,
                    drop_z_below=ht,  # specify layer height from ground
                    kill=max_TIN_edge,  # trim edges in TIN > max_TIN_edge
                    step=xy_res,  # resolution of layer DEM
                    use_tile_bb=True,  # trim tile buffer
                    echo=echo,
                    wine_prefix=wine_prefix)
            else:
                proc_dem2 = self.las2dem(
                    i=infile,
                    odir=odir,
                    odix=odix,
                    obil=True,
                    drop_z_below=ht,  # specify layer height from ground
                    kill=max_TIN_edge,  # trim edges in TIN > max_TIN_edge
                    step=xy_res,  # resolution of layer DEM
                    use_tile_bb=True,  # trim tile buffer
                    echo=echo,
                    wine_prefix=wine_prefix)
            dem2_procs.append(proc_dem2)

        # merge the CHM layers into a single pit free CHM raster
        infiles = os.path.join(tmpdir, 'chm_layers', '*.bil')
        outfile = basename + '_chm_pitfree.bil'
        proc_grid = self.lasgrid(
            i=infiles,
            merged=True,
            o=outfile,
            odir=outdir,
            highest=True,
            step=xy_res,  # resolution of pit-free CHM
            echo=echo,
            wine_prefix=wine_prefix)

        if cleanup:
            shutil.rmtree(tmpdir)

        return (proc_height, proc_dem1, proc_thin, dem2_procs, proc_grid)


# def clean_buffer_polys(poly_shp, tile_shp, odir, simp_tol=None, simp_topol=None):
#     """Removes polygons within the buffer zone of a tile.
#
#     This function removes polygons from a shapefile that fall in the buffered
#     area of point cloud tile. When building footprints or tree crowns (for
#     example) are delineated from a point cloud, a buffer around the tile is
#     generally be used to avoid edge effects. This tool computes the centroid of
#     each polygon and determines whether it falls within the bounds of the
#     unbuffered tile. It outputs a new shapefile containing only those polygons
#     whose centroids fall within the unbuffered tile.
#
#     The polygons may be simplified using optional arguments simp_tol and
#     simp_topol to reduce the number of points that define their boundaries.
#
#     Parameters
#     ----------
#     polygons_shp: string, path to shapefile (required)
#         A shapefile containing the polygons delineated within a buffered tile.
#     tile_shp: string, path to shapefile (required)
#         A shapefile containing the bounds of the tile WITHOUT buffers
#     odir: string, path to directory (required)
#         Path to the output directory for the new shapefile
#     simp_tol = numeric,
#         Tolerance level for simplification. All points within a simplified
#         geometry will be no more than simp_tol from the original.
#     simp_topol = boolean (optional)
#         Whether or not to preserve topology of polygons. If False, a quicker
#         algorithm will be used, but may produce self-intersecting or otherwise
#         invalid geometries.
#     """
#     fname = os.path.basename(poly_shp)
#     outfile = os.path.join(odir, fname)
#     os.makedirs(odir, exist_ok=True)
#
#     tile_boundary = gpd.read_file(tile_shp)
#     polys = gpd.read_file(poly_shp)
#
#     # boolean indicator of whether each polygon falls within tile boundary
#     clean_polys_ix = polys.centroid.within(tile_boundary.loc[0].geometry)
#     # retrieve the polygons within the boundary
#     clean_polys = polys[clean_polys_ix]
#
#     if simp_tol:
#         clean_polys = clean_polys.simplify(simp_tol, simp_topol)
#
#     if len(clean_polys) > 0:
#         clean_polys.to_file(outfile)
#
