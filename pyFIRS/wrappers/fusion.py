import os
import subprocess
import platform
import warnings
from pyFIRS.utils import listlike, PipelineError


# helper functions for formatting command line arguments
def format_fusion_kws(**kwargs):
    '''Formats keyword arguments for FUSION command line usage.'''
    kws = []
    for key, value in kwargs.items():
        # catch and replace kwarg names python doesn't like (class, ascii)
        if key == 'las_class':  # can't specify 'class' as a kwarg
            key = 'class'
        if key == 'asc':  # can't specify 'ascii' as a kwarg
            key = 'ascii'

        if isinstance(value, bool):
            kws.append('/{}'.format(key))
        elif listlike(value):
            kws.append('/{}:'.format(key) + ','.join(str(x) for x in value))
        else:
            kws.append('/{}:'.format(key) + str(value).replace('/', '\\'))
    return kws


def format_fusion_args(arg):
    '''Formats positional arguments for FUSION command line usage'''
    if listlike(arg):
        return " ".join(str(x) for x in arg)
    else:
        return str(arg).replace('/', '\\')


# Pythonic wrappers for FUSION command line tools
class useFUSION(object):
    "A class for executing FUSION functions as methods"

    def __init__(self, src='C:\\FUSION'):
        "Initialize with a path to the FUSION executables"
        self.src = src
        self.system = platform.system()

    def run(self, cmd, *params, **kwargs):
        "Formats and executes a FUSION command line call using subprocess"
        # prepend the path to FUSION tools to the user-specified command
        cmd = os.path.join(self.src, cmd)

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
            del kwargs['odir']

        if 'wine_prefix' in kwargs:
            wine_prefix = kwargs['wine_prefix']
            del kwargs['wine_prefix']
        else:
            wine_prefix = None

        # format kwargs as FUSION 'switches'
        switches = format_fusion_kws(**kwargs)

        # format the required parameters for each function as strings
        params = [format_fusion_args(param) for param in params]

        cmd = os.path.join(self.src, cmd)

        if self.system == 'Linux':
            # if we're on a linux system, execute the commands using WINE
            if wine_prefix:  # if we're using specific WINE server
                proc = subprocess.run(
                    'WINEPREFIX={} wine {}.exe {} {}'.format(
                        wine_prefix, cmd, ' '.join(switches),
                        ' '.join(params)),
                    stderr=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    shell=True)
            else:  # no wine_prefix defined
                try:
                    proc = subprocess.run(
                        ['wine', cmd + '.exe', *switches, *params],
                        stderr=subprocess.PIPE,
                        stdout=subprocess.PIPE)
                except OSError:  # we're probably running Windows Subsystem for Linux
                    # or don't have wine installed
                    proc = subprocess.run([cmd + '.exe', *switches, *params],
                                          stderr=subprocess.PIPE,
                                          stdout=subprocess.PIPE)

        else:
            proc = subprocess.run([cmd + '.exe', *switches, *params],
                                  stderr=subprocess.PIPE,
                                  stdout=subprocess.PIPE)

        if echo:
            print(proc.stdout.decode())
            print(proc.stderr.decode())

        if proc.returncode != 0:
            cmd_name = os.path.basename(cmd)
            error_msg = proc.stderr.decode()
            raise PipelineError(
                '''{} failed on with the following error message
                {}'''.format(cmd_name, error_msg))

        return proc

    def ascii2dtm(self, surfacefile, xyunits, zunits, coordsys, zone,
                  horizdatum, vertdatum, gridfile, **kwargs):
        """Converts raster data stored in ESRI ASCII raster format into a PLANS
        format data file.

        Data in the input ASCII raster file can represent a surface or raster
        data. ASCII2DTM converts areas containing NODATA values into areas with
        negative elevation values in the output data file.

        Parameters (required)
        ----------
        surfacefile: string, path to file
            Name for output canopy surface file (stored in PLANS DTM format
            with .dtm extension).
        xyunits: string
            Units for LIDAR data XY
                'M' for meters
                'F' for feet.
        zunits: string
            Units for LIDAR data elevations:
                'M' for meters
                'F' for feet
        coordsys: int
            Coordinate system for the canopy model:
                0 for unknown
                1 for UTM
                2 for state plane
        zone: int
            Coordinate system zone for the canopy model (0 for unknown).
        horizdatum: int
            Horizontal datum for the canopy model:
                0 for unknown
                1 for NAD27
                2 for NAD83
        vertdatum: int
            Vertical datum for the canopy model:
                0 for unknown
                1 for NGVD29
                2 for NAVD88
                3 for GRS80
        gridfile: string, path to file
            Name of the ESRI ASCII raster file containing surface data.

        **Kwargs (optional)
        ------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        multiplier: numeric
            Multiply all data values in the input surface by the constant.
        offset: numeric
            Add the constant to all data values in the input surface. The
            constant can be negative.
        nan: boolean
            Use more robust (but slower) logic when reading values from the
            input file to correctly parse NAN values (not a number).
        """
        cmd = 'ascii2dtm'
        params = [
            surfacefile, xyunits, zunits, coordsys, zone, horizdatum,
            vertdatum, gridfile
        ]
        self.run(cmd, *params, **kwargs)

    def asciiimport(self, paramfile, inputfile, outputfile=None, **kwargs):
        """ASCIIImport allows you to use the configuration files that describe
        the format of ASCII data files to convert data into FUSION’s LDA format.
        The configuration files are created using FUSION’s Tools... Data
        conversion... Import generic ASCII LIDAR data... menu option. This
        option allows you to interactively develop the format specifications
        needed to convert and ASCII data file into LDA format.

        Parameters
        ----------
        paramfile: string, path to file (required)
            Name of the format definition parameter file (created in FUSION's
            Tools... Data conversion... Import generic ASCII LIDAR data... menu
            option.
        inputfile: string, path to file (required)
            Name of the ASCII input file containing LIDAR data.
        outputfile: string, path to file (optional)
            Name for the output LDA or LAS file (extension will be provided
            depending on the format produced). If outputfile is omitted, the
            output file is named using the name of the input file and the
            extension appropriate for the format (.lda for LDA, .las for LAS).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        las: boolean
            Output file is stored in LAS version 1.0 format.
        -------
        additional boolean kwargs available for most FUSION commands include:
        quiet, verbose, newlog, version, locale, nolaszipdll

        Progress information for the conversion is displayed when the /verbose
        switch is used.
        -------
        """
        cmd = 'asciiimport'
        params = [paramfile, inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def canopymaxima(self, inputfile, outputfile, **kwargs):
        """Uses a canopy height model to identify local maxima using a
        variable-size evaluation window.

        The window size is based on the canopy height. For some forest types,
        this tool can identify individual trees. However, it does not work in
        all forest types and it can only identify dominant and codominant trees
        in the upper canopy. The local maxima algorithm in CanopyMaxima is
        similar to that reported in Popescu et al. (2002) and Popescu and Wynn
        (2004) and implemented in the TreeVAW software (Kini and Popescu, 2004).

        Parameters (required)
        ----------
        inputfile: string, path to file
            Name for the input canopy height model file.
        outputfile: string, path to file
            Name for the output CSV file containing the maxima.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        ground: string, path to file
            Use the specified surface model(s) to represent the ground surface.
            File may be wildcard or text list file (extension .txt).
        threshold: numeric (#)
            Limit analysis to areas above a height of # units (default: 10.0).
        wse: 4- or 6-tuple or list-like (A,B,C,D) or (A,B,C,D,E,F)
            Constant and coefficients for the variable window size equation used
            to compute the window size given the canopy surface height window:
            width = A + B*ht + C*ht^2 + D*ht^3 + E*ht^4 + F*ht^5
            Defaults values are for metric units:
            A = 2.51503, B = 0, C = 0.00901, D = 0, E = 0, F = 0.
            Use A = 8.251, B = 0, C = 0.00274, D = E = F = 0 when using
            imperial units.
        mult: numeric
            Window size multiplier (default: 1.0).
        res: numeric
            Resolution multiplier for intermediate grids (default: 2.0). A value
            of 2 results in intermediate grids with twice the number of rows and
            columns
        outxy: 4-tuple or list-like (xmin,ymin,xmax,ymax)
            Restrict output of tree located outside of the extent defined by
            (xmin,ymin) and (xmax,ymax). Tree on the left and bottom edges will
            be output, those on the top and right edges will not.
        crad: boolean
            Output 16 individual crown radii for each tree. Radii start at 3
            o'clock and are in counter-clockwise order at 22.5 degree intervals.
        shape: boolean
            Create shapefile outputs for the canopy maxima points and the
            perimeter of the area associated with each maxima.
        img8: boolean
            Create an 8-bit image showing local maxima and minima (use when 24
            bit image fails due to large canopy model).
        img24: boolean
            Create an 24-bit image showing local maxima and minima.
        new: boolean
            Create a new output file (erase output file if one exists).
        summary: boolean
            Produce a summary file containing tree height summary statistics.
        projection: string, path to file
            Associate the specified projection file with shapefile and raster
            data products.
        minmax: integer (#)
            Change the calculation method for the min/max crown width.
            Options:
                0 = report the max and min crown radii * 2
                1 = report the max and min diameters
                2 = report N-S diameter and the E-W diameter
                3 = report max diameter and the diameter perpendicular to the
                    max diameter line and the rotation to the max line
        """
        cmd = 'canopymaxima'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def canopymodel(self, surfacefile, cellsize, xyunits, zunits, coordsys,
                    zone, horizdatum, vertdatum, datafiles, **kwargs):
        """Creates a canopy surface model using a LIDAR point cloud.

        By default, the algorithm used by CanopyModel assigns the elevation of
        the highest return within each grid cell to the grid cell center.
        CanopyModel provides for smoothing of the generated surface using a
        median or a mean filter or both. Specialized logic, activated using the
        /peaks switch, preserves local maxima in the surface while smoothing to
        force the surface to adhere to the tops of trees. CanopyModel provides
        options to compute a texture metric (coefficient of variation of surface
        values within an n by n window), slope, or aspect for the canopy model
        and output them as the final surface. When used with a bare-earth model,
        CanopyModel subtracts the ground elevations from the return elevations
        to produce a canopy height model. Output from CanopyModel is a PLANS
        format DTM file that uses floating point elevation values and contains
        coordinate projection information.

        Parameters (required)
        ----------
        surfacefile: string
            Name for output canopy surface file (stored in PLANS DTM format
            with .dtm extension).
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data.
        xyunits: string
            Units for LIDAR data XY:
                'M' for meters
                'F' for feet.
        zunits: string
            Units for LIDAR data elevations:
                'M' for meters
                'F' for feet
        coordsys: int
            Coordinate system for the canopy model:
                0 for unknown
                1 for UTM
                2 for state plane
        zone: int
            Coordinate system zone for the canopy model (0 for unknown).
        horizdatum: int
            Horizontal datum for the canopy model:
                0 for unknown
                1 for NAD27
                2 for NAD83
        vertdatum: int
            Vertical datum for the canopy model:
                0 for unknown
                1 for NGVD29
                2 for NAVD88
                3 for GRS80
        datafiles: string or list-like of strings of paths to file(s)
            LIDAR data file (LDA, LAS, ASCII LIDARDAT formats)...may be wildcard
            or name of text file listing the data files. If wildcard or text
            file is used, no other datafiles will be recognized.

        **Kwargs (optional)
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        median: int
            Apply median filter to model using # by # neighbor window.
        mean: int
            Apply mean filter to model using # by # neighbor window.
        texture: int
            Calculate the surface texture metric using # by # neighbor window.
        slope: boolean
            Calculate surface slope for the final surface.
        aspect: boolean
            Calculate surface aspect for the final surface.
        outlier: 2-tuple or list-like -- (low, high)
            Omit points with elevations below low and above high if used with a
            bare-earth surface this option will omit points with heights below
            low or above high.
        multiplier: numeric
            Multiply the output values by the constant (#).
        return: string or int
            Specifies the returns to be included in the sample (can include
            A,1,2,3,4,5,6,7,8,9,F,L,O) Options are specified without commas
            (e.g. return=123) For LAS files only: F indicates first and only
            returns, L indicates last of many returns.
        las_class: string or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be used when creating the canopy
            surface. If defined as a string, classification values should be
            separated by a comma e.g. (2,3,4,5) and can range from 0 to 31.
            If the first character of string is “~”, all classes except those
            listed will be used.
        ground: string or path to file
            Use the specified bare-earth surface model(s) to normalize the LIDAR
            data. The file specifier can be a single file name, a “wildcard”
            specifier, or the name of a text file containing a list of model
            files (must have “.txt” extension). In operation, CanopyModel will
            determine which models are needed by examining the extents of the
            input point data.
        asc: boolean
            Write the output surface in ASCII raster format in addition to
            writing the surface in DTM format.
        grid: string, 4-tuple or list-like (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize. If defined as a string, use format 'X,Y,W,H'
        gridxy: string, 4-tuple or list-like (X1,Y1,X2,Y2)
            Force the origin of the output grid (lower left corner) to be
            (X1,Y1) instead of computing an origin from the data extents and
            force the upper right corner to be (X2, Y2). X2 and Y2 will be
            rounded up to a multiple of cellsize. If defined as a string, use
            format 'X1,Y1,X2,Y2'.
        align: string, path to file
            Force alignment of the output grid to use the origin (lower left
            corner), width and height of the specified dtmfile. Behavior is the
            same as /gridxy except the X1,Y1,X2,Y2 parameters are read from the
            dtmfile.
        extent: string, path to file
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file but
            adjust the origin to be an even multiple of the cell size and the
            width and height to be multiples of the cell size.
        rasterorigin: boolean
            Offset the origin and adjust the extent of the surface so raster
            data products created using the surface will align with the extent
            specified with the /grid or /gridxy options. /rasterorigin is only
            used in conjunction with the /grid or /gridxy option.
        surface: boolean
            Use the bare-earth surface model in conjunction with values
            specified in /outlier to omit points based on their height above the
            ground surface but create a surface that is not normalized relative
            to the bare-earth surface (the surface uses the point elevations).
        peaks: boolean
            Preserve localized peaks in the final surface. Only useful with
            /median or /smooth.
        pointcount: boolean
            Output the number of data points in each cell in addition to the
            canopy surface/height values. Counts are output in .DTM format. If
            there are no points for a cell, the elevation/height value for the
            cell is set to 999.0.
        nofill: boolean
            Don’t fill holes in the surface model. In general, holes result from
            a lack of data within a cell. The default behavior is to fill holes
            in the interior of the surface model.
        --------
        additional boolean kwargs available for most FUSION commands include:
        quiet, verbose, newlog, version, locale, nolaszipdll
        --------
        """
        cmd = 'canopymodel'
        params = [
            surfacefile, cellsize, xyunits, zunits, coordsys, zone, horizdatum,
            vertdatum, datafiles
        ]
        self.run(cmd, *params, **kwargs)

    def catalog(self, datafile, catalogfile=None, **kwargs):
        """Produces a set of descriptive reports describing several important
        characteristics of LIDAR data sets.

        It is most often used to evaluate a new acquisition for internal
        quality, completeness of data coverage and return or pulse density. The
        primary output of Catalog is a web page that contains a summary of all
        data tiles evaluated including attribute summaries for each tile and
        overall summaries for the entire data set. Catalog provides options that
        will create the index files needed to use the LIDAR data with FUSION
        making it the logical first step in any analysis procedure. In addition
        to the web page, Catalog can produce images representing the coverage
        area, pulse and return densities, and intensity values for the entire
        acquisition. When data are stored in LAS format, Catalog includes a
        summary of points by classification code from the LAS data. All images
        produced by Catalog have associated world files so they can be used
        within FUSION to provide a frame-of-reference for analysis. Catalog also
        produces a FUSION hotspot file that provides specific details for each
        data tile in the FUSION environment.

        Parameters
        ----------
        datafile: string, path to file (required)
            LIDAR data file template or name of a text file containing a list of
            file names (list file must have .txt extension).
        catalogfile: string, path to file (optional)
            Base name for the output catalog file (extensions will be added).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        image: boolean
            Create image files showing the coverage area for each LIDAR file.
        index: boolean
            Create LIDAR data file indexes if they don't already exist.
        newindex: boolean
            Create new LIDAR data file indexes for all files (even if they
            already exist).
        drawtiles: boolean
            Draw data file extents and names on the intensity image.
        coverage: boolean
            Create one image that shows the nominal coverage area for all data
            files included in the catalog. Also creates a FUSION hotspot file
            that provides details for each file in the catalog.
        countreturns: boolean
            Adds columns in the CSV and HTML output to show the number of
            returns by return number for each data file and all data files
            combined. Runs that use this option can take much longer to process
            because Catalog has to read every point in the data files to count
            up the different returns. In theory, LAS files have this information
            in their header. However, file produced by some version of TerraScan
            do not have these fields populated with the actual number of data
            points by return number.
        uselascounts: boolean
            Use the return counts from the header of LAS files instead of
            scanning the entire data file to count the returns. Many LAS files
            produced by TerraScan do not contain valid data for the return
            counts so make sure your data has good numbers before using this
            switch
        rawcounts: boolean
            Outputs the number of returns (or first returns) in each cell. Used
            in conjunction with the /density and /firstdensity options. The
            output is in PLANS DTM format.
        density: 3-tuple or list-like (area,min,max)
            Creates an image for all data files that shows the return density
            for the area represented by each pixel. area is the pixel area, min
            is the minimum acceptable point density, and max is the upper limit
            for the acceptable density range. Cells with point densities falling
            within the min-max range are colored green, cells with point
            densities below the minimum are colored red, and cells with
            densities above the maximum are colored blue
        firstdensity: 3-tuple or list-like (area,min,max)
            Creates an image for all data files that shows the density of first
            returns for the area represented by each pixel. area is the pixel
            area, min is the minimum acceptable point density, and max is the
            upper limit for the acceptable density range. Cells with first
            return densities falling within the min-max range are colored green,
            cells with point densities below the minimum are colored red, and
            cells with densities above the maximum are colored blue
        intensity: 3-tuple or list-like (area,min,max)
            Creates an intensity image for all data files using the average
            intensity for all first returns within each pixel. area is the pixel
            area, min is the minimum intensity value, and max is the maximum
            intensity value. A black to white color ramp is mapped to the range
            of intensity values defined by min and max. Ideally, min and max
            correspond to the range of intensity values present in the data.
            However, you may not always know the range of values for a given
            data set.
        imageextent: 4-tuple or list-like of numerics (xmin,ymin,xmax,ymax)
            Limit the area covered by image products to the specified extent.
        bmp: boolean
            Save second copy of intensity, return density, and pulse density
            images in BMP format with associated world file.
        outlier: numeric (multiplier)
            Performs a simple analysis to identify data tiles that might contain
            elevation outliers. The analysis marks tiles where the minimum,
            maximum, or range of elevations are outside the range defined by:
            mean value +- multiplier * std dev
            The default multiplier is 2.0
        las_class: string
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the class you DO NOT want included in the subsample.
            e.g. /class:~2,3 would include all class values EXCEPT 2 and 3.
            Older versions used lasclass for this switch name. /lasclass will
            still work but new scripts should use class.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        noclasssummary: boolean
            Do not create a summary table showing the number of points by LAS
            classification values. Only valid for LAS format files.
        validate: numeric (maxreturn)
            Produce report describing potential errors in point data files.
            Report will contain files with errors the might cause problems for
            other FUSION programs.
        projection: string, path to file
            Associate the specified projection file with image products.
        """
        cmd = 'catalog'
        params = [datafile, catalogfile]
        self.run(cmd, *params, **kwargs)

    def clipdata(self,
                 inputfile,
                 samplefile,
                 xmin=None,
                 ymin=None,
                 xmax=None,
                 ymax=None,
                 **kwargs):
        """Creates sub-samples of LIDAR data for various analysis tasks.

        The sub-sample can be round or rectangular and can be large or small.
        ClipData provides many of the same sampling options found in FUSION but
        it is not used by FUSION to perform subsampling of LIDAR data sets
        (FUSION has its own logic to accomplish this task). ClipData is often
        used to create sample of LIDAR returns around a specific point of
        interest such as a plot center or GPS measurement point. Subsequent
        analyses using programs like CloudMetrics facilitate comparing field
        data to LIDAR point cloud metrics. ClipData can extract a single sample
        or multiple samples using a single command. When creating several
        samples, it is much more efficient to use the optional syntax to clip
        several samples using a single command line.

        ClipData can also sub-sample data within the sample area using the
        elevation values for the returns. When used in conjunction with a
        bare-earth surface model, this logic allows for sampling a range of
        heights above ground within the sample area.

        ClipData can extract specific returns (1st, 2nd, etc) or first and last
        returns (LAS files only) for the sample area. This capability, when used
        with a large sample area, can extract specific returns from an entire
        data file.

        As part of the sampling process, ClipData can add (or subtract) a fixed
        elevation from each return elevation effecting adjusting the entire
        sample up or down. This capability, when used with a large sample area,
        can adjust entire data files up or down to help align data from
        different LIDAR acquisitions.

        Parameters
        ----------
        inputfile: string (required)
            LIDAR data file template, name of a text file containing a list of
            file names (must have .txt extension), or a FUSION Catalog CSV file.
        samplefile: string (required)
            Name for subsample file (extension will be added) or a text file
            containing sample information for 1 or more samples. Each line in
            the text file should have the subsample filename and the xmin ymin
            xmax ymax values for the sample area separated by spaces or commas.
            The output filename cannot contain spaces.
        xmin, ymin, xmax, ymax: int (optional)
            Defines the lower left corner (xmin, ymin) and upper right corner
            (xmax, ymax) of the sample area bounding box.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        shape: int
            Shape of the sample area:
                0 rectangle
                1 circle
        decimate: int
            Skip # points between included points (must be > 0).
        ground: string, path to ground dtmfile
            Use the specified bare-earth surface model to normalize the LIDAR
            data (subtract the bare-earth surface elevation from each lidar
            point elevation). Use with /zmin to include points above zmin or
            with /zmax to include points below zmax (file must be FUSION/PLANS
            format). file may be wildcard or text list file (extension .txt
            only) that specifies more than one ground surface model. In
            operation, only the models that cover the sample area will be used
            to normalize point data.
        zmin: numeric
            Include points above # elevation. Use with /dtm to include points
            above # height.
        zmax: numeric
            Include points below # elevation. Use with /dtm to include points
            below # height.
        zpercent: numeric
            Include only the upper # percent of the points. If # is (-) only the
            lower # percent of the points. # can be -100 to +100.
        height: boolean
            Convert point elevations into heights above ground using the
            specified DTM file. Always Used with /dtm.
        timemin: numeric
            Include points with GPS times greater than # (LAS only).
        timemax: numeric
            Include points with GPS times less than or equal to # (LAS only).
            Interpretation of # depends on the GPS time recorded in the LAS
            point records.
        anglemin: numeric
            Include points with scan angles greater than # (LAS only).
        anglemax: numeric
            Include points with scan angles less than or equal to # (LAS only).
        zero: boolean
            Save subsample files that contain no data points. This is useful
            when automating conversion and analysis tasks and expecting a
            subsample file every time ClipData is executed.
        biaselev: numeric
            Add an elevation offset to every LIDAR point: # can be + or -.
        return: string or int
            Specifies the returns to be included in the sample. String can
            include A,1,2,3,4,5,6,7,8,9,F,L. A includes all returns. For LAS
            files only: F indicates first and only returns, L indicates last of
            many returns. F and L will not work with non-LAS files.
        las_class: string or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be included in the subsample.
            Classification values should be separated by a comma e.g. (2,3,4,5)
            and can range from 0 to 31. If the first character of string is “~”,
            all classes except those listed will be used.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        line: numeric
            LAS files only: Only include returns from the specified flight line.
            Line numbering varies by acquisition so you need to know your data
            to specify values for the flight line number.
        noindex: boolean
            Do not use the data index files to access the data. This is useful
            when the order of the data points is important or when all returns
            for a single pulse need to stay together in the subsample file.
        index: boolean
            Create FUSION index files for the samplefile.
        lda: boolean
            Write output files using FUSION's LDA format when using LAS input
            files. The default behavior after FUSION version 3.00 is to write
            data in LAS format when the input data are in LAS format. When using
            input data in a format other than LAS, sample files are written in
            LDA format.
        nooffset: boolean
            Produce an output point file that no longer has the correct
            geo-referencing. This is used when you need to work with the point
            cloud but cannot reveal the actual location of the features
            represented in the point cloud. This option modifies the header
            values in the LAS header for the output files.
        precision: 3-tuple or list-like (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'clipdata'
        params = [inputfile, samplefile, xmin, ymin, xmax, ymax]
        self.run(cmd, *params, **kwargs)

    def clipdtm(self, inputdtm, outputdtm, xmin, ymin, xmax, ymax, **kwargs):
        """Clips a portion of the gridded surface model and stores it in a new
        file. The extent of the clipped model is specified using the lower left
        and upper right corner coordinates.

        Parameters (required)
        ----------
        inputdtm: string, path to file
            Name of the existing PLANS format DTM file to be clipped.
        outputdtm: string, path to file
            Name for the new PLANS format DTM file.
        xmin, ymin: numeric
            Lower left corner for the output DTM.
        xmax, ymax: numeric
            Upper right corner for the output DTM.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        shrink: boolean
            Shrink the extent of the input model by the amounts specified by
            xmin ymin xmax ymax. xmin is removed from left side, ymin is removed
            from bottom, xmax is removed from right side, and ymax is removed
            from top.
        multiplier: numeric
            Multiply the output values by the constant (#).
        """
        cmd = 'clipdtm'
        params = [inputdtm, outputdtm, xmin, ymin, xmax, ymax]
        self.run(cmd, *params, **kwargs)

    def cloudmetrics(self, inputfile, outputfile, **kwargs):
        """Computes a variety of statistical parameters describing a LIDAR data
        set.

        Metrics are computed using point elevations and intensity values (when
        available). In operation, CloudMetrics produces one record of output for
        each data file processed. Input can be a single LIDAR data file, a file
        template that uses DOS file specifier rules, a simple text file
        containing a list of LIDAR data file names, or a LIDAR data catalog
        produced by the Catalog utility. Output is appended to the specified
        output file unless the /new switch is used to force the creation of a
        new output data file. Output is formatted as a comma separated value
        (CSV) file that can be easily read by database, statistical, and
        MS-Excel programs.

        CloudMetrics is most often used with the output from the ClipData
        program to compute metrics that will be used for regression analysis in
        the case of plot-based LIDAR samples or for tree classification in the
        case of individual tree LIDAR samples.

        Parameters (required)
        ----------
        inputfile: string, path to file
            LIDAR data file template, name of text file containing a list of
            LIDAR data file names (must have .txt extension), or a catalog file
            produced by the Catalog utility.
        outputfile: string, path to file
            Name for output file to contain cloud metrics (using a .csv will
            associate the files with MS-Excel).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        above: numeric
            Compute various cover estimates using the specified heightbreak (#).
            See the technical detail for specific cover metrics that are
            computed.
        new: boolean
            Creates a new output file and deletes any existing file with the
            same name. A header is written to the new output file.
        firstinpulse: boolean
            Use only the first return for a pulse to compute metrics. Such
            returns may not always be labeled as return 1.
        firstreturn: boolean
            Use only first returns to compute metrics.
        highpoint: boolean
            Produce a limited set of metrics that includes only the highest
            return within the data file.
        subset: boolean
            Produce a limited set of metrics ([ID], #pts, Mean ht, Std dev ht,
            75th percentile, cover). Must be used with the /above:# option.
        id: boolean
            Parse the data file name to create an identifier for the output
            record. Data file names should include a number (e.g. sample003.lda)
            or the default identifier of 0 will be assigned to the file. The
            identifier is placed in the first column of the output record before
            the input file name.
        rid: boolean
            Parse the data file name to create an identifier for the output
            record but start at the end of the filename. Data file names can
            include any characters but the end of the name should include a
            number preceded by a non-numeric character
            (e.g. 2017_01_13_sample003.las). The identifier is placed in the
            first column of the output record before the input file name.
        minht: numeric
            Use only returns above # (use when data in the input data files have
            been normalized using a ground surface model. In older versions of
            CloudMetrics this switch was htmin.
        maxht: numeric
            Use only returns below # (use when data is normalized to ground) to
            compute metrics. The maxht is not used when computing metrics
            related to the /strata or /intstrata options.
        outlier: 2-tuple or list-like of numerics (low,high)
            Omit points with elevations below low and above high. When used with
            data that has been normalized using a ground surface, low and high
            are interpreted as heights above ground. You should use care when
            using /outlier:low,high with /minht and /maxht options. If the low
            value specified with /outlier is above the value specified with
            /minht, the value for /outlier will override the value specified for
            /minht. Similarly, if the high value specified with /outlier is less
            than the value specified for /maxht, the /outlier value will
            override the value for /maxht.
        strata: list-like of numerics [#,#,#,...]
            Count returns in various height strata. # gives the upper limit for
            each strata. Returns are counted if their height is >= the lower
            limit and < the upper limit. The first strata contains points < the
            first limit. The last strata contains points >= the last limit.
            Default strata: 0.15, 1.37, 5, 10, 20, 30.
        intstrata: list-like of numerics [#,#,#,…]
            Compute metrics using the intensity values in various height strata.
            Strata for intensity metrics are defined in the same way as the
            /strata option. Default strata: 0.15, 1.37.
        kde: 2-tuple or list-like [window,mult]
            Compute the number of modes and minimum and maximum node using a
            kernal density estimator. Window is the width of a moving average
            smoothing window in data units and mult is a multiplier for the
            bandwidth parameter of the KDE. Default window is 2.5 and the
            multiplier is 1.0
        rgb: string
            Compute intensity metrics using the color value from the RGB color
            for the returns. Valid with LAS version 1.2 and newer data files
            that contain RGB information for each return (point record types 2
            and 3). Valid color values are R, G, or B.
        relcover: boolean
            Compute the proportion of first (or all) returns above the mean and
            mode values.
        alldensity: boolean
            Use all returns when computing density (percent cover, cover above
            the mean and cover above the mode) default is to use only first
            returns when computing density.
        """
        if 'relcover' in kwargs:
            warnings.warn(
                """relcover is obsolete as of CloudMetrics version
            2.0. Metrics are computed as part of the default set of metrics.""",
                DeprecationWarning)
        if 'alldensity' in kwargs:
            warnings.warn(
                """alldensity is obsolete as of CloudMetrics version
            2.0. Metrics are computed as part of the default set of metrics.""",
                DeprecationWarning)

        cmd = 'cloudmetrics'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def cover(self, groundfile, coverfile, heightbreak, cellsize, xyunits,
              zunits, coordsys, zone, horizdatum, vertdatum, datafile,
              **kwargs):
        """Computes estimates of canopy closure using a grid.

        Output values for cover estimates range from 0.0 to 100.0 percent.
        Canopy closure us defined as the number of returns over a specified
        height threshold divided by the total number of returns within each
        cell. In addition, Cover can compute the proportion of pulses that are
        close to a bare-ground surface model to help assess canopy penetration
        by the laser scanner. Wit the addition of an upper height limit, Cover
        can compute the proportion of returns falling within specific height
        ranges providing estimates of relative vegetation density for various
        height strata.

        Parameters (required)
        ----------
        groundfile: string, path to file
            File specifier for the bare-ground surface model used to normalize
            all return elevations. The file specifier can be a single file name,
            a “wildcard” specifier, or the name of a text file containing a list
            of model files (must have “.txt” extension). In operation, Cover
            will determine which models are needed by examining the extents of
            the input point data.
        coverfile: string, path to file
            Name for the cover data file. The cover data is stored in the PLANS
            DTM format using floating point values.
        heightbreak: numeric
            Height break for the cover calculation.
        cellsize: numeric
            Grid cell size for the cover data.
        xyunits: string
            Units for LIDAR data XY:
                'M' for meters,
                'F' for feet.
        zunits: string
            Units for LIDAR data elevations:
                'M' for meters,
                'F' for feet.
        coordsys: int
            Coordinate system for the cover data:
                0 for unknown
                1 for UTM
                2 for state plane
        zone: int
            Coordinate system zone for the cover data (0 for unknown).
        horizdatum: int
            Horizontal datum for the cover data:
                0 for unknown
                1 for NAD27
                2 for NAD83
        vertdatum: int
            Vertical datum for the cover data:
                0 for unknown
                1 for NGVD29
                2 for NAVD88
                3 for GRS80
        datafiles: string, or list-like of strings, path to file(s)
            LIDAR data file (LDA, LAS, ASCII LIDARDAT formats)... may be
            wildcard or name of text file listing the data files. If wildcard or
            text file is used, no other datafile# parameters will be recognized.
            Several data files can be specified. The limit depends on the length
            of each file name. When using multiple data files, it is best to use
            a wildcard for datafile1 or create a text file containing a list of
            the data files and specifying the list file as datafile1.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        all: boolean
            Use all returns to calculate the cover data. The default is to use only first returns.
        las_class: string or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be included in the subsample.
            Classification values should be separated by a comma e.g. (2,3,4,5)
            and can range from 0 to 31. If the first character of string is “~”,
            all classes except those listed will be used.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        penetration: boolean
            Compute the proportion of returns close to the ground surface by
            counting the number of returns within +-heightbreak units of the
            ground.
        upper: numeric
            Use an upperlimit when computing the cover value. This allows you to
            calculate the proportion of returns between the heightbreak and
            upperlimit.
        """
        cmd = 'cover'
        params = [
            groundfile, coverfile, heightbreak, cellsize, xyunits, zunits,
            coordsys, zone, horizdatum, vertdatum, datafile
        ]
        self.run(cmd, *params, **kwargs)

    def csv2grid(self, inputfile, column, outputfile, **kwargs):
        """Converts data stored in comma separated value (CSV) format into ASCII
        raster format.

        In operation, users specify the column from the CSV file to convert.
        CSV2Grid expects a header file that corresponds to the input file. The
        header file name is formed from the input file name by appending the
        text “_ascii_header” and changing the extension to “.txt”. Normally, the
        CSV files used with CSV2Grid are produced by GridMetrics.

        Parameters (required)
        ----------
        inputfile: string, path to file
            Name of the input CSV file. This file is normally output from
            GridMetrics.
        column: int
            Column number for the values to populate the grid file (column
            numbers start with 1).
        outputfile: string, path to file
            Name for the output ASCII raster file.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        multiplier: numeric
            Multiply all data values by the constant (#).
        ndzero: numeric
            If the value in the target column is NODATA, look at the value in
            column # and, if it is a valid value (not NODATA), change the value
            for the target column to 0.0 for output to the ASCII grid file. This
            option is useful when the ASCII grid file will be used for further
            analysis in GIS or statistical packages.
        """
        cmd = 'csv2grid'
        params = [inputfile, column, outputfile]
        self.run(cmd, *params, **kwargs)

    def densitymetrics(self, groundfile, cellsize, slicethickness, outputfile,
                       datafiles, **kwargs):
        """Produces a series of grids where each grid contains density
        information for a specific range of heights above ground.

        Densities are reported as the proportion of the returns within the
        layer. Output consists of a CSV file with columns that correspond to the
        layers and PLANS format DTM files (one for each layer) containing the
        point density information.

        Parameters (required)
        ----------
        groundfile: string, path to file
            File specifier for the bare-ground surface model used to normalize
            all return elevations. The file specifier can be a single file name,
            a “wildcard” specifier, or the name of a text file containing a list
            of model files (must have “.txt” extension). In operation,
            DensityMetrics will determine which models are needed by examining
            the extents of the input point data.
        cellsize: numeric
            Desired grid cell size for the point density data in the same units
            as the point data.
        slicethickness: numeric
            Thickness for each “slice” in the same units as the point
            elevations.
        outputfile: string, path to file
            Base file name for output. Metrics are stored in CSV format with the
            extension .csv unless the /nocsv switch is specified, Other outputs
            are stored in files named using the base name and additional
            descriptive information.
        datafiles: string, or list-like of strings, path(s) to files
            LIDAR data file(s) (LDA, LAS, ASCII LIDARDAT formats)... may be
            wildcard or name of text file listing the data files. If wildcard or
            text file is used, no other datafile# parameters will be recognized.
            Several data files can be specified. The limit depends on the length
            of each file name. When using multiple data files, it is best to use
            a wildcard for datafile1 or create a text file containing a list of
            the data files and specifying the list file as datafile1.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        outlier: 2-tuple or list-like of numerics (low,high)
            Ignore points with elevations below low and above high. Low and high
            are interpreted as heights above ground as defined by the
            groundfile.
        maxsliceht: numeric
            Limit the range of height slices to 0 to high.
        las_class: string, or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be included when computing
            density metrics. Classification values should be separated by a
            comma e.g. (2,3,4,5) and can range from 0 to 31. If the first
            character of string is “~”, all classes except those listed will be
            used.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        nocsv: boolean
            Do not create a CSV output file for cell metrics.
        first: boolean
            Use only first returns to compute all metrics. The default is to use
            all returns to compute the metrics.
        slices: list-like of numerics [#,#,#,…]
            Use specific slice height breaks rather that evenly spaced breaks
            based on the range of heights in the data. You can specify a maximum
            of 64 slice heights. The first slice always starts at 0.0. Slice
            heights must be specified in ascending order. The highest slice will
            contain the count of all points with heights greater than or equal
            to the last height break. Slice height ranges are defined as: lower
            ht <= point height < upper height.
        grid: 4-tuple or list-like of numerics (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize.
        gridxy: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Force the origin of the output grid (lower left corner) to be
            (X1,Y1) instead of computing an origin from the data extents and
            force the upper right corner to be (X2, Y2). X2 and Y2 will be
            rounded up to a multiple of cellsize.
        align: string, path to dtmfile
            Force alignment of the output grid to use the origin (lower left
            corner), width and height of the specified dtmfile. Behavior is the
            same as /gridxy except the X1,Y1,X2,Y2 parameters are read from the
            dtmfile.
        buffer: numeric
            Add an analysis buffer of the specified width (same units as LIDAR
            data) around the data extent when computing metrics but only output
            metrics for the area specified via /grid, /gridxy, or /align. When
            /buffer is used without one of these options, metrics are output for
            an area that is inside the actual extent of the return data as
            metrics within the buffer area are not output.
        cellbuffer: int
            Add an analysis buffer specified as the number of extra rows and
            columns around the data extent when computing metrics but only
            output metrics for the area specified via /grid, /gridxy, or /align.
            When /cellbuffer is used without one of these options, metrics are
            output for an area that is inside the actual extent of the return
            data as metrics for the buffer area are not output.
        """
        cmd = 'densitymetrics'
        params = [groundfile, cellsize, slicethickness, outputfile, datafile]
        self.run(cmd, *params, **kwargs)

    def dtm2ascii(self, inputfile, outputfile=None, **kwargs):
        """Converts data stored in the PLANS DTM format into ASCII raster files.

        Such files can be imported into GIS software such as ArcInfo. DTM2ASCII
        provides the same functionality as the Tools... Terrain model... Export
        model... menu option in FUSION.

        Parameters (required)
        ----------
        inputfile: string, path to file
            Name of the PLANS DTM file to be converted into ASCII raster format.
        outputfile: string, path to file
            Name for the converted file. If outputfile is omitted, the output
            file name will be constructed from the inputfile name and the
            extension .asc. When the csv switch is specified, the default
            extension used to construct an output file name is .csv.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        csv: boolean
            Output data values in comma separated value format. Header is the
            column number. Data are arranged in rows with the northern-most row
            first in the file.
        raster: boolean
            Interpret the DTM points as the attribute for a cell and adjust the
            origin of the ASCII grid file so that the lower left data point is
            the center of the lower left grid cell. For almost all applications,
            you should use the /raster option.
        multiplier: numeric
            Multiply the output values by the constant (#).
        """
        cmd = 'dtm2acsii'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def dtm2envi(self, inputfile, outputfile=None, **kwargs):
        """Converts data stored in the PLANS DTM format into ENVI standard
        format raster files.

        Such files can be imported into GIS software such as ENVI and ArcInfo.

        Parameters
        ----------
        inputfile: string, path to file (required)
            Name of the PLANS DTM file to be converted into ASCII raster format.
        outputfile: string, path to file (optional)
            Name for the converted file. If outputfile is omitted, the output
            file name will be constructed from the inputfile name and the
            extension .nvi. The associated ENVI header file is named by
            appending “.hdr” to the inputfile name.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        south: boolean
            Specifies that data are located in the southern hemisphere
        """
        cmd = 'dtm2envi'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def dtm2tif(self, inputfile, outputfile=None, **kwargs):
        """Converts data stored in the PLANS DTM format into a TIFF image and
        creates a world file that provides coordinate system reference data for
        the image.

        Such images can be imported into GIS software or used in other analysis
        processes.

        Parameters
        ----------
        inputfile: string, path to file (required)
            Name of the PLANS DTM file to be converted into ASCII raster format.
        outputfile: string, path to file (optional)
            Name for the converted file. If outputfile is omitted, the output
            file name will be constructed from the inputfile name and the
            extension .tif.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        mask: boolean
            Produces a mask image showing the areas in the DTM with valid data
            values. In the mask image, a value of 0 indicates a cell with
            invalid data (NODATA value) and a value of 255 indicates a cell with
            a valid data value.
        """
        cmd = 'dtm2tif'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def dtm2xyz(self, inputfile, outputfile=None, **kwargs):
        """Converts data stored in the PLANS DTM format into ASCII text files
        containing XYZ points.

        Such files can be imported into GIS software as point data with the
        elevation as an attribute or used in other analysis processes.

        Parameters
        ----------
        inputfile: string, path to file (required)
            Name of the PLANS DTM file to be converted into ASCII raster format.
        outputfile: string, path to file (optional)
            Name for the converted file. If outputfile is omitted, the output
            file name will be constructed from the inputfile name and the
            extension .xyz. If the /csv switch is used, the extension will be
            .csv.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        void: boolean
            Output points from DTM with NODATA value (default is to omit).
            NODATA value is -9999.0 for the elevation.
        csv: boolean
            Output XYZ points in comma separated value format (CSV). If /csv is
            used with no outputfile, an extension of .csv will be used to form
            the output file name.
        noheader: boolean
            Do not include the column headings in CSV output files. Ignored if
            /csv is not used
        """
        cmd = 'dtm2xyz'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def dtmdescribe(self, inputfile, outputfile, **kwargs):
        """Reads header information for PLANS format DTM files and outputs the
        information to an ASCII text file compatible with most spreadsheet and
        database programs. DTMDescribe can provide information for a single file
        or multiple files.

        Parameters (required)
        ----------
        inputfile: string, path to file
            DTM file name, DTM file template, or name of a text file containing
            a list of file names (must have .txt extension).
        outputfile: string, path to file
            Name for the output ASCII CSV file. If no extension is provided, an
            extension (.csv) will be added.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        stats: boolean
            Compute descriptive statistics for the data values in the DTM
        """
        cmd = 'dtmdescribe'
        params = [inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def dtmheader(self, filename=None):
        """Launches DTMHeader, an interactive program to examine and modify
        PLANS DTM file header information.

        DTMHeader allows you to easily view and change the header information
        for a PLANS DTM file. To make it most convenient, associate the .dtm
        extension with DTMHeader so you can simply double-click a .dtm file to
        view the header.

        The values in the header that can be modified are:
            Planimetric units
            Elevation units
            Descriptive name
            Coordinate system and zone
            Horizontal datum
            Vertical datum

        Parameters (optional)
        ----------
        filename: string, path to file
            Name of the PLANS DTM file to be examined.
        """
        cmd = 'dtmheader'
        params = [filename]
        self.run(cmd, *params)

    def filterdata(self, filtertype, filterparms, windowsize, outputfile,
                   datafile, **kwargs):
        """Applies various filters to return data files to produce new return
        data files with only the returns that meet the filter requirements.

        The most common application for filtertype is to remove “outliers” from
        return data files. Other filter options overlay the return data with a
        user-specified grid and produce output return files that contain only
        the returns with the minimum or maximum elevation for each grid cell.

        Parameters (required)
        ----------
        filtertype: 'outlier', 'outlier2', 'minimum' or 'maximum'
            Filtering algorithm used to remove returns from the datafile(s).
            The following options (by name) are supported:
                outlier: removes returns above or below the mean elevation plus
                    or minus filterparms * standard deviation of the elevations
                outlier2: More robust outlier removal (experimental)
                minimum:  removes all returns except the return with the minimum
                    elevation
                maximum: removes all returns except the return with the maximum
                    elevation
        filterparms: numeric
            Parameter specific to the filtering method. For outlier this is the
            multiplier applied to the standard deviation. For minimum and
            maximum, filterparms is ignored (but a value must be included on the
            command line...use 0)
        windowsize: numeric
            Size of the window used to compute the standard deviation of
            elevations or the minimum/maximum return
        outputfile: string, path to file
            Name of the output file. If any part of the name includes spaces,
            include the entire name in double quotation marks
        datafile: string, path to file
            LIDAR data file name or template or name of a text file containing a
            list of file names (list file must have .txt extension).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        lda: boolean
            Write output files using FUSION's LDA format when using LAS input
            files. The default behavior after FUSION version 3.00 is to write
            data in LAS format when the input data are in LAS format. When using
            input data in a format other than LAS, sample files are written in
            LDA format.
        layers: boolean
            Output intermediate raster data layers
        index: boolean
            Create FUSION index files for the output file.
        minsd: numeric
            Minimum standard deviation for points within a comparison window for
            filtering to take place. Default is 1.0 (same units as elevation
            data). This switch is only useful when using the outlier filter.
        minpts: int
            Minimum number of points in the comparison window for filtering to
            take place. This option can be used with all filters but must
            specify at least 3 points when used with the outlier filter.
        minrange: numeric
            Minimum range in elevations within a window for outlier filtering to
            take place. Default is 150.0 elevation units Used only with the
            outlier2 filter.
        mingap: numeric
            Minimum vertical distance that define a gap. Used to isolate points
            above the majority of points in the filter window. Used only with
            the outlier2 filter.
        gapratio: numeric
            Proportion of points in window that can be above a vertical gap.
            Ranges from 0.0 to 1.0 Used only with the outlier2 filter.
        las_class: string or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be included in the subsample.
            Classification values should be separated by a comma e.g. (2,3,4,5)
            and can range from 0 to 31. If the first character of string is “~”,
            all classes except those listed will be used.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        precision: 3-tuple or list-like of numerics (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        reclass: int
            Change the classification code for points identified as outliers and
            write them to the output file. The optional value is the
            classification code assigned to the points. Only valid when used
            with the outlier and outlier2 filters.
        -------
        additional boolean kwargs available for most FUSION commands include:
        quiet, verbose, newlog, version, locale, nolaszipdll
        -------
        """
        cmd = 'filterdata'
        params = [filtertype, filterparms, windowsize, outputfile, datafile]
        self.run(cmd, *params, **kwargs)

    def firstlastreturn(self, outputfile, datafile, **kwargs):
        """Extracts first and last returns from a LIDAR point cloud.

        It is most commonly used when the point cloud data are provided in a
        format that does not identify the last return of a pulse.
        FirstLastReturn provided two definitions of last returns: the last
        return recorded for each pulse and the last return recorded for pulse
        with more than one return. The former includes first returns that are
        also the last return recorded for a pulse and the latter does not.

        Parameters (required)
        ----------
        outputfile: string, path to file
            Base file name for output data files. First and last returns are
            written to separate files that are named by appending
            “_first_returns” and “_last_returns” to the base file name.
        datafile: string, path to file
            LIDAR data file template or name of a text file containing a list of
            file names (list file must have .txt extension).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        index: boolean
            Create FUSION index files for the files containing the first and
            last returns.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        lastnotfirst: boolean
            Do not included first returns that are also last returns in the last
            returns output file.
        uselas: boolean
            Use information stored in the LAS point records to determine which
            returns are first and last returns.
        lda: boolean
            Write output files using FUSION's LDA format when using LAS input
            files. The default behavior after FUSION version 3.00 is to write
            data in LAS format when the input data are in LAS format. When using
            input data in a format other than LAS, sample files are written in
            LDA format.
        precision: 3-tuple or list-like of numerics (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'firstlastreturn'
        params = [outputfile, datafile]
        self.run(cmd, *params, **kwargs)

    def gridmetrics(self, groundfile, heightbreak, cellsize, outputfile,
                    datafiles, **kwargs):
        """Computes a series of descriptive statistics for a LIDAR data set.

        Output is a raster (grid) represented in database form with each record
        corresponding to a single grid cell. GridMetrics is similar to
        CloudMetrics except it computes metrics for all returns within each cell
        in the output grid. Cloudmetrics computes a single set of metrics for
        the entire data set. The default output from GridMetrics is an ASCII
        text file with comma separated values (CSV format). Field headings are
        included and the files are easily read into database and spreadsheet
        programs. Optionally, GridMetrics can output raster layers stored in
        PLANS DTM format. GridMetrics compute statistics using both elevation
        and intensity values in the same run. GridMetrics can apply the fuel
        models developed to predict canopy fuel characteristics in Douglas-fir
        forest type in Western Washington (Andersen, et al. 2005). Application
        of the fuel models to other stand types or geographic regions may
        produce questionable results.

        Parameters (required)
        ----------
        groundfile: string, path to file
            Name for ground surface model(s) (PLANS DTM with .dtm extension)...
            may be wildcard or name of text file listing the data files.
            Multiple ground models can be used to facilitate processing of large
            areas where a single model for the entire acquisition is too large
            to hold in memory.
        heightbreak: numeric
            Height break for cover calculation.
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data.
        outputfile: string
            Base name for output file. Metrics are stored in CSV format with
            .csv extension unless the /nocsv switch is used. Other outputs are
            stored in files named using the base name and additional descriptive
            information.
        datafiles: string or list-like of strings with path(s) to lidar file(s)
            First LIDAR data file (LDA, LAS, ASCII LIDARDAT formats)...may be
            wildcard or name of text file listing the data files. If wildcard or
            text file is used, no other datafile# parameters will be recognized.
            Several data files can be specified. The limit depends on the length
            of each file name. When using multiple data files, it is best to use
            a wildcard for datafile1 or create a text file containing a list of
            the data files and specifying the list file as datafile1.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        outlier: 2-tuple (low,high)
            Omit points with elevations below low and above high. low and high
            are interpreted as heights above ground.
        las_class: string or list-like
            Used with LAS format files only. Specifies that only points with
            classification values listed are to be included in the subsample.
            Classification values should be separated by a comma e.g. (2,3,4,5)
            and can range from 0 to 31. If the first character of string is “~”,
            all classes except those listed will be used.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        id: string
            Include the identifier string as the last column in every record in
            the outputfile. The identifier will be included in all files
            containing metrics (elevation, intensity, and topo). The identifier
            cannot include spaces.
        minpts: int
            Minimum number of points in a cell required to compute metrics
            (default is 4 points).
        minht: numeric
            Minimum height used for points used to compute metrics. Density is
            always computed using all points including those with heights below
            the specified minimum height.
        nocsv: boolean
            Do not create an output file for cell metrics.
        noground: boolean
            Do not use a ground surface model. When this option is specified,
            the groundfile parameter should be omitted from the command line.
            Cover estimates, densitytotal, densityabove, and densitycell metrics
            are meaningless when no ground surface model is used unless the
            LIDAR data have been normalize to the ground surface using some
            other process.
        nointdtm: boolean
            Do not create an internal DTM to use when normalizing point
            elevations. The default behavior is to create an internal model that
            corresponds to the extent of the point data (or the area specified
            using the /grid, /gridxy, or /align switches). In some cases,
            creating the internal model causes problems with processing. Most
            often this caused problems for small areas with metrics being
            computed for a large cell size. The internal model was created to
            cover a slightly larger area than the data extent resulting in bad
            metrics along the top and right sides of the data extent.
        diskground: boolean
            Do not load ground surface models into memory. When this option is
            specified, larger areas can be processed but processing may be 4 to
            5 times slower. Ignored when /noground option is used.
        alldensity: boolean
            This switch is obsolete as of GridMetrics version 3.0.
            Use all returns when computing density (percent cover) default is to
            use only first returns when computing density.
        first: boolean
            Use only first returns to compute all metrics. Default is to use all
            returns for metrics.
        intensity: boolean
            This switch is obsolete as of GridMetrics version 3.0.
            Compute metrics using intensity values (default is elevation).
        nointensity: boolean
            Do not compute metrics using intensity values (default is to
            compute metrics using both intensity and elevation values).
        rgb: string
            Compute intensity metrics using the color value from the RGB color
            for the returns. Valid with LAS version 1.2 and newer data files
            that contain RGB information for each return (point record types 2
            and 3). Valid color values are R, G, or B.
        fuel: boolean
            Apply fuel parameter models (cannot be used with /first switch).
        grid: 4-tuple or list-like of numerics (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize.
        gridxy: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Force the origin of the output grid (lower left corner) to be
            (X1,Y1) instead of computing an origin from the data extents and
            force the upper right corner to be (X2, Y2). X2 and Y2 will be
            rounded up to a multiple of cellsize.
        align: string, path to dtmfile
            Force alignment of the output grid to use the origin (lower left
            corner), width and height of the specified dtmfile. Behavior is the
            same as /gridxy except the X1,Y1,X2,Y2 parameters are read from the
            dtmfile.
        extent: string, path to dtmfile
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file but
            adjust the origin to be an even multiple of the cell size and the
            width and height to be multiples of the cell size.
        buffer: numeric
            Add an analysis buffer of the specified width (same units as LIDAR
            data) around the data extent when computing metrics but only output
            metrics for the area specified via /grid, /gridxy, or /align. When
            /buffer is used without one of these options, metrics are output for
            an area that is inside the actual extent of the return data as
            metrics within the buffer area are not output.
        cellbuffer: numeric
            Add an analysis buffer specified as the number of extra rows and
            columns around the data extent when computing metrics but only
            output metrics for the area specified via /grid, /gridxy, or /align.
            When /cellbuffer is used without one of these options, metrics are
            output for an area that is inside the actual extent of the return
            data as metrics for the buffer area are not output.
        strata: list-like of numerics [#,#,#,…]
            Count returns in various height strata. # gives the upper limit for
            each strata. Returns are counted if their height is >= the lower
            limit and < the upper limit. The first strata contains points < the
            first limit. The last strata contains points >= the last limit.
            Default strata: 0.15, 1.37, 5, 10, 20, 30.
        intstrata: list-like of numerics [#,#,#,…]
            Compute metrics using the intensity values in various height strata.
            Strata for intensity metrics are defined in the same way as the
            /strata option. Default strata: 0.15, 1.37.
        kde: 2-tuple or list-like (window,mult)
            Compute the number of modes and minimum and maximum node using a
            kernel density estimator. Window is the width of a moving average
            smoothing window in data units and mult is a multiplier for the
            bandwidth parameter of the KDE. Default window is 2.5 and the
            multiplier is 1.0
        asc: boolean
            Store raster files in ASCII raster format for direct import into
            ArcGIS. Using this option preserves metrics with negative values.
            Such values are lost when raster data are stored using the PLANS DTM
            format. This switch is ignored unless it is used with the /raster
            switch.
        topo: 2-tuple or list-like (dist,lat)
            Compute topographic metrics using the groundfile(s) and output them
            in a separate file. Distance is the cell size for the 3 by 3 cell
            analysis area and lat is the latitude (+north, -south).
        raster: string or list-like of strings
            Create raster files containing the point cloud metrics. layers is a
            list of metric names separated by commas. Raster files are stored in
            PLANS DTM format or ASCII raster format when the /ascii switch is
            used. Topographic metrics are not available with the /raster:layers
            switch.

            Available metrics are:
                count: Number of returns above the minimum height
                densitytotal: total returns used for calculating cover
                densityabove: returns above heightbreak
                densitycell: Density of returns used for calculating cover
                min: minimum value for cell
                max: maximum value for cell
                mean: mean value for cell
                mode: modal value for cell
                stddev: standard deviation of cell values
                variance: variance of cell values
                cv: coefficient of variation for cell
                cover: cover estimate for cell
                abovemean: proportion of first (or all) returns above the mean
                abovemode: proportion of first (or all) returns above the mode
                skewness: skewness computed for cell
                kurtosis: kurtosis computed for cell
                AAD: average absolute deviation from mean for the cell
                p01: 1st percentile value for cell (must be p01, not p1)
                p05: 5th percentile value for cell (must be p05, not p5)
                p10: 10th percentile value for cell
                p20: 20th percentile value for cell
                p25: 25th percentile value for cell
                p30: 30th percentile value for cell
                p40: 40th percentile value for cell
                p50: 50th percentile value (median) for cell
                p60: 60th percentile value for cell
                p70: 70th percentile value for cell
                p75: 75th percentile value for cell
                p80: 80th percentile value for cell
                p90: 90th percentile value for cell
                p95: 95th percentile value for cell
                p99: 99th percentile value for cell
                iq: 75th percentile minus 25th percentile for cell
                90m10: 90th percentile minus 10th percentile for cell
                95m05: 95th percentile minus 5th percentile for cell
                r1count: Count of return 1 points above the minimum height
                r2count: Count of return 2 points above the minimum height
                r3count: Count of return 3 points above the minimum height
                r4count: Count of return 4 points above the minimum height
                r5count: Count of return 5 points above the minimum height
                r6count: Count of return 6 points above the minimum height
                r7count: Count of return 7 points above the minimum height
                r8count: Count of return 8 points above the minimum height
                r9count: Count of return 9 points above the minimum height
                rothercount: Count of other returns above the minimum height
                allcover: (all returns above cover ht) / (total returns)
                afcover: (all returns above cover ht) / (total first returns)
                allcount: number of returns above cover ht
                allabovemean: (all returns above mean ht) / (total returns)
                allabovemode: (all returns above ht mode) / (total returns)
                afabovemean: (all returns above mean ht) / (total first returns)
                afabovemode: (all returns above ht mode) / (total first returns)
                fcountmean: number of first returns above mean ht
                fcountmode: number of first returns above ht mode
                allcountmean: number of returns above mean ht
                allcountmode: number of returns above ht mode
                totalfirst: total number of 1st returns
                totalall: total number of returns
        -------
        additional boolean kwargs available for most FUSION commands include:
        quiet, verbose, newlog, version, locale, nolaszipdll
        -------
        """
        if 'raster' in kwargs:
            warnings.warn(
                """It is likely that the /raster option will be removed at some
            point. The amount of code required to implement this option is quite
            large and the DTM files produced by the option cannot support
            negative numbers. There are better ways to get individual metrics
            for input into other analyses. For example, you can use the CSV2Grid
            utility to extract specific columns from the .CSV files produced by
            GridMetrics. Use of the /raster option is discouraged.""",
                PendingDeprecationWarning)

        cmd = 'gridmetrics'
        params = [groundfile, heightbreak, cellsize, outputfile, datafiles]
        self.run(cmd, *params, **kwargs)

    def gridsample(self, gridfile, inputfile, outputfile, windowsize,
                   **kwargs):
        """Produces a comma separated values (CSV) file that contains values for
        the grid cells surrounding a specific XY location.

        Input is a file containing a list of XY coordinate pairs, one pair per
        line of the file. The user specifies the size of the sample window on
        the command line. Output is a CSV file containing the original XY
        location and the grid values from the sample window.

        Parameters
        ----------
        gridfile: string, path to file
            Name for the input grid file (PLANS DTM format)
        inputfile: string, path to file
            Name of file containing XY locations (space or comma delimited)
        outputfile: string, path to file
            Name for the output ASCII CSV file
        windowsize: int
            Size of the sample window in grid cells (must be an odd number)

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        center: boolean
            Include the location of the center of the cell that contains the
            sample point in the output CSV file
        """
        cmd = 'gridsample'
        params = [gridfile, inputfile, outputfile, windowsize]
        self.run(cmd, *params, **kwargs)

    def gridsurfacecreate(self, surfacefile, cellsize, xyunits, zunits,
                          coordsys, zone, horizdatum, vertdatum, datafile,
                          **kwargs):
        """Creates a gridded surface model using collections of random points.

        The surface model is stored in PLANS DTM format using floating point
        elevation values. Individual cell elevations are calculated using the
        average elevation of all points within the cell. GridSurfaceCreate is
        most often used with bare-earth point sets obtained from LIDAR vendors
        or by using the GroundFilter program.

        Parameters
        ----------
        surfacefile: string, path to file
            Name for output surface file (stored in PLANS DTM format with .dtm
            extension)
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data
        xyunits: string
            Units for LIDAR data XY
                'M' for meters
                'F' for feet.
        zunits: string
            Units for LIDAR data elevations
                'M' for meters
                'F' for feet.
        coordsys: int
            Coordinate system for LIDAR data:
                0 for unknown
                1 for UTM
                2 for state plane)
        zone: int
            Coordinate system zone for LIDAR data (0 for unknown)
        horizdatum: int
            Horizontal datum:
                 0 for unknown
                 1 for NAD27
                 2 for NAD83
        vertdatum: int
            Vertical datum:
                 0 for unknown
                 1 for NGVD29
                 2 for NAVD88
                 3 for GRS80
        datafiles: string, or list-like of strings, path(s) to file(s)
            First LIDAR file (LDA, LAS, ASCII XYZ, ASCII XYZI formats)... may be
            wildcard and omit other data file parameters. If wildcard or text
            file is used, no other datafile# parameters will be recognized.
            Several data files can be specified. The limit depends on the length
            of each file name. When using multiple data files, it is best to use
            a wildcard for datafile1 or create a text file containing a list of
            the data files and specifying the list file as datafile1.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        The order of the switches is important...the first filter specified on
        the command line will be the first filter applied to the model (median
        or smooth)

        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        median: int (#)
            Apply median filter to model using # by # neighbor window
        smooth: int (#)
            Apply mean filter to model using # by # neighbor window
        slope: numeric (#)
            Filter areas from the surface with slope greater than # percent.
            Slope filtering takes place after all other smoothing operations
        spike: numeric (#)
            Filter to remove spikes with slopes greater than # percent. Spike
            filtering takes place after slope filtering
        residuals: boolean
            Compute residual statistics for all points
        filldist: int (#)
            Maximum search radius (in cells) used when filling holes in the
            surface. Default is 99 cells.
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        minimum: boolean
            Use the minimum elevation value in cells to create the surface
        maximum: boolean
            Use the maximum elevation value in cells to create the surface
        grid: 4-tuple or list-like of numerics (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize
        gridxy: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Force the origin of the output grid to be (X1,Y1) instead of
            computing an origin from the data extents and force the grid to use
            (X2,Y2) as the upper right corner of the coverage area. The actual
            upper right corner will be adjusted to be a multiple of cellsize
        align: string, path to dtmfile
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file
        extent: string, path to dtmfile
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file but
            adjust the origin to be an even multiple of the cell size and the
            width and height to be multiples of the cell size.
        """
        cmd = 'gridsurfacecreate'
        params = [
            surfacefile, cellsize, xyunits, zunits, coordsys, zone, horizdatum,
            vertdatum, datafile
        ]
        self.run(cmd, *params, **kwargs)

    def gridsurfacestats(self, inputfile, outputfile, samplefactor, **kwargs):
        """Computes the surface area and volume under a surface (or between the
        surface and a ground surface).

        When used with a canopy height or surface model, it provides information
        useful for describing canopy surface roughness and the volume occupied
        by tree canopies.

        Parameters
        ----------
        inputfile: string, path to file
            File specifier for the input surface file. This can be a single
            file, a wildcard specifier, or a text list file (extension .txt
            only).
        outputfile: string, path to file
            Base name for the output files containing the surface statistics
            including the .dtm extension. If the /ascii switch is used, the
            output files will all be in ASCII raster format with the .asc
            extension.
        samplefactor: numeric
            Multiplier for outputfile cell size. outputfile cells will represent
            samplefactor * samplefactor cells from the inputfile. When multiple
            input files are used, the cell size of the first file is used to
            compute the outfile cell size.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        ground: string, path to file
            Use the specified surface model to represent the ground surface file
            may be wildcard or text list file (extension .txt only)
        asc: boolean
            Output all files in ASCII raster format with the .asc extension
        area: boolean
            Compute the surface area of inputfile instead of the surface area
            divided by the flat cell area
        halfcell: boolean
            Force alignment of the output grid to match the grid used by
            GridMetrics. This option cannot be used with the /grid, /gridxy,
            /extent, or /align switches
        svonly: boolean
            Output only the surface volume metric layer.
        grid: 4-tuple or list-like of numerics (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize
        gridxy: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Force the origin of the output grid to be (X1,Y1) instead of
            computing an origin from the data extents and force the grid to use
            (X2,Y2) as the upper right corner of the coverage area. The actual
            upper right corner will be adjusted to be a multiple of cellsize
        align: string, path to dtmfile
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file
        extent: string, path to dtmfile
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file but
            adjust the origin to be an even multiple of the cell size and the
            width and height to be multiples of the cell size.
        """
        cmd = 'gridsurfacestats'
        params = [inputfile, outputfile, samplefactor]
        self.run(cmd, *params, **kwargs)

    def groundfilter(self, outputfile, cellsize, datafile, **kwargs):
        """Filters a cloud of LIDAR returns to identify those returns that lie
        on the probable ground surface (bare-earth points).

        GroundFilter does not produce a perfect set of bare-earth returns in
        that it does not completely remove returns from large, relatively flat,
        elevated surface such as building roofs. Most vegetation returns are
        removed with the appropriate coefficients for the weight function and
        sufficient iterations. Experimentation with GroundFilter has shown that
        the default coefficients for the weight function produce good results in
        high-density point clouds (> 4 returns/sq m). The program can be used
        with low-density point clouds but some experimentation may be needed to
        select appropriate coefficients. In general, GroundFilter produces point
        sets that result in surface models that are adequate for calculating
        vegetation heights. The point set and resulting models may not be
        adequate when the bare-earth surface is the primary product.

        The output from GroundFilter is a file containing only the points
        classified as ground returns stored in LDA format. The output file can
        be used with the GridSurfaceCreate or TINSurfaceCreate utilities to
        produce a ground surface model.

        Parameters
        ----------
        outputfile: string, path to file
            Name for the output point data file (stored in LDA format)
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data. This is used
            for intermediate surfaces and is not the cell size for the final
            ground surface model created using GridSurfaceCreate.
        datafiles: string or list-like of strings, path(s) to file(s)
            LIDAR data file(s) (LDA, LAS, ASCII XYZ, ASCII XYZI format)... may
            be wildcard or text file listing the data files (if using wildcard
            or text file, omit other data file parameters). You can specify as
            many data files as you need but a wildcard specifier may be more
            efficient

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        The order of the median and smooth switches is important...the first
        filter specified will be the first filter applied to the intermediate
        surface model

        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        lda: boolean
            Store output data file in FUSION's LDA format
        surface: boolean
            Create a surface model using the final ground points
        median: numeric (#)
            Apply median filter to intermediate surface model using # by #
            window (default is no median filter)
        smooth: numeric (#)
            Apply mean filter to intermediate surface model using # by # window
            (default is no smoothing)
        finalsmooth: boolean
            Apply smoothing after the final iteration before selecting bare-
            earth points... only used when smooth or median switch is used
        outlier: 2-tuple or list-like of numerics (low,high)
            Omit points with elevations below low and above high
        gparam: numeric (#)
            Value of g parameter (default is -2.0)
        wparam: numeric (#)
            Value of w parameter (default is 2.5)
        aparam: numeric (#)
            Value of a parameter (default is 1.0)
        bparam: numeric (#)
            Value of b parameter (default is 4.0)
        tolerance: numeric (#)
            Tolerance value for final filtering of ground points. If not
            specified, weight values are used to determine final ground points.
            When specified, points with elevations less than or equal to the
            final surface elevation plus the tolerance are included in the
            ground point data set.
        iterations: int (#)
            Number of iterations for the filtering logic (default is 5)
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        extent: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Only consider points within the extent defined by (X1,Y1) and
            (X2, Y2) for filtering. /extent determines which points are
            considered for filtering.
        trim: 4-tuple or list-like of numerics (X1,Y1,X2,Y2)
            Only output points within the extent defined by (X1,Y1) and
            (X2, Y2). /trim is used along with /extent to allow filtering using
            points within a larger extent but only output points within the
            smaller extent. This minimizes edge artifacts in the final point set
            and surface created using the points. /trim determines which
            filtered points are output.
        diagnostics: boolean
            Display diagnostic information and produce diagnostic files
        precision: 3-tuple or list-like of numerics (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'groundfilter'
        params = [outputfile, cellsize, datafile]
        self.run(cmd, *params, **kwargs)

    def imagecreate(self, imagefilename, pixelsize, datafiles, **kwargs):
        """Creates an image from LIDAR data using the intensity value or
        elevation of the highest return within an image pixel.

        Optionally uses the height above a surface model to create the image.
        The output image is geo-referenced using a world file. The extent of the
        image is computed so that the image origin is a multiple of the pixel
        size. The default image file format is JPEG.

        Parameters
        ----------
        imagefilename: string, path to file
            Name for the output image file. The file will be stored in the
            specified format regardless of the extension.
        pixelsize: numeric
            Size (in the same units as the LIDAR data) of each pixel in the
            output image.
        datafile: string, path to file
            LIDAR data files stored in binary LDA or LAS formats or ASCII
            LIDARDAT format (LIDARDAT format may not be supported in future
            versions).

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        bmp: boolean
            Store the output image file in Windows BMP format and set output
            image file extension to “.bmp”.
        jpeg: boolean
            Store the output image file in JPEG format and set output image file
            extension to “.jpg”.
        coloroption: int
            Method used to assign color to each image pixel.
            Valid values are:
                0: Assign color using intensity
                1: Assign color using elevation
                2: Assign color using height above surface
        dtm: string, path to dtmfile
            Name of the surface file used to compute heights. Only PLANS format
            surface models are recognized.
        minimum: numeric
            Minimum value used to constrain the color assigned to a pixel.
            Returns with values less than value will be colored using the
            starting color.
        maximum: numeric
            Maximum value used to constrain the color assigned to a pixel.
            Returns with values greater than value will be colored using the
            ending color.
        startcolor: numeric or 3-tuple or list-like of numerics
            Starting color in the color ramp used to assign pixel colors. Value
            can be a single number representing a combined RGB color or a series
            of three values separated by commas representing the R, G, and B
            color components.
        stopcolor: numeric or 3-tuple or list-like of numerics
            Ending color in the color ramp used to assign pixel colors. Value
            can be a single number representing a combined RGB color or a series
            of three values separated by commas representing the R, G, and B
            color components.
        hsv: boolean
            Use the HSV color model to create the color ramp.
        rgb: boolean
            Use the RGB color model to create the color ramp.
        backgroundcolor: numeric or 3-tuple or list-like of numerics
            Background color for areas of the image not covered by LIDAR data.
            Value can be a single number representing a combined RGB color or a
            series of three values separated by commas representing the R, G,
            and B color components.
        """
        cmd = 'imagecreate'
        params = [imagefilename, pixelsize, datafile]
        self.run(cmd, *params, **kwargs)

    def intensityimage(self, cellSize, imagefile, datafile, **kwargs):
        """Creates images using the intensity values from a point cloud.

        Airborne laser scanning (commonly referred to as LIDAR) data have proven
        to be a good source of information for describing the ground surface and
        characterizing the size and extent of man-made features such as road
        systems and buildings. The technology has gained a strong foothold in
        mapping operations traditionally dominated by photogrammetric
        techniques. In the forestry context, airborne laser scanning data have
        been used to produce detailed bare-earth surface models and to
        characterize vegetation structure and spatial extent. Hyyppä et al.
        (2004) provide an overview of LIDAR applications for forest
        measurements. One often overlooked component of LIDAR data, return
        intensity, is seldom used in analysis procedures. LIDAR return intensity
        is related to the ratio of the amount of the laser energy detected by
        the receiver for a given reflection point to the amount of total energy
        emitted for the laser pulse (Baltsavias, 1999; Wehr and Lohr, 1999).
        Because this ratio is quite small (Baltsavias, 1999), intensity values
        reported in LIDAR data are scaled to a more useful range (8-bit values
        are common). Intensity values are collected by most sensors in use today
        and providers include the intensity values in point cloud data files
        stored in the standard LAS format (ASPRS, 2005). Flood (2001) points out
        that while intensity data have been available for some time, their use
        in commercial data processing workflows is limited. Song et al. (2002)
        evaluated the potential for identifying a variety of surface materials
        (asphalt, grass, house roof, and trees) using LIDAR intensity in an
        urban environment. Charaniya et al. (2004) used LIDAR point data, LIDAR
        intensity, USGS 10m-resolution digital elevation model, and
        black-and-white ortho-photographs to classify LIDAR points into the same
        categories. Hasegawa (2006) conducted experiments to investigate the
        effects of target material, scan geometry and distance-to-target on
        intensity values using a typical airborne scanner attached to a fixed
        mount on the ground. He also evaluated the use of airborne LIDAR data to
        identify a variety of materials. He concluded that some materials were
        easily separated (soil, gravel, old asphalt, and grass) while others
        were not easy to separate (cement, slate, zinc, brick, and trees).
        Brennan and Webster (2006) utilized a rule-based object-oriented
        approach to classify a variety of land cover types using LIDAR height
        and intensity data. They found that both height and intensity
        information were needed to separate and classify ten land cover types.
        Brandtberg (2007) also found that use of intensity data significantly
        improved LIDAR detection and classification of leaf-off eastern
        deciduous forests.

        Parameters
        ----------
        cellSize: numeric
            Pixel size for the intensity image (same units as LIDAR data)
        imagefile: string, path to file
            Name for the image file
        datafile: string, path to file
            LIDAR data file template or name of a text file containing a list of
            file names (must have .txt extension)

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        minint: numeric
            Minimum intensity percentile used for the image (default: 2.0)
        maxint: numeric
            Maximum intensity percentile used for the image (default: 98.0)
        intrange: 2-tuple or list-like of numerics (min,max)
            Force the scaling of intensity values using the specified min and
            max values. Setting the min value to -1 should force the output
            range to start with 1 instead of 0. Doing this in combination with
            /void:0,0,0 will allow you to identify areas with no data in the
            output image as they will have a value of 0 for all color bands.
        intcell: numeric
            Cell size multiplier for intermediate images
        void: 3-tuple of ints (R,G,B)
            Color for areas with no data (default is red (255,0,0))
        allreturns: boolean
            Use all returns to create the intensity image
        lowest: boolean
            Use lowest return in pixel to assign intensity value...should be
            used with /allreturns for best effect
        lowall: boolean
            Combines the /lowest and /allreturns switches...will have no effect
            if used with either the /lowest or /allreturns
        saveint: boolean
            Save the intermediate image files
        rasterorigin: boolean
            Force alignment to match other raster products generated from point
            data (offsets the origin of the image by 1/2 pixel)
        diskonly: boolean
            Do not attempt to read all returns into memory for processing
        hist: boolean
            Produce the intensity histogram data files
        jpg: boolean
            Save the intensity image using the JPEG format (default is BMP)
        projection: string, path to file
            Associate the specified projection file with image products.
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpreted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        """
        cmd = 'intensityimage'
        params = [cellSize, imagefile, datafile]
        self.run(cmd, *params, **kwargs)

    def joindb(self, basefile, basefield, addfile, addfield, startfield,
               outputfile, **kwargs):
        """Merges two CSV files into a single file using a key field.

        This is a JOIN operation in database programs.

        Parameters
        ----------
        basefile: string, path to file
            Primary data file in CSV format.
        basefield: string
            Field in the primary data file that will be matched to records in
            the secondary data file.
        addfile: string, path to file
            Secondary data file in CSV format.
        addfield: string
            Field in the secondary data file that will be matched to records in
            the primary data file.
        startfield: string
            Starting field in the secondary data file. All fields (columns)
            starting with the startfield will be added to records in the
            outputfile.
        outputfile: string, path to file
            Name for the new data file. The extension specifies the desired
            format (.csv). outputfile can be the same as basefile. outputfile
            cannot be the same as addfile.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        noheader: boolean
            Treat the first line of the basefile and addfile as data. Default
            behavior assumes the first line of each file contains column names.
            Valid for CSV format files only.
        """
        cmd = 'joindb'
        params = [
            basefile, basefield, addfile, addfield, startfield, outputfile
        ]
        self.run(cmd, *params, **kwargs)

    def lda2ascii(self,
                  inputfile,
                  outputfile,
                  format,
                  identifier=None,
                  noheader=None):
        """Converts LIDAR data files into ASCII text files.

        Provides simple conversion capabilities for all LIDAR formats supported
        by FUSION. It was originally developed as a way to check the values
        being stored in LDA format files but still has utility as a conversion
        tool. It provides capabilities similar to the Tools...Data conversion...
        Export data from LAS or LDA formats to other formats... menu option.

        Parameters
        ----------
        inputfile: string, path to file (required)
            LDA LIDAR data file
        outputfile: string, path to file (required)
            ASCII XYZ output file
        format: int (required)
            format identifier:
                0=XYZ
                1=Pulse Return X Y Z Nadir Intensity
                2=LTK CSV
                3=LTK CSV with identifier
        identifier: int (optional)
            identifier for format 3 output, cannot include spaces
        noheader: int (optional)
            if zero, suppress the heading in the output file. Can only be used
            when [identifier] is used
        """
        cmd = 'lda2ascii'
        params = [inputfile, outputfile, format, identifier, noheader]
        self.run(cmd, *params)

    def mergedata(self, datafile, outputfile, **kwargs):
        """MergeData combines several point cloud files into a single output
        file.

        The merge is accomplished by sequentially reading each input file and
        writing the point data to the output file.

        Parameters (required)
        ----------
        datafile: string, path to file
            LIDAR data file template or name of a text file containing a list of
            file names (list file must have .txt extension).
        outputfile: string, path to file
            Name for the output data file with extension.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        lda: boolean
            Store output data file in FUSION's LDA format
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        index: boolean
            Create FUSION index files for output data file
        precision: 3-tuple or list-like (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'mergedata'
        params = [datafile, outputfile]
        self.run(cmd, *params, **kwargs)

    def mergedtm(self, outputfile, inputfile, **kwargs):
        """MergeDTM combines several PLANS format DTM files into a single output
        file.

        MergeDTM can merge hundreds of DTM files and is not limited by the
        amount of memory available on the computer, only the amount of disk
        space on the output device. MergeDTM provides the same capability as the
        Tools... Terrain model... Combine... menu option in FUSION except
        MergeDTM does not automatically fill voids areas in the final output
        model (FUSION provides an option to do this).

        Parameters (required)
        ----------
        outputfile: string, path to file
            Name for the output .DTM file containing the merged data
        inputfile: string, path to file
            File name template for .DTM files to be merged, a list of .DTM file
            names, a text file containing a list of .DTM files

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        cellsize: numeric (size)
            Resample the data using a size by size pixel
        overlap: string (operator)
            Specify how overlap areas should be treated.
            Operators are: average, min, max, add, new
            The new operator populates a pixel using the value from the last
            .DTM file containing valid data for the pixel
        disk: boolean
            Merge the .DTM files to a disk file. The default behavior is to try
            to hold the merged model into memory but there can be problems when
            there is not quite enough memory for the model.
        precision: int
            Override the default precision for the merged output file. The
            default behavior is to use the highest precision of the input models
            to establish the precision of the output model.
            Valid values for precision are:
                0     2-byte integer
                1     4-byte integer
                2     4-byte floating point (C type: float)
                3     8-byte floating point (C type: double)
        exactextent: boolean
            Preserve the exact extent of the input models in the output model.
            The default behavior is to expand the extent to the nearest multiple
            of the output cell size.
        nofill: boolean
            Do not fill holes in the merged DTM.
        """
        cmd = 'mergedtm'
        params = [outputfile, inputfile]
        self.run(cmd, *params, **kwargs)

    def mergeraster(self, outputfile, inputfile, **kwargs):
        """MergeRaster combines several ASCII Raster format files into a single
        output file.

        MergeRaster can merge hundreds of ASCII Raster files but is limited by
        the amount of memory available on the computer.

        Parameters (required)
        ----------
        outputfile: string, path to file
            Name for the output ASCII raster file containing the merged data.
        inputfile: string, path to file
            File name template for ASCII raster files to be merged, a list of
            ASCII raster file names, or text file containing a list of ASCII
            raster file names with the ".txt" extension. If you are specifying a
            single input file on the command line, the file extension cannot be
            ".txt".

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        overlap: string
            Specify how overlap areas should be treated.
            Operators are: average, min, max, add, new
            The new operator populates a pixel using the value from the last
            ASCII raster file containing valid data for the pixel.
        compare: boolean
            Reports if values in cells common to two or more input files are
            different by more than 0.001.
        precision: int (#)
            Output value with # digits to the right of the decimal point.
            Default precision is 4 digits.
        nodata: numeric (#)
            Use value (#) to indicate no data in the output file. Default NODATA
            value is -9999.
        """
        cmd = 'mergeraster'
        params = [outputfile, inputfile]
        self.run(cmd, *params, **kwargs)

    def pdq(self, datafile=None, **kwargs):
        """PDQ is a simple, fast data viewer for .LDA, .LAS, and .DTM files.

        PDQ supports drag-and-drop so you can drag data files onto an icon
        (shortcut) for PDQ or you can drop data files onto a running instance of
        PDQ. For point cloud data, PDQ automatically applies a color ramp to the
        data using the point elevations. The color ramp runs from brown (lowest)
        to green (highest).

        Parameters (optional)
        ----------
        datafile: string, path to file
            Name of the input data file. File must be in LDA, LAS, or PLANS DTM
            format.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        m: boolean
            Allows multiple instances of PDQ. To use this option, you must run
            PDQ from a DOS command line.
        s: boolean
            Synchronize multiple instances of PDQ so data manipulation in one
            instance will be mirrored in all other instances. To use this
            option, you must run PDQ from a DOS command line.
        """
        cmd = 'pdq'
        params = [datafile]
        self.run(cmd, *params, **kwargs)

    def polyclipdata(self, polyfile, outputfile, datafile, **kwargs):
        """PolyClipData clips point data using polygons stored in ESRI
        shapefiles.

        The default behavior of PolyClipData is to produce a single output file
        that contains all points that are inside all of the polygons in the
        shapefile. Optional behaviors include including only points outside the
        polygons, producing individual files containing the points within each
        polygon in the shapefile, and clipping points within a single polygon
        specified using a field from the shapefile database.

        Parameters
        ----------
        polyfile: string, path to file
            Name of the polygon file used for clipping. Format should be Arc
            shapefile.
        outputfile: string, path to file
            Base name for output data files. Default behavior is to create one
            output file for all polygons in polyfile. See below for use of
            /outside, /shape, and /multifile switches to modify this behavior.
        datafile: string, path to file
            LIDAR data file template, list of data files, or name of a text file
            containing a list of file names (must have .txt extension)

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        lda: boolean
            Store output data file in FUSION's LDA format
        index: boolean
            Create FUSION index files for output data file
        outside: boolean
            Create output file containing all points outside of all polygons in
            polyfile. When used with /shape switch, output file will contain all
            points outside the specified polygon.
        multifile: boolean
            Create separate output data files for each polygon in polyfile.
            NOTE: The /multifile switch can result in thousands of files for
            large polygon coverages.
        shape: 2-tuple or list-like (field#,value)
            Use the feature in polyfile identified by value in field field#.
            Output file will contain all points within the specified polygon.
            field# is a 1-based index that refers to fields in the DBASE file
            associated with the shapefile. The /shape switch is ignored for
            other formats. If polygon identifier contains a space, enclose the
            identifier in quotes. Use a * for value to force indexing of the
            polygon file and parse polygon identifiers from field#. If you do
            not use the /shape switch in conjunction with the /multifile switch,
            output files will be identified using the sequential number
            associated with the polygon rather than a value from a database
            field.
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        precision: 3-tuple or list-like (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'polyclipdata'
        params = [polyfile, outputfile, datafile]
        self.run(cmd, *params, **kwargs)

    def repairgriddtm(self, groundfiles, extraspace, **kwargs):
        """RepairGridDTM creates a new set of DTM tiles for datasets where the
        set of DTM tiles do not provide full coverage for the area.

        This occurs most often when DTMs are delivered in a raster format where
        the vendor has not considered the correct alignment of the cells in
        adjacent tiles. Basically, the DTMs have been treated as raster data
        where the entire cell is represented by the value for the cell instead
        of lattice data where the values represent points located at the center
        of each cell. When the individual files are interpreted as surfaces,
        there is a 1 cell gap between adjacent tiles.

        Parameters
        ----------
        groundfiles: string, path to file
            File specifier for ground surface models in PLANS DTM format. May be
            wildcard or text list file (extension .txt only). Normally this will
            specify all DTM files for a project area
        extraspace: numeric
            Distance added around all sides of a DTM tile.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        """
        cmd = 'repairgriddtm'
        params = [groundfiles, extraspace]
        self.run(cmd, *params, **kwargs)

    def returndensity(self, outputfile, cellsize, datafiles, **kwargs):
        """ReturnDensity produces raster outputs containing the number of
        returns in each cell.

        This type of output can also be produced using the Catalog utility but
        ReturnDensity is more efficient and only produces the return count
        layer. Output from ReturnDensity is useful in the LTKProcessor and
        AreaProcessor workflow tools to help subdivide large acquisitions into
        manageable chunks for processing. ReturnDensity also produces percentile
        values for first return intensity values. These are useful for
        establishing scaling values when creating intensity images or comparing
        intensity values for different acquistitions.

        Parameters
        ----------
        outputfile: string, path to file
            Name for output raster file (stored in PLANS DTM format with .dtm
            extension)
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data
        datafile: string, or list-like of strings, path to file(s)
            LIDAR data file (LDA, LAS, ASCII LIDARDAT formats)... may be
            wildcard or name of text file listing the data files. If wildcard or
            text file is used, no other datafile# parameters will be recognized.
            Several data files can be specified. The limit depends on the length
            of each file name. When using multiple data files, it is best to use
            a wildcard for datafile1 or create a text file containing a list of
            the data files and specifying the list file as datafile1.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        """
        cmd = 'returndensity'
        params = [outputfile, cellsize, datafiles]
        self.run(cmd, *params, **kwargs)

    def splitdtm(self, inputdtm, outputdtm, columns, rows, **kwargs):
        """SplitDTM divides a .DTM file into smaller tiles.

        It is most often used when bare-ground surface models are delivered in
        large tiles that cannot be managed efficiently. The command-line tools
        in FUSION all provide the ability to work with multiple ground surface
        files. In operation, the tools either build new, temporary grids
        covering the area of interest or access the ground models from disk.
        When creating workflows for processing data covering large areas,
        LTKProcessor allows users to specify the maximum number of returns in a
        processing tile. For datasets where the ground model tiles are large,
        some tools have problems loading the large models while trying to build
        a temporary model to support the analysis tasks. The result is that the
        programs start to page data in and out of memory and performance suffers
        by several orders of magnitude. By subdividing large ground tiles into
        smaller tiles, memory is used more efficiently and the paging behavior
        is avoided. SplitDTM offers two ways to specify the arrangement of the
        new tiles. In the first, the user specifies the number of rows and
        columns of new tiles and for the second the user specifies the maximum
        number of cells in a new tile. SplitDTM provides a default maximum
        number of cells (25 million) that produces tiles that are about 100Kb
        and cover an area that is 5000 by 5000 cells. This size seems to work
        well with all processing tools and improves the overall processing
        efficiency.

        Parameters (required)
        ----------
        inputdtm: string, path to file
            Name of the existing PLANS format DTM file to be subdivided.
        outputdtm: string, path to file
            Base name for the new PLANS format DTM tiles. The actual tile name
            will have the column and row appended to it.
        columns: int
            Number of columns of tiles.
        rows: int
            Number of rows of tiles.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        maxcells: int
            Maximum number of cells in the output tiles. Columns and Rows are
            ignored (but values are needed on the command line). New values for
            Columns and Rows will be calculated. The default maximum number of
            cells is 25,000,000.
        """
        cmd = 'splitdtm'
        params = [inputdtm, outputdtm, columns, rows]
        self.run(cmd, *params, **kwargs)

    def surfacesample(self, surfacefile, inputfile, outputfile, **kwargs):
        """SurfaceSample produces a comma separated values (CSV) file that
        contains a value interpolated from the surface at a specific XY
        location.

        Input is a file containing a list of XY coordinate pairs, one pair per
        line of the file. Output is a CSV file containing the original XY
        location and the surface value. Optionally SurfaceSample can generate a
        network of radial profiles given center point, the number of radial
        lines, length of each line and the desired point spacing. SurfaceSample
        can also generate a series of evenly spaced sample points along a line
        specified by two endpoints. The formats of the input and output files
        vary depending on the options used.

        Parameters
        ----------
        surfacefile: string, path to file
            Name for the input surface file (PLANS DTM format) surfacefile may
            be wildcard or text list file (extension .txt)
        inputfile: string, path to file
            Name of file containing XY locations (space or comma delimited). If
            the /pattern option is used with type 3, the inputfile should
            contain two coordinate pairs that specify the endpoint of the line.
            If the /id option is used, the inputfile should contain a point
            identifier in the first column.
        outputfile: string, path to file
            Name for the output ASCII CSV file

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        pattern: 4-tuple or list-like (type,p1,p2,p3)
            Generate a test pattern of sample points centered on the XY location
            from the inputfile.

            Pattern type 1 is a radial network of p1 lines that are at
            least p2 long with sample points every p3 units. The first radial
            is at 3 o'clock. Radials are generated in counter-clockwise order

            Pattern type 2 is a radial network of p1 lines that are at
            least p2 long with sample points every p3 units ON AVERAGE. The
            sample point spacing decreases as the distance from the XY
            location increases to provide sample point locations that
            represent a uniform area The first radial is at 3 o'clock.
            Radials are generated in counter-clockwise order.

            Pattern type 3 expects two coordinate pairs in the input data.
            The pairs specify the endpoints of a line. Lines with points
            spaced p1 units apart are created and written to the output.
        topo: 2-tuple or list-like of numerics (dist,lat)
            Compute solar radiation index (SRI) as defined in:
            Keating, K. A.,P. J. P. Gogan, J. M. Vore, and L. R. Irby. 2007.
            A simple solar radiation index for wildlife habitat studies. Journal
            of Wildlife Management 71(4):1344-1348.

            Algorithm uses a 3- by 3-cell grid with cells that are dist by dist
            units to compute topographic attributes. Latitude values are in
            degrees with north positive and south negative. Output from this
            option includes slope, aspect, SRI, and other topographic indices.
            The range for SRI is 0.0 to 2.0 (differs from the -1.0 to 1.0 range
            in Keating et al. 2007).
        noheader: boolean
            Suppress the header line in the output file. This option is useful
            when you want to use PDQ to view the outputfile
        novoid: boolean
            Suppress output for XY sample points ouside the surface extent or
            for points with invalid surface elevations
        id: boolean
            Read a point identifier from the first field of each line from the
            inputfile and output it as the first field of the outputfile. If id
            is used with pattern, a separate output file is created for each
            point in the inputfile. outputfiles are named using outputfile as
            the base name with the point identifier appended to the filename.
            Even when inputfile contains a single point, the outputfile name is
            modified to reflect the point identifier.
        """
        cmd = 'surfacesample'
        params = [surfacefile, inputfile, outputfile]
        self.run(cmd, *params, **kwargs)

    def thindata(self, outputfile, density, cellsize, datafile, **kwargs):
        """Thins LIDAR data to specific pulse densities.

        This capability is useful when comparing analysis results from several
        LIDAR acquisitions that were collected using different pulse densities.
        ThinData is also useful when the density within a single LIDAR data set
        is not uniform. This is often the case with data collected from a slow-
        flying helicopter or when flightline overlap was not closely monitored.
        ThinData has also been used in simulation experiments to assess the
        effect of LIDAR pulse density on the accuracy of estimated forest
        inventory metrics such as overall tree height.

        Parameters
        ----------
        outputfile: string, path to file
            Name for output data file...will contain new dataset thinned to the
            desired density
        density: numeric
            Desired data density in pulses per unit area
        cellsize: numeric
            Cell size used to compute data density specified in square units
        datafile: string, path to file
            LIDAR data file template or name of a text file containing a list of
            file names (must have .txt extension)

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            Suppress the use of the LASzip dll (c) Martin Isenburg...
            Removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        rseed: int (#)
            Use random number stream # (100 streams are available 0-99)
        index: boolean
            Create FUSION index files for the OutputFile
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3. All points will be used to determine the pulse arrangement
            but only points that match the specified classification values will
            be written to the output file.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        lda: boolean
            Write sample files using FUSION's LDA format when using LAS input
            files. The default behavior of ClipData (after version 2.35) is to
            write data in LAS format when the input data are in LAS format. When
            using input data in a format other than LAS, sample files are
            written in LDA format.
        precision: 3-tuple or list-like (scaleX,scaleY,scaleZ)
            Control the scale factor used for X, Y, and Z values in output LAS
            files. These values will override the values in the source LAS
            files. There is rarely any need for the scale parameters to be
            smaller than 0.001.
        """
        cmd = 'thindata'
        params = [outputfile, density, cellsize, datafile]
        self.run(cmd, *params, **kwargs)

    def tiledimagemap(self, outputhtml, indeximage, tiletemplate, **kwargs):
        """TiledImageMap creates a web page consisting of a single image map
        that corresponds to a mosaic of image tiles.

        The image map is linked to individual tile images allowing the user of
        the page to browse the coverage area switching between the larger
        overview image and the higher-resolution individual image tiles. For
        colored overview images, a legend image that describes the image can be
        included. TiledImageMap is most often used to organize intensity images
        created from LIDAR data but it can be used to provide web-ready display
        of any spatial information that is organized into tiles.

        TiledImageMap is particularly useful when LIDAR data have been delivered
        in “tiles” and subsequent data products have been produced using files
        representing individual “tiles”. Creating web pages with TiledImageMap
        makes it easy to browse analysis results and facilitates access to
        analysis products without using GIS.

        Parameters
        ----------
        outputhtml: string, path to file
            Name of output file name for HTML page (extension not needed)
        indeximage: string, path to file
            Name of large image used to create the image map. The image must
            have a corresponding world file.
        tiletemplate: string, path to file
            File template for image tiles or name of a text file containing a
            list of image tiles (must have .txt extension)

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific numeric
            formats (e.g. use a comma for the decimal separator)
        legend: string, path to file
            Add a legend file to the HTML page to the left of the image map file
            is the name of the image to use for the legend
        """
        cmd = 'tiledimagemap'
        params = [outputhtml, indeximage, tiletemplate]
        self.run(cmd, *params, **kwargs)

    def tinsurfacecreate(self, surfacefile, cellsize, xyunits, zunits,
                         coordsys, zone, horizdatum, vertdatum, datafiles,
                         **kwargs):
        """Creates a gridded surface model from point data.

        The algorithm used in TINSurfaceCreate first creates a TIN surface model
        using all of the points and then interpolates a gridded model from the
        TIN surface. TINSurfaceCreate works well when all points in a dataset
        are to be used as surface points. For example, after filtering a LIDAR
        point cloud to identify bare-ground points, TINSurfaceCreate can be used
        to create a gridded surface model. However, if any non-ground points
        remain in the dataset, they will be incorporated into the TIN ground
        surface model and will, most likely, affect the resulting gridded
        surface model. In general, TINSurfaceCreate should be used only when you
        know all points in the dataset are surface points. If this is not the
        case, use GridSurfaceCreate to create the gridded surface model.

        Parameters
        ----------
        surfacefile: string, path to file
            Name for output surface file (stored in PLANS DTM format with .dtm
            extension)
        cellsize: numeric
            Desired grid cell size in the same units as LIDAR data
        xyunits: string
            Units for LIDAR data XY
                'M' for meters
                'F' for feet.
        zunits: string
            Units for LIDAR data elevations:
                'M' for meters
                'F' for feet
        coordsys: int
            Coordinate system for the canopy model:
                0 for unknown
                1 for UTM
                2 for state plane
        zone: int
            Coordinate system zone for the canopy model (0 for unknown).
        horizdatum: int
            Horizontal datum for the canopy model:
                0 for unknown
                1 for NAD27
                2 for NAD83
        vertdatum: int
            Vertical datum for the canopy model:
                0 for unknown
                1 for NGVD29
                2 for NAVD88
                3 for GRS80
        datafiles: string or list-like of strings of paths to file(s)
            LIDAR data file (LDA, LAS, ASCII LIDARDAT formats)...may be wildcard
            or name of text file listing the data files. If wildcard or text
            file is used, no other datafiles will be recognized.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        return:string Specifies the returns to be included in the sample (can
              include A,1,2,3,4,5,6,7,8,9,F,L,O) Options are specified without
              commas (e.g. /return:123) For LAS files only: F indicates first
              and only returns, L indicates last of many returns.
        las_class: string or list-like
            LAS files only: Specifies that only points with classification
            values listed are to be included in the subsample. Classification
            values should be separated by a comma. e.g. (2,3,4,5) and can range
            from 0 to 31. If the first character in string is ~, the list is
            interpretted as the classes you DO NOT want included in the
            subsample. e.g. /class:~2,3 would include all class values EXCEPT 2
            and 3.
        ignoreoverlap: boolean
            Ignore points with the overlap flag set (LAS V1.4+ format)
        """
        cmd = 'tinsurfacecreate'
        params = [
            surfacefile, cellsize, xyunits, zunits, coordsys, zone, horizdatum,
            vertdatum, datafiles
        ]
        self.run(cmd, *params, **kwargs)

    def topometrics(self, surfacefile, cellsize, topopointspacing, latitude,
                    tpiwindowsize, outputfile, **kwargs):
        """Computes topographic metrics using surface models.

        The logic it uses is exactly the same as that used in GridMetrics except
        TopoMetrics computes a topographic position index (TPI) based on methods
        described by Weiss (2001) and Jenness (2006).

        Parameters
        ----------
        surfacefile: string, path to file
            Name for the input surface file (PLANS DTM format). surfacefile may
            be wildcard or text list file (extension .txt)
        cellsize: numeric
            Size of the cell used to report topographic metrics
        topopointspacing: numeric
            The spacing between points in the 3 by 3 array used to compute the
            basic topographic metrics.
        latitude: numeric
            Latitude for the center of the data area. North latitude is
            positive, South is negative. Latitude is used to compute the solar
            radiation index.
        tpiwindowsize: numeric
            The size of the window used to compute the topographic position
            index (TPI). When the /square option is used, the TPI window will be
            tpiwindowsize by tpiwindowsize. For round windows, the diameter will
            be TPIWindowSize.
        outputfile: string, path to file
            Base name for the output metrics (CSV format). "_topo_metrics" will
            be appended to the name provided on the command line.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        grid: 4-tuple or list-like (X,Y,W,H)
            Force the origin of the output grid to be (X,Y) instead of computing
            an origin from the data extents and force the grid to be W units
            wide and H units high...W and H will be rounded up to a multiple of
            cellsize
        gridxy: 4-tuple or list-like (X1,Y1,X2,Y2)
            Force the origin of the output grid to be (X1,Y1) instead of
            computing an origin from the data extents and force the grid to use
            (X2,Y2) as the upper right corner of the coverage area. The actual
            upper right corner will be adjusted to be a multiple of cellsize
        align: string, path to file
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file
        extent: string, path to file
            Force the origin and extent of the output grid to match the lower
            left corner and extent of the specified PLANS format DTM file but
            adjust the origin to be an even multiple of the cell size and the
            width and height to be multiples of the cell size.
        square: boolean
            Use a square-shaped mask when computing the topographic position
            index. The default mask shape is a circle.
        annulusinnerdia: numeric
            Use a donut-shaped mask when computing the topographic position
            index. The outer diameter is topodistance and the inner diameter is
            dia. When using this option to define a narrow ring, use the
            /verbose option to display the computed mask to ensure it meets your
            expectations.
        annuluswidth: numeric
            Use a donut-shaped mask when computing the topographic position
            index. The outer diameter is tpiwindowsize and the inner diameter is
            tpiwindowsize-(width*2). When using this option to define a narrow
            ring, use the /verbose option to display the computed mask to ensure
            it meets your expectations.
        diskground: boolean
            Do not load ground surface models into memory or create a temporary
            surface in memory. When this option is specified, larger areas can
            be processed but processing will be very much slower.
        nointernalground: boolean
            Do not create a temporary surface in memory. When this option is
            specified, larger areas can be processed but processing will be much
            slower. This option has no effect when used with /diskground.
        lockdtmcellsize: boolean
            Force the cell size used when creating an internal ground model to
            be the same as the highest resolution input ground surface model.
            The default behavior will increase the cell size until the internal
            model will fit in the available memory.
        """
        cmd = 'topometrics'
        params = [
            surfacefile, cellsize, topopointspacing, latitude, tpiwindowsize,
            outputfile
        ]
        self.run(cmd, *params, **kwargs)

    def treeseg(self, chm, ht_threshold, outputfile, **kwargs):
        """The TreeSeg program applies a watershed segmentation algorithm to a
        canopy height model to produce “basins” that correspond to dominate
        clumps of tree foliage and branches.

        In some cases, the segments represent individual trees but it is common
        for segments to encompass several tree crowns. The resulting segments
        also represent dominant and co-dominant trees better that mid-story and
        under-story/suppressed trees. Output can consist of several products:
            * Basin list in CSV format containing the basin number (first basin
              is 2), the location of the high point, number of canopy height
              model cells within the basin, maximum surface height for the
              basin, and the row and column within the canopy height model for
              the basin high point.
            * Basin map in ASCII raster format: Each basin is assigned a unique
              number starting with 2. Areas not part of any basin are given a
              value of 1.
            * Maximum basin height map in ASCII raster format. Each basin is
              assigned a value corresponding to the maximum height for the
              basin.
            * High point list in shapefile format with the same fields as the
              basin list.
            * Basin/crown perimeter polygons in shapefile format. Polygon
              attributes include the location of the highest point (from the
              basin), the actual area of the polygon, and the maximum height of
              the basin.

        Parameters
        ----------
        chm: string, path to file
            Name for canopy height model (CHM), (PLANS DTM with .dtm extension).
            May be wildcard or text list file (extension .txt only). This can be
            a canopy surface model if the /ground option is used to specify a
            ground surface for normalization.
        ht_threshold: numeric
            Minimum height for object segmentation. Portions of the CHM below
            this height are not considered in the segmentation.
        outputfile: string, path to file
            Base name for output file. Metrics are stored in CSV format with
            .csv extension. Other outputs are stored in files named using the
            base name and additional descriptive information.

        **Kwargs (optional), aka "Switches" in FUSION
        --------
        interactive: boolean
            Present a dialog-based interface
        quiet: boolean
            Suppress all output during the run
        verbose: boolean
            Display all status information during the run
        version: boolean
            Report version information and exit with no processing
        newlog: boolean
            Erase the existing log file and start a new log
        log: string, path to file
            Use the name specified for the log file
        locale: boolean
            Adjust program logic to input and output locale-specific
            numericformats (e.g. use a comma for the decimal separator)
        nolaszipdll: boolean
            suppress the use of the LASzip dll (c) Martin Isenburg...
            removes support for compressed LAS (LAZ) files. This option is only
            useful for programs that read or write point files.
        height: boolean
            Normalize height model(s) using ground model(s)
        ptheight: boolean
            Normalize point heights using ground model(s)
        maxht: numeric
            Force the maximum height for the segmentation. This will override
            the actual maximum value in the CHM. Use this option to force equal
            vertical resolution across areas with varying maximum canopy
            heights.
        grid: 4-tuple or list-like (X,Y,W,H)
            Force the origin of the analysis area to be (X,Y) instead of
            computing an origin from the CHM extent and force the width and
            height to be W and H.
        gridxy: 4-tuple or list-like (X1,Y1,X2,Y2)
            Force the origin of the analysis area to be (X1,Y1) instead of
            computing an origin from the CHM extent and force the area to use
            (X2,Y2) as the upper right corner.
        align: string, path to file
            Force the origin and extent of the analysis area to match the lower
            left corner and extent of the specified PLANS format DTM file
        buffer: numeric
            Add a buffer to the data extent specified by /grid, /gridxy or
            /align when segmenting but only output data for the segments located
            within the extent.
        ground: string, path to file
            Use a surface file to normalize the canopy surface (PLANS DTM with
            .dtm extension). May be wildcard or text list file (extension .txt
            only)
        points: string, path to file
            LIDAR point data file(s) in LDA or LAS format. May be wildcard or
            text list file (extension .txt only). Points are assigned to
            individual basins or crown polygons and a separate file (in LDA
            format) is output for each basin or polygon.
        segmentpts: boolean
            Output points for the raster segments. Default is to output points
            for crown polygons when the /shape option is used and for raster
            segments when /shape is not used. Used only with the /points option.
        clipfolder: string, path to folder
            Folder name where point files for individual clips are stored. Used
            only with the /points option. If not specified, point files are
            stored in the same folder with other outputs. The folder name must
            end with a trailing backslash and must already exist.
        shape: boolean
            Create a shapefile containing the high points and basin metrics
        cleantile: boolean
            Output an ASCII raster map that only includes basins within the
            reporting extent defined by the /grid, /gridxy, and /align options.
        htmultipler: numeric (#)
            Multiply the high point and surface heights by # for output
            products. Also multiply individual point heights by # before writing
            point files (see /points option).
        projection: string, path to file
            Associate the specified projection file with shapefile and raster
            data products.
        """
        cmd = 'treeseg'
        params = [hm, ht_threshold, outputfile]
        self.run(cmd, *params, **kwargs)
