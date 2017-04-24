# extractETAS

## MATLAB function:

    [x, y, z, t, data] = extractETAS(xs, ys, zs, date_start, date_end, variable)
  
  Extract ETAS ocean model (Oliver et al. 2016) data over specified spatial and temporal domains. [Daily](http://data.imas.utas.edu.au/portal/search?uuid=47450e8b-da02-463d-aca4-f00b2077bd56) and [monthly](http://data.imas.utas.edu.au/portal/search?uuid=43dee465-0bc1-4f5a-b855-4cab14eaa9f1) ETAS data are available on the [IMAS data portal](http://data.imas.utas.edu.au/) and [THREDDS server]``(http://thredds.imas.utas.edu.au/thredds/catalog/IMAS/catalog.html).

###  Inputs:

    xs         - Longitude(s) of requested data [degrees East]
    ys         - Latitude(s) of requested data [degrees North, so 42S = -42]
    zs         - Depth(s) of requested data, [metres, positive downward]
    date_start - Start date of requested data [year month day] (earliest: [1993 1 1])
    date_end   - End date of requested data [year month day] (latest: [2016 12 31])
    variable   - Which model variable is requested

   xs, ys, zs can each be either a scalar value e.g. ys = -42.5
   or a range e.g. xs = [148.5 149] or zs = [30 100]. If xs and ys
   are both scalars, then the data is provided at the nearest neighbour
   model grid cell, if xs or ys are provided as a range [min max] then
   the data are interpolated onto a regular 0.01 x 0.01 degree grid and
   and provided as a rectangular matrix. If zs is a scalar then then
   the data are provided on the nearest vertical level, if zs is a
   range then the data are provided over the range of model levels
   that encompass the requested range.

###   Possible variables are

    'temp'    - Temperature [deg. C]
    'salt'    - Salinity [PSU]
    'u'       - Zonal velocity [m/s]
    'v'       - Meridional velocity [m/s]
    'w'       - Vertical velocity [m/s]
    'dens'    - Density [kg/m^3]
    'eta'     - Sea surface height [m]
    'nhf'     - Net surface heat flux [W/m^2]
    'botz'    - Water depth [m]
    'u1'      - Velocity along model coordiate x1 [m/s]
                (approximately along-shore, positive east/northwards)
    'u2'      - Velocity along model coordiate x2 [m/s]
                (approximately cross-shore, positive towards land)

  All variables are 4-dimensional (x,y,z,t) except 'eta' and 'nhf'
  which are 3-D (x,y,t) and 'botz' which is 2-D (x,y). Therefore,
  the vertical coordinate (zs) and date ranges (date_start, end)
  are ignored if not required.

###  Optional inputs:

    [x, y, z, t, data] = extractETAS(xs, ys, zs, date_start, date_end, variable, dt, resolution)

    dt         - Time resolution ('daily', 'monthly', 'annual')
                 Output provides data values with that sample rate. 'daily'
                 provides daily snapshots while 'monthly' and 'annual' are
                 averages over that time scale.
                 For 'monthly', the specified start and end days are ignored
                 (i.e. whole months are used) and for 'annual' the months are
                 also ignored (whole years are used)
                 [DEFAULT = 'daily'; any value other than 'daily', 'monthly'
                 or 'annual' provides daily means]
    resolution - Horizontal grid resolution [degrees lat/lon] onto which the
                 data is interpolated. Only relevant if xs and/or ys are a
                 range rather than single point location
                 [DEFAULT = 0.01 degrees; this is also the finest resolution
                 permitted]

   Note that extraction of daily data over a large area or over many years
   can be slow due to the interpolation calculations. The interpolation
   can also be memory intensive so ensure you have free memory or choose
   a coarser resolution onto which to interpolate the data. Calculations
   on 'monthly' and 'annual' data are significantly faster and less 
   memory-intensive than on 'daily' data.

   The data is extracted from the online IMAS THREDDS server*, and so an
   internet connection is required. If you have the data files locally, then
   change the values of the 'header' variables in the extractETAS.m script
   to point to the data location.
   http://thredds.imas.utas.edu.au/thredds/catalog/IMAS/catalog.html

###  Outputs:

    'x', 'y', 'z'  - The spatial coordinates of the output data
    't'            - The time coordinates of the output data (datenum format)
    'data'         - The requested data, an array of size (X, Y, Z, T) where
                     the dimensions are for the x, y, z, t dimensions resp.

###  Examples:

   To extract a daily time series of temperature at [148.5E, 41S], 5 m depth
   between 3 March 1995 and 20 September 1998:

    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp');

   To extract the same data as monthly means:

    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp', 'monthly');

   or as annual means:

    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp', 'annual');

   To extract the daily data as above over an (x,z) section in the range [148,149]E and [0 100]m

    [x, y, z, t, data] = extractETAS([148 149], -41, [0 100], [1995 3 3], [1998 9 20], 'temp');

   To extract monthly data over a rectangular area, at one depth, interpolated to 0.02 degrees
   grid instead of the default 0.01 degree grid:

    [x, y, z, t, data] = extractETAS([148 149], [-43 -42], 5, [1995 3 3], [1998 9 20], 'temp', 'monthly', 0.02);

###  References:

   Oliver, E. C. J., M. Herzfeld and N. J. Holbrook (2016), Modelling the shelf
   circulation off eastern Tasmania, Continental Shelf Research, 130, pp. 14-33,
   doi: 10.1016/j.csr.2016.10.005. url:
   http://passage.phys.ocean.dal.ca/~olivere/docs/OliverHerzfeldHolbrook_2016_CSR_ETAS.pdf
