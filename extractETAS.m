function [x, y, z, t, data] = extractETAS(varargin)

%
%  [x, y, z, t, data] = extractETAS(xs, ys, zs, date_start, date_end, variable)
%  
%  Extract ETAS ocean model (Oliver et al. 2016) data over specified
%  spatial and temporal domains.
%
%  Inputs:
%
%   xs         - Longitude(s) of requested data [degrees East]
%   ys         - Latitude(s) of requested data [degrees North, so 42S = -42]
%   zs         - Depth(s) of requested data, [metres, positive downward]
%   date_start - Start date of requested data [year month day] (earliest: [1993 1 1])
%   date_end   - End date of requested data [year month day] (latest: [2016 12 31])
%   variable   - Which model variable is requested
%
%   xs, ys, zs can each be either a scalar value e.g. ys = -42.5
%   or a range e.g. xs = [148.5 149] or zs = [30 100]. If xs and ys
%   are both scalars, then the data is provided at the nearest neighbour
%   model grid cell, if xs or ys are provided as a range [min max] then
%   the data are interpolated onto a regular 0.01 x 0.01 degree grid and
%   and provided as a rectangular matrix. If zs is a scalar then then
%   the data are provided on the nearest vertical level, if zs is a
%   range then the data are provided over the range of model levels
%   that encompass the requested range.
%
%   Possible variables are
%
%    'temp'    - Temperature [deg. C]
%    'salt'    - Salinity [PSU]
%    'u'       - Zonal velocity [m/s]
%    'v'       - Meridional velocity [m/s]
%    'w'       - Vertical velocity [m/s]
%    'dens'    - Density [kg/m^3]
%    'eta'     - Sea surface height [m]
%    'nhf'     - Net surface heat flux [W/m^2]
%    'botz'    - Water depth [m]
%    'u1'      - Velocity along model coordiate x1 [m/s]
%                (approximately along-shore, positive east/northwards)
%    'u2'      - Velocity along model coordiate x2 [m/s]
%                (approximately cross-shore, positive towards land)
%
%    All variables are 4-dimensional (x,y,z,t) except 'eta' and 'nhf'
%    which are 3-D (x,y,t) and 'botz' which is 2-D (x,y). Therefore,
%    the vertical coordinate (zs) and date ranges (date_start, end)
%    are ignored if not required.
%
%  Optional inputs:
%
%  [x, y, z, t, data] = extractETAS(xs, ys, zs, date_start, date_end, variable, dt, resolution)
%
%   dt         - Time resolution ('daily', 'monthly', 'annual')
%                Output provides data values with that sample rate. 'daily'
%                provides daily snapshots while 'monthly' and 'annual' are
%                averages over that time scale.
%                For 'monthly', the specified start and end days are ignored
%                (i.e. whole months are used) and for 'annual' the months are
%                also ignored (whole years are used)
%                [DEFAULT = 'daily'; any value other than 'daily', 'monthly'
%                or 'annual' provides daily means]
%   resolution - Horizontal grid resolution [degrees lat/lon] onto which the
%                data is interpolated. Only relevant if xs and/or ys are a
%                range rather than single point location
%                [DEFAULT = 0.01 degrees; this is also the finest resolution
%                permitted]
%
%   Note that extraction of daily data over a large area or over many years
%   can be slow due to the interpolation calculations. The interpolation
%   can also be memory intensive so ensure you have free memory or choose
%   a coarser resolution onto which to interpolate the data. Calculations
%   on 'monthly' and 'annual' data are significantly faster and less 
%   memory-intensive than on 'daily' data.
%
%   The data is extracted from the online IMAS THREDDS server*, and so an
%   internet connection is required. If you have the data files locally, then
%   change the values of the 'header' variables in the extractETAS.m script
%   to point to the data location.
%   * http://thredds.imas.utas.edu.au/thredds/catalog/IMAS/catalog.html
%
%  Outputs:
%
%   'x', 'y', 'z'  - The spatial coordinates of the output data
%   't'            - The time coordinates of the output data (datenum format)
%   'data'         - The requested data, an array of size (X, Y, Z, T) where
%                    the dimensions are for the x, y, z, t dimensions resp.
%
%  Examples:
%
%   To extract a daily time series of temperature at [148.5E, 41S], 5 m depth
%   between 3 March 1995 and 20 September 1998:
%
%    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp');
%
%   To extract the same data as monthly means:
%
%    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp', 'monthly');
%
%   or as annual means:
%
%    [x, y, z, t, data] = extractETAS(148.5, -41, 5, [1995 3 3], [1998 9 20], 'temp', 'annual');
%
%   To extract the daily data as above over an (x,z) section in the range [148,149]E and [0 100]m
%
%    [x, y, z, t, data] = extractETAS([148 149], -41, [0 100], [1995 3 3], [1998 9 20], 'temp');
%
%   To extract monthly data over a rectangular area, at one depth, interpolated to 0.02 degrees
%   grid instead of the default 0.01 degree grid:
%
%    [x, y, z, t, data] = extractETAS([148 149], [-43 -42], 5, [1995 3 3], [1998 9 20], 'temp', 'monthly', 0.02);
%
%  References:
%
%   Oliver, E. C. J., M. Herzfeld and N. J. Holbrook (2016), Modelling the shelf
%   circulation off eastern Tasmania, Continental Shelf Research, 130, pp. 14-33,
%   doi: 10.1016/j.csr.2016.10.005. url:
%   http://passage.phys.ocean.dal.ca/~olivere/docs/OliverHerzfeldHolbrook_2016_CSR_ETAS.pdf
%

% grab inputs, nargs = 5 then output spatial average, otherwise output on grid
[cax, args, nargs] = axescheck(varargin{:});
xs = args{1};
ys = args{2};
zs = args{3};
date_start = args{4};
date_end = args{5};
var = args{6};
res = -1;
if nargs > 6
  % Set temporal mean time scale
  if sum(strcmp('daily', args(7:end))) > 0
    dt = 'daily';
  elseif sum(strcmp('monthly', args(7:end))) > 0
    dt = 'monthly';
  elseif sum(strcmp('annual', args(7:end))) > 0
    dt = 'annual';
  else
    dt = 'daily';
  end
  % Set horizontal resolution
  for i = [7:length(args)]
    if isnumeric(args{i})
      res = args{i};
    end
  end
  % Default if left unset
  if res == -1
    res = 0.01;
  end
else % Defaults
  res = 0.01;
  dt = 'daily';
end

% Point location or extended in (x,y,z)
Nx = length(xs);
Ny = length(ys);
Nz = length(zs);

if (Nx < 1) | (Nx > 2)
    error('Longitudes (xs) must be scalar, or a range indicated by an array of length two')
elseif (Nx == 2) & (xs(1) > xs(2))
    error('Longitude range (xs) must be [East West] i.e. West > East')
end
if (Ny < 1) | (Ny > 2)
    error('Latitudes (zs) must be scalar, or a range indicated by an array of length two')
elseif (Ny == 2) & (ys(1) > ys(2))
    error('Latitude range (ys) must be [South North] i.e. North > South')
end
if (Nz < 1) | (Nz > 2)
    error('Depths (zs) must be scalar, or a range indicated by an array of length two')
elseif (Nz == 2) & (zs(1) > zs(2))
    error('Depth range (zs) must be [Shallow Deep] i.e. Deep > Shallow (depth is defined positive donwward)')
end
if res < 0.01
    error('Resolution (res) must be at least 0.01 degrees')
end

% Start and end times
ts = datenum(date_start);
te = datenum(date_end);
t = ts:te; % time vector
% Filter by available time
t_ETAS_start = datenum([1993 1 1]);
t_ETAS_end = datenum([2016 12 31]);
t_ETAS = t_ETAS_start:t_ETAS_end;
t = t(ismember(t, t_ETAS));
% Monthly or annual time index
if strcmp(dt, 'monthly')
  dates = datevec(t);
  years = unique(dates(:,1));
  t0 = [];
  for year = years'
    months = unique(dates(dates(:,1)==year,2));
    for month = months'
      t0 = [t0; datenum([year month 15])]; % Monthly mean centred mid-month
    end
  end
  t = t0;
elseif strcmp(dt, 'annual')
  dates = datevec(t);
  years = unique(dates(:,1));
  t0 = [];
  for year = years'
    for month = [1:12] % All months required for annual mean
      t0 = [t0; datenum([year month 15])];
    end
  end
  t = t0;
end
T = length(t);
%
dates = datevec(t);
years = unique(dates(:,1));
t_SHOC = t - datenum([1990 1 1]);

% File locations
if strcmp(dt, 'daily')
  %header = '/media/ecoliver/Insect/data/SHOC/etas/runs/historical/rivers/no_tides/combined/etas_';
  header = 'http://thredds.imas.utas.edu.au/thredds/dodsC/IMAS/EOliver_ETAS_V2_ocean_model_output_daily/etas_';
  footer = '_combined.nc';
else
  %header = '/media/ecoliver/Insect/data/SHOC/etas/runs/historical/rivers/no_tides/combined/monthly/etas_';
  header = 'http://thredds.imas.utas.edu.au/thredds/dodsC/IMAS/EOliver_ETAS_V2_ocean_model_output_monthly/etas_';
  footer = '_combined_monthly.nc';
end
% Load in some basic grif information
% Output variables
if strcmp(var, 'u1') % Defined on left grid face
  xc = ncread([header '2014' footer], 'x_left');
  yc = ncread([header '2014' footer], 'y_left');
elseif strcmp(var, 'u2') % Defined on back grid face
  xc = ncread([header '2014' footer], 'x_back');
  yc = ncread([header '2014' footer], 'y_back');
else % All other variables defined on grid centre
  xc = ncread([header '2014' footer], 'x_centre');
  yc = ncread([header '2014' footer], 'y_centre');
end
zc = ncread([header '2014' footer], 'z_centre');
z_edges = ncread([header '2014' footer], 'z_grid'); z_edges(end) = 0;
dz = diff(z_edges);

% Initialize data variable
% xy-grid
if (Nx == 1) & (Ny == 1) % Point location
  X = 1; X1 = 1;
  Y = 1; X2 = 1;
  d = (xc - xs).^2 + (yc - ys).^2;
  [M, I] = min(d(:));
  [i1, i2] = ind2sub(size(d), I);
  x = xc(i1, i2);
  y = yc(i1, i2);
  pointLoc = 1;
else % Fine grid
  x = xs(1):res:xs(end);
  y = ys(1):res:ys(end);
  X = length(x);
  Y = length(y);
  i1 = 1;
  i2 = 1;
  [X1, X2] = size(xc);
  pointLoc = 0;
  % Full grid
  [xx, yy] = meshgrid(x, y);
  xc_vec = xc(:);
  yc_vec = yc(:);
  valid = ~isnan(xc_vec);
  xc_vec = xc_vec(valid);
  yc_vec = yc_vec(valid);
end

% z-grid
if Nz == 1 % Single level
  Z = 1;
  [M, k1] = min(abs(zc - (-1*zs)));
  z = abs(zc(k1));
else % Vertical range
  kk = find((zc < -zs(1)) & (zc > -zs(end)));
  k1 = max(kk(1) - 1, 1);
  k2 = min(kk(end) + 1, length(zc));
  z = abs(zc(k1:k2));
  Z = length(z);
end

% Output variables
if sum(strcmp(var, {'temp', 'salt', 'u', 'v', 'w', 'dens', 'u1', 'u2'})) > 0
  DIM = 4;
  data = zeros(X,Y,Z,T);
elseif sum(strcmp(var, {'eta', 'nhf'})) > 0
  DIM = 3;
  data = zeros(X,Y,T);
else % botz
  DIM = 2;
end

% Output status message
if pointLoc
  disp(['Extracting ' var ' at lon = ' num2str(xs) ', lat = ' num2str(ys) ', z = ' num2str(zs) ', ' datestr(t(1)) ' to ' datestr(t(end)) ', as ' dt ' means'])
else
  disp(['Extracting ' var ' over lon = [' num2str(xs(1)) ',' num2str(xs(end)) '], lat = [' num2str(ys(1)) ',' num2str(ys(end)) '], z = [' num2str(zs(1)) ',' num2str(zs(end)) '], ' datestr(t(1)) ' to ' datestr(t(end)) ', as ' dt ' means, and interpolating to ' num2str(res) ' degree horizontal resolution'])
end

% Loop over time to load in data
if DIM == 2 % 'botz'
    file = [header num2str(years(1)) footer];
    tmp = ncread(file, var, [i1 i2], [X1 X2]);
    tmp(tmp>0) = NaN;
    tmp = abs(tmp);
    if pointLoc
      data = tmp;
    else
      tmp3 = tmp(:);
      F = TriScatteredInterp(xc_vec, yc_vec, tmp3(valid));
      data = F(xx, yy)';
    end
    t = NaN;
else
  tt0 = 1;
  for year = years'
    disp(['Processing... ' num2str(year)])
    file = [header num2str(year) footer];
    t_year = ncread(file, 't');
    tt_which = find((t_year >= t_SHOC(1)) & (t_year <= t_SHOC(end)));
    tt1 = tt_which(1);
    TT = length(tt_which);
    % Load in data
    if DIM == 4 % 'temp', 'salt', 'u', 'v', 'w', 'dens'
      tmp = ncread(file, var, [i1 i2 k1 tt1], [X1 X2 Z TT]);
    elseif DIM == 3 % 'eta', 'nhf'
      tmp = ncread(file, var, [i1 i2 tt1], [X1 X2 TT]);
    end
    % Grid data or just deal with point location
    if pointLoc
      if DIM == 4 % 'temp', 'salt', 'u', 'v', 'w', 'dens'
        data(:,:,:,tt0:(tt0+TT-1)) = tmp;
      elseif DIM == 3 % 'eta', 'nhf'
        data(:,:,tt0:(tt0+TT-1)) = tmp;
      end
    else % Interpolate onto grid
      if DIM == 4
        for tt = 1:size(tmp,DIM)
          for k = 1:size(tmp,DIM-1)
            tmp2 = tmp(:,:,k,tt);
            tmp3 = tmp2(:);
            F = TriScatteredInterp(xc_vec, yc_vec, tmp3(valid));
            data(:,:,k,tt0+tt-1) = F(xx, yy)';
          end
        end
      elseif DIM == 3
        for tt = 1:size(tmp,DIM)
            tmp2 = tmp(:,:,tt);
            tmp3 = tmp2(:);
            F = TriScatteredInterp(xc_vec, yc_vec, tmp3(valid));
            data(:,:,tt0+tt-1) = F(xx, yy)';
        end
      end
    end
    tt0 = tt0 + TT;
  end
end

% Calculate annual means if necessary
if strcmp(dt, 'annual') & DIM > 2
  % New time vector
  dates = datevec(t);
  years = unique(dates(:,1));
  t0 = [];
  for year = years'
    t0 = [t0; datenum([year 7 1])]; % Centre annual means mid-year
  end
  t = t0;
  % Annual means of monthly data
  if DIM == 4;
    data0 = zeros(X,Y,Z,length(years));
    for year = years'
      tt = (dates(:,1) == year);
      data0(:,:,:,year-years(1)+1) = mean(data(:,:,:,tt), 4);
    end
  elseif DIM == 3;
    data0 = zeros(X,Y,length(years));
    for year = years'
      tt = (dates(:,1) == year);
      data0(:,:,year-years(1)+1) = mean(data(:,:,tt), 3);
    end
  end
  data = data0;
end

% Vertical coordinate for non-z data
if DIM == 4
  z = z;
elseif DIM == 3 % 'eta', 'nhf'
  z = 0;
else % botz
  z = NaN;
end

% Squeeze data matrix
data = squeeze(data);

