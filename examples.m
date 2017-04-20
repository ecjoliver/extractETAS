%
% Script of examples of loading and plotting some ETAS data
%

% To extract a time series of temperature at [148.5E, 43S], 20 m depth between 3 March 1995 and 20 September 1998:
figure; clf; hold on;
% The following extracts daily data (since dt is not specified, 'daily' is the default)
[x, y, z, t, data] = extractETAS(148.5, -43, 20, [1995 3 3], [1998 9 20], 'temp');
plot(t, data, 'k.')
% The following extracts daily data (since dt is explicitly specified as 'daily')
[x, y, z, t, data] = extractETAS(148.5, -43, 20, [1995 3 3], [1998 9 20], 'temp', 'daily');
plot(t, data, 'k-')
% The following extracts monthly data (since dt is explicitly specified as 'monthly')
[x, y, z, t, data] = extractETAS(148.5, -43, 20, [1995 3 3], [1998 9 20], 'temp', 'monthly');
plot(t, data, 'b-o')
% The following extracts annual data (since dt is explicitly specified as 'annual')
[x, y, z, t, data] = extractETAS(148.5, -43, 20, [1995 3 3], [1998 9 20], 'temp', 'annual');
plot(t, data, 'r-o')
ylabel('Temperature [deg. C]');
datetick('x', 'mmmyy')

% Extract at the same location and time range, but across a range of depths, and plot them all together
[x, y, z, t, data] = extractETAS(148.5, -43, [0 100], [1995 3 3], [1998 9 20], 'temp');
figure; clf; hold on;
contourf(t, -z, data); shading flat; colorbar;
ylabel('Depth [m]');
datetick('x', 'mmmyy')

% Extract daily data as above over an (x,z) section in the range [148,149]E and [0 100]m
% Note that this can be quite slow depending on your internet connectiona and computer speed...
[x, y, z, t, data] = extractETAS([148 149], -41, [0 100], [1995 3 3], [1998 9 20], 'temp');
figure; clf;
subplot(2,2,1); hold on;
contourf(x, -z, data(:,:,1)'); shading flat; colorbar;
title(datestr(t(1)))
subplot(2,2,2); hold on;
contourf(x, -z, data(:,:,400)'); shading flat; colorbar;
title(datestr(t(400)))
subplot(2,2,3); hold on;
contourf(x, -z, data(:,:,800)'); shading flat; colorbar;
title(datestr(t(800)))
subplot(2,2,4); hold on;
contourf(x, -z, data(:,:,1200)'); shading flat; colorbar;
title(datestr(t(1200)))

% Extract monthly data over the same domain/time range as above
% Extracting monthly data is significantly faster than daily...
% This time let's plot salinity ('salt')
[x, y, z, t, data] = extractETAS([148 149], -41, [0 100], [1995 3 3], [1998 9 20], 'salt', 'monthly');
figure; clf;
subplot(2,2,1); hold on;
contourf(x, -z, data(:,:,1)'); shading flat; colorbar;
title(datestr(t(1)))
subplot(2,2,2); hold on;
contourf(x, -z, data(:,:,10)'); shading flat; colorbar;
title(datestr(t(10)))
subplot(2,2,3); hold on;
contourf(x, -z, data(:,:,20)'); shading flat; colorbar;
title(datestr(t(20)))
subplot(2,2,4); hold on;
contourf(x, -z, data(:,:,30)'); shading flat; colorbar;
title(datestr(t(30)))

% Etract monthly data over a rectangular area, at one depth, interpolated to 0.02 degrees grid instead of the default 0.01 degree grid. Get temp and (u,v)
[x, y, z, t, temp] = extractETAS([148 149], [-43 -42], 5, [1995 3 3], [1998 9 20], 'temp', 'monthly', 0.02);
[x, y, z, t, u] = extractETAS([148 149], [-43 -42], 5, [1995 3 3], [1998 9 20], 'u', 'monthly', 0.02);
[x, y, z, t, v] = extractETAS([148 149], [-43 -42], 5, [1995 3 3], [1998 9 20], 'v', 'monthly', 0.02);
figure; clf; d = 4;
subplot(2,2,1); hold on;
contourf(x, y, temp(:,:,1)'); colorbar;
quiver(x(1:d:end), y(1:d:end), u(1:d:end,1:d:end,1)', v(1:d:end,1:d:end,1)', 'k');
title(datestr(t(1)))
subplot(2,2,2); hold on;
contourf(x, y, temp(:,:,10)'); colorbar;
quiver(x(1:d:end), y(1:d:end), u(1:d:end,1:d:end,10)', v(1:d:end,1:d:end,10)', 'k');
title(datestr(t(10)))
subplot(2,2,3); hold on;
contourf(x, y, temp(:,:,20)'); colorbar;
quiver(x(1:d:end), y(1:d:end), u(1:d:end,1:d:end,20)', v(1:d:end,1:d:end,20)', 'k');
title(datestr(t(20)))
subplot(2,2,4); hold on;
contourf(x, y, temp(:,:,30)'); colorbar;
quiver(x(1:d:end), y(1:d:end), u(1:d:end,1:d:end,30)', v(1:d:end,1:d:end,30)', 'k');
title(datestr(t(30)))

% Plot the water depth over the large portion of the domain
[x, y, z, t, data] = extractETAS([146.5 149.5], [-44.5 -40.5], 30, [2005 1 1], [2006 7 1], 'botz');
pcolor(x, y, data'); shading flat; colorbar;
