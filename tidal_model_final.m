%% 1. Setup and Data Loading
% Remove existing output files to avoid conflicts
if exist('tidal_animation.mp4','file'), delete('tidal_animation.mp4'); end
if exist('tidal_animation.avi','file'), delete('tidal_animation.avi'); end
if exist('tidal_data.mat','file'), delete('tidal_data.mat'); end

% clearing and cleaning previous loaded data and variables 
clear; close all; clc;

% 1.1 Loading the shapefile for model domain
model_domain = 'mohesh.shp';  % path to the shapefile
coast = shaperead(model_domain);
coastX = [coast.X]; coastY = [coast.Y];

% 1.2 Loading bathymetry data from CSV as lon, lat, elevation
bathy = readtable('bathymetry.csv');
lon = bathy.lon;
lat = bathy.lat;
z = bathy.elevation;

% 1.3 Creating spatial grid
dx = 0.001;  % degrees (~102 m at chittagong )  https://en.wikipedia.org/wiki/Decimal_degrees 
[xGrid, yGrid] = meshgrid(min(lon):dx:max(lon), min(lat):dx:max(lat));
zGrid = griddata(lon, lat, z, xGrid, yGrid, 'natural');

%% 2. Definig Tidal Constituents and Their Respective Parameters
% 2.1 Parameters for tidal constituents(M2, S2, N2, K2, K1, O1, and P1): amplitudes (m), phases (deg), and periods (h)
tidal_constituents = struct( ...
    'M2', struct('amplitude',0.88618882,'phase',87.530685,'period',12.42), ...
    'S2', struct('amplitude',0.4187011,'phase',117.06008,'period',12.00),  ...
    'N2', struct('amplitude',0.17735326,'phase',79.006669,'period',12.66), ...
    'K2', struct('amplitude',0.11927057,'phase',118.57742,'period',12.00), ...
    'K1', struct('amplitude',0.14572086,'phase',245.62559,'period',23.93), ...
    'O1', struct('amplitude',0.057382711,'phase',236.10362,'period',25.82), ...
    'P1', struct('amplitude',0.043141832,'phase',243.8674,'period',24.07)   ...
);

% 2.2 Building arrays for model computation
names = fieldnames(tidal_constituents);
numC = numel(names);
A   = zeros(1,numC);
omega = zeros(1,numC);
phi = zeros(1,numC);
for k = 1:numC
    name = names{k};
    A(k) = tidal_constituents.(name).amplitude;
    T_h = tidal_constituents.(name).period;
    omega(k) = 2*pi/(T_h*3600);             % rad/s
    phi(k)   = tidal_constituents.(name).phase/180*pi;  % rad
end

% 2.3 Mean sea level
H0 = 0.55;   % m above geoid (from CMEMS) 

%% 3. Setting up Time Vector
epoch = datetime(2025,05,2,0,0,0,'TimeZone','UTC');  % user-defined epoch
model_duration = days(7);              % 4days simulation
t_end = epoch + model_duration;
dt = 1800;      % 1800s for 0.5h time step
temporal_grid = 0:dt:seconds(t_end-epoch);
t = temporal_grid;
nt = length(t);

%% 4. Computing Tidal Elevations
elevations = zeros(size(xGrid,1), size(xGrid,2), nt);
for k = 1:nt
    ht = H0;
    for j = 1:numC
        ht = ht + A(j) * sin(omega(j)*t(k) + phi(j));
    end
    elevations(:,:,k) = ht;
end

%% 5. Setting up Visualization & Video
fig = figure('Color','w','Position',[100 100 1280 720]);
set(fig, 'Resize', 'off');  % Prevents resizing between frames; it's a bug fix, don't touch it.

% Video writer
profiles = VideoWriter.getProfiles;
profileNames = {profiles.Name};
if any(strcmp(profileNames,'MPEG-4'))
    v = VideoWriter('tidal_animation.mp4','MPEG-4');
elseif any(strcmp(profileNames,'Motion JPEG AVI'))
    v = VideoWriter('tidal_animation.avi','Motion JPEG AVI');
else
    error('No supported VideoWriter profiles found.');
end
v.FrameRate = 10;
open(v);

% 5.1 Defining custom colormap
% Colormap: dark blue (deep water), light blue (shallow water), light brown (low land), dark brown (high land)

nLevels = 128;
blueColors  = [linspace(0,0.7,nLevels)', linspace(0.3,0.85,nLevels)', ones(nLevels,1)];
brownColors = [linspace(1,0.55,nLevels)', linspace(0.9,0.27,nLevels)', linspace(0.7,0.07,nLevels)'];
cmap = [blueColors; brownColors];
colormap(cmap);
minElev = min(zGrid(:)); maxElev = max(zGrid(:));

%% 6. Generating Frames & Writing Video (HD with margin for labels)
ax = gca;
set(ax, 'Units', 'normalized', 'Position', [0.08 0.1 0.88 0.85]);  % Leaving margin to prevent label cut off

lonVec = min(lon):dx:max(lon);
latVec = min(lat):dx:max(lat);
for k = 1:nt
    displayData = zGrid - elevations(:,:,k);     % main visualization equation; don't confuse, it's correct here.
    imagesc(lonVec, latVec, displayData, 'Parent', ax);
    axis(ax, 'equal', 'tight'); set(ax,'YDir','normal');
    title(ax, sprintf('Sea Surface Hight at: t=%.1f h', t(k)/3600), 'FontSize', 12);
    xlabel(ax, 'Longitude', 'FontSize', 11); ylabel(ax, 'Latitude', 'FontSize', 11);
    colormap(ax, cmap); colorbar(ax); caxis(ax, [minElev, maxElev]);
    hold(ax, 'on'); plot(ax, coastX, coastY, 'k','LineWidth',0.7); hold(ax, 'off');

    frame = getframe(fig);  % Full figure frame
    writeVideo(v, frame);
end

%% 7. Saving Model Results & Cleaning up
close(v); clear v;
close(fig);
save('tidal_data.mat','xGrid','yGrid','elevations','t');

%% Notes:
