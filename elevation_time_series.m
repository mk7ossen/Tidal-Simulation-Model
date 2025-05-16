load('tidal_data.mat'); % loads xGrid, yGrid, elevations, t

% Find index of a point of interest (e.g., near the center)
[~, ix] = min(abs(xGrid(1,:) - mean(xGrid(1,:))));
[~, iy] = min(abs(yGrid(:,1) - mean(yGrid(:,1))));

% Extract time series of elevation at this point
elev_point = squeeze(elevations(iy, ix, :));

figure('Color','w');
plot(t/3600, elev_point, 'LineWidth',1.5, 'Color',[0 0.4470 0.7410]); % MATLAB default blue
xlabel('Time (hours)');
ylabel('Sea Surface Height (m)');
title('Tidal Elevation Time Series Data');
grid on;
set(gca,'FontSize',12);
