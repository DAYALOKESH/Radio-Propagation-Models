% -------------------------------------------------------------------------
% MATLAB Implementation of Single Knife-Edge Diffraction Analysis
% -------------------------------------------------------------------------
%
% Description:
% This script models the effect of single knife-edge diffraction on radio
% wave propagation. It reads a terrain profile, identifies the most
% significant obstacle (the "knife edge"), and calculates the resulting
% path loss and electric field strength.
%
% The script generates three plots for comprehensive analysis:
%   1. Terrain Profile: Shows the terrain, antenna locations, and the
%      line-of-sight (LOS) path, highlighting the knife-edge obstacle.
%   2. Field Strength Comparison: Compares the electric field strength with
%      and without diffraction effects.
%   3. Path Loss Comparison: Compares the path loss in free space versus
%      the total path loss including diffraction.
%
% Required File:
%   X.txt - A text file containing space-separated distance (meters) and
%           height (meters) data for the terrain profile. Each line
%           represents a point on the terrain.
%
% Example X.txt content:
% 0 20
% 1000 45
% 2000 60
% 3000 110
% 4000 80
% 5000 55
% 6000 30
% 7000 25
% 8000 40
% 9000 35
% 10000 30
%
% -------------------------------------------------------------------------

clear;
clc;
close all;

%% 1. Define Physical Constants and Parameters
% -------------------------------------------------------------------------
f = 970;          % Frequency in MHz
c = 3e8;          % Speed of light in m/s
lambda = c / (f * 1e6); % Wavelength in meters

% Define transmitter and receiver heights above ground
h_tx_gnd = 52;    % Transmitter height above ground (meters)
h_rx_gnd = 10;    % Receiver height above ground (meters)

%% 2. Read Terrain Elevation Data
% -------------------------------------------------------------------------
% Load terrain data from X.txt. The file should have two columns:
% Column 1: Distance from transmitter (meters)
% Column 2: Ground elevation (meters)
try
    terrain_data = load('X.04');
catch
    error('Error: X.txt not found. Please create the terrain data file in the same directory.');
end

d = terrain_data(:, 1); % Distance vector (m)
h_gnd = terrain_data(:, 2); % Ground height vector (m)

num_points = length(d);
d_total = d(end); % Total path distance

%% 3. Calculate Effective Antenna Heights
% -------------------------------------------------------------------------
% Effective height = Antenna height above ground + Ground elevation
h_tx_eff = h_tx_gnd + h_gnd(1);   % Effective transmitter height
h_rx_eff = h_rx_gnd + h_gnd(end); % Effective receiver height

%% 4. Identify the Knife Edge and Calculate Obstruction Height
% -------------------------------------------------------------------------
% Calculate the height of the direct line-of-sight (LOS) path at each point
los_height = h_tx_eff + (h_rx_eff - h_tx_eff) * (d / d_total);

% Calculate the obstruction height (h) of each terrain point relative to the LOS path
h_obs = h_gnd - los_height;

% Find the point of maximum obstruction, which defines the knife edge
[h_ke, ke_idx] = max(h_obs);

% If h_ke is not positive, there is no obstruction (full LOS)
if h_ke <= 0
    disp('Clear line-of-sight path exists. No diffraction loss calculated.');
    v = -inf; % Set Fresnel parameter to a large negative value
    G_db = 0;   % No diffraction gain/loss
    ke_idx = -1; % Flag to indicate no knife edge
else
    disp(['Knife-edge obstacle identified at distance ', num2str(d(ke_idx)/1000), ' km.']);
    disp(['Obstacle height above LOS: ', num2str(h_ke), ' m.']);
    
    % Get distances to and from the knife edge
    d1 = d(ke_idx);
    d2 = d_total - d1;

    %% 5. Compute Fresnel Parameter (v) for the Knife Edge
    % ---------------------------------------------------------------------
    % Standard formula for the Fresnel parameter
    v = h_ke * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));
    disp(['Calculated Fresnel parameter (v): ', num2str(v)]);

    %% 8. Calculate Diffraction Gain (G_db) using Lee's Approximation
    % ---------------------------------------------------------------------
    % Lee's approximation provides a good estimate for diffraction gain/loss in dB.
    % Note: The result is a gain, which is typically negative (i.e., a loss).
    if v <= -1
        G_db = 0;
    elseif v > -1 && v <= 0
        G_db = 20 * log10(0.5 - 0.62 * v);
    elseif v > 0 && v <= 1
        G_db = 20 * log10(0.5 * exp(-0.95 * v));
    elseif v > 1 && v <= 2.4
        % Added a small epsilon to prevent log of a non-positive number in edge cases
        G_db = 20 * log10(0.4 - sqrt(0.1184 - (0.38 - 0.1 * v)^2 + 1e-9));
    else % v > 2.4
        G_db = 20 * log10(1 / (2.2 * v));
    end
    disp(['Calculated Diffraction Gain (G_db): ', num2str(G_db), ' dB']);
end

%% 6. Calculate Free Space Path Loss (FSPL)
% -------------------------------------------------------------------------
% FSPL (dB) = 20*log10(d_km) + 20*log10(f_MHz) + 32.45
% We calculate it for each point along the path.
% Avoid log10(0) at the first point by starting from a small distance.
d_km = d / 1000;
fl_db = zeros(size(d_km));
non_zero_dist_idx = d_km > 0;
fl_db(non_zero_dist_idx) = 20 * log10(d_km(non_zero_dist_idx)) + 20 * log10(f) + 32.45;

%% 9. Compute Total Path Loss (tl_db)
% -------------------------------------------------------------------------
% Before the knife edge, total loss is just the free space loss.
% After the knife edge, the diffraction loss is added.
tl_db = fl_db; % Initialize with free space loss
if ke_idx > 0 % Apply diffraction loss only if an obstacle exists
    % The diffraction loss is constant for all points in the shadow region
    diffraction_loss = -G_db; % Convert gain (negative) to loss (positive)
    tl_db(ke_idx:end) = tl_db(ke_idx:end) + diffraction_loss;
end

%% 10. Calculate Electric Field Strength
% -------------------------------------------------------------------------
% Field strength is calculated relative to a reference free-space field (E0).
% We normalize E0 at 1 km to be 1 V/m for simplicity.
E0_ref = 1; % Reference field strength at 1 km

% Free space electric field (decays with distance)
E_free_space = zeros(size(d_km));
E_free_space(non_zero_dist_idx) = E0_ref ./ d_km(non_zero_dist_idx);

% Electric field with diffraction
% The field is scaled by the diffraction gain factor |F(v)| = 10^(G_db/20)
E_diffracted = E_free_space;
if ke_idx > 0
    diffraction_factor = 10^(G_db / 20);
    E_diffracted(ke_idx:end) = E_diffracted(ke_idx:end) * diffraction_factor;
end

%% 11. Generate Comprehensive Comparison Visualizations
% -------------------------------------------------------------------------
figure('Name', 'Single Knife-Edge Diffraction Analysis', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

% Plot 1: Terrain Profile and Path Geometry
subplot(3, 1, 1);
plot(d/1000, h_gnd, 'k', 'LineWidth', 2, 'DisplayName', 'Terrain Profile');
hold on;
plot(d/1000, los_height, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Line-of-Sight (LOS)');

% Plot Tx and Rx locations
plot(d(1)/1000, h_tx_eff, 'pg', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Transmitter');
plot(d(end)/1000, h_rx_eff, 'pr', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Receiver');

% Highlight the knife edge
if ke_idx > 0
    plot(d(ke_idx)/1000, h_gnd(ke_idx), 'm^', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Knife Edge');
    % Draw line from LOS to knife edge peak
    plot([d(ke_idx)/1000, d(ke_idx)/1000], [los_height(ke_idx), h_gnd(ke_idx)], 'c-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(d(ke_idx)/1000, h_gnd(ke_idx) + 5, ['h = ' num2str(h_ke, '%.1f') ' m'], 'Color', 'c', 'FontWeight', 'bold');
end

title('Terrain Profile and Path Geometry');
xlabel('Distance (km)');
ylabel('Height (m)');
legend('show', 'Location', 'northwest');
grid on;
axis tight;
hold off;

% Plot 2: Electric Field Strength Comparison
subplot(3, 1, 2);
plot(d(non_zero_dist_idx)/1000, E_free_space(non_zero_dist_idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Free Space Field');
hold on;
plot(d(non_zero_dist_idx)/1000, E_diffracted(non_zero_dist_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'Field with Diffraction');
if ke_idx > 0
    xline(d(ke_idx)/1000, 'm--', 'Label', 'Knife Edge Location');
end
title('Comparison of Electric Field Strength');
xlabel('Distance (km)');
ylabel('Relative Field Strength (V/m)');
legend('show', 'Location', 'northeast');
grid on;
axis tight;
ylim([0, max(E_diffracted)*1.1]); % Set y-axis limit for better visualization
hold off;

% Plot 3: Path Loss Comparison
subplot(3, 1, 3);
plot(d(non_zero_dist_idx)/1000, fl_db(non_zero_dist_idx), 'b-', 'LineWidth', 2, 'DisplayName', 'Free Space Path Loss (FSPL)');
hold on;
plot(d(non_zero_dist_idx)/1000, tl_db(non_zero_dist_idx), 'r--', 'LineWidth', 2, 'DisplayName', 'Total Loss (with Diffraction)');
if ke_idx > 0
    xline(d(ke_idx)/1000, 'm--', 'Label', 'Knife Edge Location');
end
title('Comparison of Path Loss');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
legend('show', 'Location', 'southeast');
grid on;
axis tight;
hold off;

% Set overall figure properties
sgtitle(['Knife-Edge Diffraction Analysis (f = ' num2str(f) ' MHz)'], 'FontSize', 16, 'FontWeight', 'bold');
