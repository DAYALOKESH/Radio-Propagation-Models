function deygout_analysis_matlab()
    % DEYGOUT_ANALYSIS_MATLAB - A comprehensive, single-file implementation of the
    % Deygout method for multiple knife-edge diffraction analysis.
    %
    % This script prompts the user for a terrain file, performs the analysis,
    % and generates detailed visualizations of the terrain profile, obstacles,
    % path loss, diffraction loss, and electric field strength.
    %
    % To Run:
    % 1. Save this entire code as a single .m file (e.g., deygout_analysis_matlab.m).
    % 2. Place your terrain data file (e.g., 'X.04') in the same directory.
    % 3. Run the script from the MATLAB command window by typing: deygout_analysis_matlab

    % --- System Parameters ---
    params.tx_height = 52.0;       % Transmitter height (m)
    params.rx_height = 2.4;        % Receiver height (m)
    params.freq = 970e6;           % Frequency in Hz
    params.k_factor = 4 / 3;       % Effective Earth radius factor
    params.earth_radius = 6371000; % Earth radius in m
    params.tx_power_dbm = 40;      % Assumed transmitter power in dBm
    params.tx_gain_dbi = 0;        % Assumed antenna gain in dBi
    
    % --- Derived Parameters ---
    C = 299792458; % Speed of light in m/s
    params.lambda = C / params.freq;
    params.effective_earth_radius = params.k_factor * params.earth_radius;
    params.eirp_dbm = params.tx_power_dbm + params.tx_gain_dbi;

    % --- Main Execution ---
    try
        % Get terrain file from user
        terrain_file = input('Enter the name of the terrain file (e.g., X.04): ', 's');
        
        % 1. Read and Validate Terrain Data
        fprintf('Reading terrain data from %s...\n', terrain_file);
        terrain_raw = read_terrain_data(terrain_file);

        % 2. Add Antenna Heights and Apply Earth Curvature
        fprintf('Applying antenna heights and Earth curvature corrections...\n');
        terrain_with_antennas = terrain_raw;
        terrain_with_antennas.height(1) = terrain_with_antennas.height(1) + params.tx_height;
        terrain_with_antennas.height(end) = terrain_with_antennas.height(end) + params.rx_height;
        terrain_corrected = apply_earth_curvature(terrain_with_antennas, params);

        % 3. Generate Full Propagation Profile
        fprintf('Analyzing path and calculating propagation profiles...\n');
        profile = generate_full_profile(terrain_corrected, params);
        
        % 4. Display Key Results
        fprintf('\n--- Analysis Complete ---\n');
        fprintf('Total Path Loss at Receiver: %.2f dB\n', profile.total_path_loss(end));
        fprintf('Maximum Diffraction Loss: %.2f dB\n', max(profile.diffraction_loss));
        fprintf('Minimum Field Strength: %.2f dBuV/m\n', min(profile.field_strength));
        fprintf('%d principal obstacles identified.\n', length(profile.obstacles));

        % 5. Visualize Results
        fprintf('Generating plots...\n');
        visualize_results(terrain_corrected, profile, params);
        
    catch ME
        fprintf(2, 'An error occurred: %s\n', ME.message);
        fprintf(2, 'Error in function %s, line %d.\n', ME.stack(1).name, ME.stack(1).line);
    end
end

% -------------------------------------------------------------------------
% HELPER FUNCTIONS
% -------------------------------------------------------------------------

function terrain = read_terrain_data(filename)
    % Reads terrain data from a space-separated file (distance, height).
    if ~exist(filename, 'file')
        error('Terrain file "%s" not found.', filename);
    end
    data = load(filename);
    if size(data, 2) ~= 2
        error('Terrain file must have exactly 2 columns.');
    end
    terrain.distance = data(:, 1)'; % Use row vectors
    terrain.height = data(:, 2)';
end

function corrected_terrain = apply_earth_curvature(terrain, params)
    % Corrects terrain height for Earth's curvature.
    corrected_terrain = terrain;
    total_distance = terrain.distance(end);
    d1 = terrain.distance;
    d2 = total_distance - d1;
    
    curvature_correction = (d1 .* d2) / (2 * params.effective_earth_radius);
    % Do not correct start and end points relative to the line between them
    corrected_terrain.height(2:end-1) = corrected_terrain.height(2:end-1) - curvature_correction(2:end-1);
end

function h = get_clearance_height(obs_d, obs_h, p1, p2)
    % Calculates the clearance height of an obstacle above the LOS path.
    d1 = obs_d - p1.dist;
    d2 = p2.dist - obs_d;
    total_dist = p2.dist - p1.dist;
    if total_dist <= 0
        h = -inf;
        return;
    end
    line_height = p1.height + (p2.height - p1.height) * d1 / total_dist;
    h = obs_h - line_height;
end

function v = calculate_v(h, d1, d2, lambda)
    % Calculates the Fresnel-Kirchoff diffraction parameter v.
    if d1 <= 0 || d2 <= 0 || lambda <= 0
        v = 0;
        return;
    end
    v = h * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));
end

function J = calculate_J_v(v)
    % Calculates single knife-edge diffraction loss J(v) using ITU-R P.526.
    if v <= -0.78
        J = 0;
    else
        J = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
    end
end

function result = deygout_recursive(start_point, end_point, obstacles, params)
    % Recursively calculates Deygout diffraction loss and identifies principal obstacles.
    % Returns a struct with total loss and the list of identified obstacles.
    
    result.loss = 0;
    result.obstacles = [];

    if isempty(obstacles)
        return;
    end

    max_v = -inf;
    principle_edge_idx = -1;

    for i = 1:length(obstacles)
        obs = obstacles(i);
        h = get_clearance_height(obs.dist, obs.height, start_point, end_point);
        if h > 0
            d1 = obs.dist - start_point.dist;
            d2 = end_point.dist - obs.dist;
            v = calculate_v(h, d1, d2, params.lambda);
            if v > max_v
                max_v = v;
                principle_edge_idx = i;
            end
        end
    end

    if principle_edge_idx == -1 || max_v <= -0.78
        return;
    end

    loss_p = calculate_J_v(max_v);
    principle_edge = obstacles(principle_edge_idx);

    obstacles_left = obstacles(1:principle_edge_idx-1);
    obstacles_right = obstacles(principle_edge_idx+1:end);

    left_result = deygout_recursive(start_point, principle_edge, obstacles_left, params);
    right_result = deygout_recursive(principle_edge, end_point, obstacles_right, params);

    result.loss = loss_p + left_result.loss + right_result.loss;
    result.obstacles = [principle_edge, left_result.obstacles, right_result.obstacles];
end

function profile = generate_full_profile(terrain, params)
    % Generates the complete propagation profile.
    tx = struct('dist', terrain.distance(1), 'height', terrain.height(1));
    rx = struct('dist', terrain.distance(end), 'height', terrain.height(end));
    n_points = length(terrain.distance);

    % Step 1: Find all principal obstacles for the entire path
    all_potential_obstacles = [];
    for i = 2:n_points-1
        h = get_clearance_height(terrain.distance(i), terrain.height(i), tx, rx);
        if h > 0
            all_potential_obstacles = [all_potential_obstacles, struct('dist', terrain.distance(i), 'height', terrain.height(i))];
        end
    end
    full_path_analysis = deygout_recursive(tx, rx, all_potential_obstacles, params);
    
    % Sort the identified main obstacles by distance
    if ~isempty(full_path_analysis.obstacles)
        [~, sort_idx] = sort([full_path_analysis.obstacles.dist]);
        main_obstacles = full_path_analysis.obstacles(sort_idx);
    else
        main_obstacles = [];
    end
    profile.obstacles = main_obstacles;

    % Step 2: Generate profile by calculating loss up to each point
    profile.diffraction_loss = zeros(1, n_points);
    profile.total_path_loss = zeros(1, n_points);
    profile.field_strength = zeros(1, n_points);

    for i = 2:n_points
        current_point = struct('dist', terrain.distance(i), 'height', terrain.height(i));
        
        % Consider only obstacles that are before the current point
        obstacles_for_subpath_idx = [main_obstacles.dist] < current_point.dist;
        obstacles_for_subpath = main_obstacles(obstacles_for_subpath_idx);
        
        sub_path_analysis = deygout_recursive(tx, current_point, obstacles_for_subpath, params);
        diff_loss = sub_path_analysis.loss;
        profile.diffraction_loss(i) = diff_loss;
        
        dist_m = current_point.dist - tx.dist;
        fs_loss = 20*log10(dist_m) + 20*log10(params.freq) - 147.55;
        
        total_loss = fs_loss + diff_loss;
        profile.total_path_loss(i) = total_loss;

        e_strength = params.eirp_dbm - total_loss + 20*log10(params.freq / 1e6) + 77.2;
        profile.field_strength(i) = e_strength;
    end
    profile.field_strength(1) = profile.field_strength(2); % Pad start for plotting
end

function visualize_results(terrain, profile, params)
    % Creates two plots to visualize the analysis results.
    distances_km = terrain.distance / 1000;
    tx = struct('dist', terrain.distance(1), 'height', terrain.height(1));
    rx = struct('dist', terrain.distance(end), 'height', terrain.height(end));
    
    % Calculate Line of Sight path
    los_path = tx.height + (rx.height - tx.height) * (terrain.distance - tx.dist) / (rx.dist - tx.dist);

    % --- Plot 1: Terrain Profile, LOS, and Obstacles ---
    figure('Name', 'Terrain Profile and Obstacles', 'Position', [100, 500, 900, 400]);
    hold on;
    
    % Plot terrain as a filled area
    area(distances_km, terrain.height, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', 'none', 'DisplayName', 'Terrain Profile');
    
    % Plot LOS
    plot(distances_km, los_path, 'r--', 'LineWidth', 2, 'DisplayName', 'Line of Sight');
    
    % Plot Obstacles
    if ~isempty(profile.obstacles)
        obs_dists_km = [profile.obstacles.dist] / 1000;
        obs_heights = [profile.obstacles.height];
        scatter(obs_dists_km, obs_heights, 100, 'p', 'filled', ...
            'MarkerFaceColor', [1.0, 0.5, 0], ...
            'MarkerEdgeColor', 'k', 'DisplayName', 'Principal Obstacles');
    end
    
    grid on;
    box on;
    hold off;
    xlabel('Distance (km)');
    ylabel('Height (m)');
    title(sprintf('Terrain Profile Analysis (f = %.0f MHz)', params.freq/1e6));
    legend('show', 'Location', 'best');
    xlim([distances_km(1), distances_km(end)]);

    % --- Plot 2: Loss and Field Strength Profiles ---
    figure('Name', 'Propagation Profiles', 'Position', [100, 100, 900, 400]);
    
    % Left Y-axis for Loss
    yyaxis left;
    hold on;
    plot(distances_km, profile.diffraction_loss, 'LineWidth', 2, 'DisplayName', 'Diffraction Loss');
    plot(distances_km, profile.total_path_loss, 'LineWidth', 2, 'DisplayName', 'Total Path Loss');
    hold off;
    ylabel('Loss (dB)');
    ax = gca;
    ax.YColor = 'k'; % Black color for left axis

    % Right Y-axis for Field Strength
    yyaxis right;
    hold on;
    plot(distances_km, profile.field_strength, 'LineWidth', 2, 'DisplayName', 'Field Strength');
    hold off;
    ylabel('Field Strength (dBuV/m)');
    ax = gca;
    ax.YColor = 'k'; % Black color for right axis
    
    % Add vertical lines for obstacles
    if ~isempty(profile.obstacles)
        obs_dists_km = [profile.obstacles.dist] / 1000;
        for i = 1:length(obs_dists_km)
            xline(obs_dists_km(i), '--', {'Obstacle', num2str(i)}, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        end
    end

    grid on;
    box on;
    xlabel('Distance (km)');
    title('Propagation Loss and Field Strength Profiles');
    legend('show', 'Location', 'best');
    xlim([distances_km(1), distances_km(end)]);
end