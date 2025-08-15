function OkumuraHata_Interactive_Model()
% OKUMURAHATA_INTERACTIVE_MODEL Comprehensive MATLAB implementation of the Okumura-Hata model.
%
% ## Description
% This function launches a self-contained interactive application for calculating
% and visualizing path loss based on the Okumura-Hata propagation model. It
% supports urban, suburban, and rural environments and allows for detailed
% frequency-based comparative analysis.
%
% ## Features
% - Interactive UI with an input field for multiple frequencies and a dropdown for terrain.
% - Compare an arbitrary number of frequencies simultaneously.
% - Support for Urban, Suburban, and Rural terrain models.
% - Option to view individual terrain plots or a combined comparison.
% - Path loss vs. Distance plots with clear legends, titles, and grids.
% - Functionality to export plots to PNG or PDF formats.
% - Functionality to export calculated path loss data to a MAT file.
% - Adherence to model constraints for frequency and antenna heights.
%
% ## Usage
% Simply run the script in MATLAB. The interactive window will appear.
%
%   >> OkumuraHata_Interactive_Model
%
% ## UI Controls
% 1.  **Frequencies (Edit Field):** Enter one or more frequencies in MHz, separated by commas (e.g., "900, 1500").
% 2.  **Terrain Selection (Dropdown):** Choose which terrain(s) to display.
% 3.  **Update Plot (Button):** Redraws the plot with the specified frequencies.
% 4.  **Save Plot (Button):** Saves the current view of the plot as a PNG or PDF.
% 5.  **Export Data (Button):** Saves the calculated path loss data for the
%     currently specified parameters into a .mat file.
%
% ## Example Scenarios
%
% 1. **Single Frequency Analysis:**
%    - Set 'Terrain Selection' to 'Urban'.
%    - Enter "900" in the 'Frequencies' box and click 'Update Plot'.
%    - The plot will show the path loss for an urban environment at 900 MHz.
%
% 2. **Multi-Frequency, Multi-Terrain Comparison:**
%    - Keep 'Terrain Selection' on 'All'.
%    - Enter "900, 1500, 1800" in the 'Frequencies' box and click 'Update Plot'.
%    - The plot will display and compare the path loss across all three
%      terrain types for each of the three specified frequencies.
%
% Author: Gemini
% Version: 2.1
% Date: 2025-08-13

% --- Model and Simulation Parameters ---
hb = 52;   % Base station antenna height in meters (hb).
hm = 2.4;  % Mobile station antenna height in meters (hm).
d = 1:0.1:40; % Distance vector from 1 km to 40 km.

% --- Parameter Validation ---
if hb < 30 || hb > 200
    error('Base station height (hb) must be between 30m and 200m.');
end
if hm < 1 || hm > 10
    error('Mobile station height (hm) must be between 1m and 10m.');
end

% --- UI Figure and Layout Setup ---
fig = uifigure('Name', 'Okumura-Hata Propagation Model', 'Position', [100 100 900 600]);
gl = uigridlayout(fig, [4, 4], 'RowHeight', {30, 30, '1x', 40}, 'ColumnWidth', {'1x', '1x', '1x', '1x'});

% --- Create UI Components ---

% Frequency Input Field
label_f = uilabel(gl, 'Text', 'Frequencies (MHz, comma-separated):');
label_f.Layout.Row = 1;
label_f.Layout.Column = 1;
edit_f = uieditfield(gl, 'text', 'Value', '900, 1500');
edit_f.Layout.Row = 1;
edit_f.Layout.Column = [2, 4];

% Terrain Selection Dropdown
label_terrain = uilabel(gl, 'Text', 'Select Terrain:');
label_terrain.Layout.Row = 2;
label_terrain.Layout.Column = 1;
dropdown_terrain = uidropdown(gl, 'Items', {'All', 'Urban', 'Suburban', 'Rural'}, 'Value', 'All');
dropdown_terrain.Layout.Row = 2;
dropdown_terrain.Layout.Column = 2;

% Plotting Axes
ax = uiaxes(gl);
ax.Layout.Row = 3;
ax.Layout.Column = [1, 4];
grid(ax, 'on');
xlabel(ax, 'Distance (km)');
ylabel(ax, 'Path Loss (dB)');

% Action Buttons
update_button = uibutton(gl, 'Text', 'Update Plot', 'ButtonPushedFcn', @updatePlot);
update_button.Layout.Row = 4;
update_button.Layout.Column = 1;

save_plot_button = uibutton(gl, 'Text', 'Save Plot', 'ButtonPushedFcn', @savePlot);
save_plot_button.Layout.Row = 4;
save_plot_button.Layout.Column = 2;

export_data_button = uibutton(gl, 'Text', 'Export Data', 'ButtonPushedFcn', @exportData);
export_data_button.Layout.Row = 4;
export_data_button.Layout.Column = 3;

% Store handles and data in a structure for easy access in callbacks
app.fig = fig;
app.ax = ax;
app.d = d;
app.hb = hb;
app.hm = hm;
app.edit_f = edit_f;
app.dropdown_terrain = dropdown_terrain;

% Set callbacks for interactive components
dropdown_terrain.ValueChangedFcn = @updatePlot;

% Perform initial plot
updatePlot([], []);


% --- Nested Callback Functions ---

    function freqs = parseFrequencies()
        % Helper function to parse and validate frequencies from the input field.
        % It uses the 'app' struct from the parent function's scope.
        freqs = [];
        freq_str = app.edit_f.Value;
        try
            % Convert comma-separated string to a numeric array
            freqs = str2double(strsplit(freq_str, ','));
            % Check for any parsing errors (e.g., non-numeric input)
            if any(isnan(freqs)) || isempty(freqs)
                uialert(app.fig, 'Invalid frequency input. Please enter numbers separated by commas.', 'Input Error');
                freqs = [];
                return;
            end
        catch ME
            uialert(app.fig, ['Error parsing frequencies: ' ME.message], 'Input Error');
            freqs = [];
        end
    end

    function updatePlot(~, ~)
        % This function is called to recalculate and redraw the path loss curves.
        % It automatically has access to the 'app' struct from its parent.
        
        freqs = parseFrequencies();
        if isempty(freqs), return; end % Stop if input is invalid

        terrain_choice = app.dropdown_terrain.Value;

        % Clear the axes for the new plot
        cla(app.ax);
        hold(app.ax, 'on');
        
        plot_handles = [];
        legend_entries = {};
        
        % Determine which terrains to plot
        terrains_to_plot = {};
        if strcmp(terrain_choice, 'All')
            terrains_to_plot = {'Urban', 'Suburban', 'Rural'};
        else
            terrains_to_plot = {terrain_choice};
        end

        % Define line styles and colors for visual distinction
        colors = lines(length(terrains_to_plot));
        line_styles = {'-', '--', ':', '-.'};
        
        % Loop through each terrain
        for i = 1:length(terrains_to_plot)
            terrain = terrains_to_plot{i};
            color = colors(i,:);
            
            % Loop through each frequency
            for j = 1:length(freqs)
                f = freqs(j);
                style = line_styles{mod(j-1, length(line_styles)) + 1};
                
                % Calculate and plot path loss
                pl = calculatePathLoss(terrain, f, app.d, app.hb, app.hm);
                h = plot(app.ax, app.d, pl, 'Color', color, 'LineWidth', 1.5, 'LineStyle', style);
                plot_handles(end+1) = h;
                legend_entries{end+1} = sprintf('%s - %.0f MHz', terrain, f);
            end
        end
        
        hold(app.ax, 'off');
        
        % Finalize plot aesthetics
        grid(app.ax, 'on');
        xlabel(app.ax, 'Distance (km)');
        ylabel(app.ax, 'Path Loss (dB)');
        title(app.ax, sprintf('Okumura-Hata Path Loss Comparison (hb=%.1fm, hm=%.1fm)', app.hb, app.hm));
        legend(app.ax, plot_handles, legend_entries, 'Location', 'southeast', 'FontSize', 8);
        axis(app.ax, 'tight');
    end

    function savePlot(~, ~)
        % Callback to save the current plot to a file.
        [file, path] = uiputfile({'*.png', 'PNG Image (*.png)'; '*.pdf', 'PDF Document (*.pdf)'}, 'Save Plot As');
        if isequal(file, 0) || isequal(path, 0)
            return; % User cancelled
        end
        
        full_path = fullfile(path, file);
        try
            exportgraphics(app.ax, full_path);
            uialert(app.fig, ['Plot saved successfully to ' full_path], 'Save Successful');
        catch ME
            uialert(app.fig, ['Error saving plot: ' ME.message], 'Save Error');
        end
    end

    function exportData(~, ~)
        % Callback to export the calculated data to a .mat file.
        freqs = parseFrequencies();
        if isempty(freqs), return; end % Stop if input is invalid

        [file, path] = uiputfile('*.mat', 'Export Data As', 'hata_model_data.mat');
        if isequal(file, 0) || isequal(path, 0)
            return; % User cancelled
        end
        
        % Gather data from current settings
        distance_km = app.d;
        
        data.description = 'Path loss data from Okumura-Hata model simulation.';
        data.parameters.base_station_height_m = app.hb;
        data.parameters.mobile_station_height_m = app.hm;
        data.parameters.frequencies_MHz = freqs;
        data.distance_km = distance_km;
        
        % Calculate all path losses for export
        for i = 1:length(freqs)
            f = freqs(i);
            freq_fieldname = sprintf('f%dMHz', round(f));
            
            data.path_loss_dB.urban.(freq_fieldname) = calculatePathLoss('Urban', f, distance_km, app.hb, app.hm);
            data.path_loss_dB.suburban.(freq_fieldname) = calculatePathLoss('Suburban', f, distance_km, app.hb, app.hm);
            data.path_loss_dB.rural.(freq_fieldname) = calculatePathLoss('Rural', f, distance_km, app.hb, app.hm);
        end

        full_path = fullfile(path, file);
        try
            save(full_path, 'data');
            uialert(app.fig, ['Data exported successfully to ' full_path], 'Export Successful');
        catch ME
            uialert(app.fig, ['Error exporting data: ' ME.message], 'Export Error');
        end
    end
end


% --- Core Model Calculation Functions ---

function PL = calculatePathLoss(terrain, f, d, hb, hm)
    % Master function to dispatch calculation based on terrain type.
    
    % Validate frequency
    if f < 150 || f > 1500
        warning('Frequency %.1f MHz is outside the valid Okumura-Hata range (150-1500 MHz). Results may be inaccurate.', f);
    end

    % Calculate the urban path loss first, as suburban/rural models are
    % corrections based on it.
    PL_urban = pathLossUrban(f, d, hb, hm);

    switch lower(terrain)
        case 'urban'
            PL = PL_urban;
        case 'suburban'
            % Suburban Correction
            C_suburban = -2 * (log10(f / 28)).^2 - 5.4;
            PL = PL_urban + C_suburban;
        case 'rural'
            % Rural Correction
            C_rural = -4.78 * (log10(f)).^2 + 18.33 * log10(f) - 40.94;
            PL = PL_urban + C_rural;
        otherwise
            error('Invalid terrain type specified.');
    end
end

function PL_urban = pathLossUrban(f, d, hb, hm)
    % Calculates path loss for a standard Urban environment.
    % This is the base formula for the Hata model.

    % Calculate a(hm), the correction factor for mobile antenna height.
    % This formula is for a medium-to-small city.
    ahm = (1.1 * log10(f) - 0.7) * hm - (1.56 * log10(f) - 0.8);
    
    % Hata formula for Urban areas
    PL_urban = 69.55 + 26.16 * log10(f) - 13.82 * log10(hb) - ahm + (44.9 - 6.55 * log10(hb)) * log10(d);
end
