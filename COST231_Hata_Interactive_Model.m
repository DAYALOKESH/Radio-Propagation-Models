function COST231_Hata_Interactive_Model()
% COST231_HATA_INTERACTIVE_MODEL Comprehensive MATLAB implementation of the COST-231 Hata model.
%
% ## Description
% This function launches a self-contained interactive application for calculating
% and visualizing path loss based on the COST-231 Hata propagation model, which
% extends the Okumura-Hata model to higher frequencies. It is typically used for
% frequencies between 1500 MHz and 2000 MHz.
%
% ## Features
% - Interactive UI with an input field for multiple frequencies and a dropdown for terrain.
% - Compare an arbitrary number of frequencies simultaneously.
% - Support for Urban (Metropolitan) and Suburban/Rural environments.
% - Option to view individual terrain plots or a combined comparison.
% - Path loss vs. Distance plots with clear legends, titles, and grids.
% - Functionality to export plots to PNG or PDF formats.
% - Functionality to export calculated path loss data to a MAT file.
% - Adherence to model constraints for frequency and antenna heights.
%
% ## Usage
% Simply run the script in MATLAB. The interactive window will appear.
%
%   >> COST231_Hata_Interactive_Model
%
% ## UI Controls
% 1.  **Frequencies (Edit Field):** Enter one or more frequencies in MHz, separated by commas (e.g., "1800, 1900").
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
%    - Enter "1800" in the 'Frequencies' box and click 'Update Plot'.
%    - The plot will show the path loss for an urban environment at 1800 MHz.
%
% 2. **Multi-Frequency, Multi-Terrain Comparison:**
%    - Keep 'Terrain Selection' on 'All'.
%    - Enter "1700, 1850, 2000" in the 'Frequencies' box and click 'Update Plot'.
%    - The plot will display and compare the path loss across both terrain
%      types for each of the three specified frequencies.
%
% Author: Gemini
% Version: 1.0
% Date: 2025-08-13

% --- Model and Simulation Parameters ---
hb = 52;   % Base station antenna height in meters (hb).
hm = 2.4;  % Mobile station antenna height in meters (hm).
d = 1:0.1:20; % Distance vector from 1 km to 20 km.

% --- Parameter Validation ---
if hb < 30 || hb > 200
    error('Base station height (hb) must be between 30m and 200m.');
end
if hm < 1 || hm > 10
    error('Mobile station height (hm) must be between 1m and 10m.');
end

% --- UI Figure and Layout Setup ---
fig = uifigure('Name', 'COST-231 Hata Propagation Model', 'Position', [150 150 900 600]);
gl = uigridlayout(fig, [4, 4], 'RowHeight', {30, 30, '1x', 40}, 'ColumnWidth', {'1x', '1x', '1x', '1x'});

% --- Create UI Components ---

% Frequency Input Field
label_f = uilabel(gl, 'Text', 'Frequencies (1500-2000 MHz):');
label_f.Layout.Row = 1;
label_f.Layout.Column = 1;
edit_f = uieditfield(gl, 'text', 'Value', '1800, 2000');
edit_f.Layout.Row = 1;
edit_f.Layout.Column = [2, 4];

% Terrain Selection Dropdown
label_terrain = uilabel(gl, 'Text', 'Select Terrain:');
label_terrain.Layout.Row = 2;
label_terrain.Layout.Column = 1;
dropdown_terrain = uidropdown(gl, 'Items', {'All', 'Urban', 'Suburban/Rural'}, 'Value', 'All');
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
        freqs = [];
        freq_str = app.edit_f.Value;
        try
            freqs = str2double(strsplit(freq_str, ','));
            if any(isnan(freqs)) || isempty(freqs)
                uialert(app.fig, 'Invalid frequency input. Please enter numbers separated by commas.', 'Input Error');
                freqs = [];
                return;
            end
            % COST-231 specific frequency validation
            if any(freqs < 1500) || any(freqs > 2000)
                uialert(app.fig, 'Frequencies must be within the COST-231 range (1500-2000 MHz).', 'Input Error');
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
        
        freqs = parseFrequencies();
        if isempty(freqs), return; end % Stop if input is invalid

        terrain_choice = app.dropdown_terrain.Value;

        cla(app.ax);
        hold(app.ax, 'on');
        
        plot_handles = [];
        legend_entries = {};
        
        % Determine which terrains to plot
        terrains_to_plot = {};
        if strcmp(terrain_choice, 'All')
            terrains_to_plot = {'Urban', 'Suburban/Rural'};
        else
            terrains_to_plot = {terrain_choice};
        end

        colors = lines(length(terrains_to_plot));
        line_styles = {'-', '--', ':', '-.'};
        
        for i = 1:length(terrains_to_plot)
            terrain = terrains_to_plot{i};
            color = colors(i,:);
            
            for j = 1:length(freqs)
                f = freqs(j);
                style = line_styles{mod(j-1, length(line_styles)) + 1};
                
                pl = calculatePathLoss(terrain, f, app.d, app.hb, app.hm);
                h = plot(app.ax, app.d, pl, 'Color', color, 'LineWidth', 1.5, 'LineStyle', style);
                plot_handles(end+1) = h;
                legend_entries{end+1} = sprintf('%s - %.0f MHz', terrain, f);
            end
        end
        
        hold(app.ax, 'off');
        
        grid(app.ax, 'on');
        xlabel(app.ax, 'Distance (km)');
        ylabel(app.ax, 'Path Loss (dB)');
        title(app.ax, sprintf('COST-231 Hata Path Loss (hb=%.1fm, hm=%.1fm)', app.hb, app.hm));
        legend(app.ax, plot_handles, legend_entries, 'Location', 'southeast', 'FontSize', 8);
        axis(app.ax, 'tight');
    end

    function savePlot(~, ~)
        [file, path] = uiputfile({'*.png', 'PNG Image (*.png)'; '*.pdf', 'PDF Document (*.pdf)'}, 'Save Plot As');
        if isequal(file, 0), return; end
        
        full_path = fullfile(path, file);
        try
            exportgraphics(app.ax, full_path);
            uialert(app.fig, ['Plot saved successfully to ' full_path], 'Save Successful');
        catch ME
            uialert(app.fig, ['Error saving plot: ' ME.message], 'Save Error');
        end
    end

    function exportData(~, ~)
        freqs = parseFrequencies();
        if isempty(freqs), return; end

        [file, path] = uiputfile('*.mat', 'Export Data As', 'cost231_hata_data.mat');
        if isequal(file, 0), return; end
        
        distance_km = app.d;
        
        data.description = 'Path loss data from COST-231 Hata model simulation.';
        data.parameters.base_station_height_m = app.hb;
        data.parameters.mobile_station_height_m = app.hm;
        data.parameters.frequencies_MHz = freqs;
        data.distance_km = distance_km;
        
        for i = 1:length(freqs)
            f = freqs(i);
            freq_fieldname = sprintf('f%dMHz', round(f));
            
            data.path_loss_dB.urban.(freq_fieldname) = calculatePathLoss('Urban', f, distance_km, app.hb, app.hm);
            data.path_loss_dB.suburban_rural.(freq_fieldname) = calculatePathLoss('Suburban/Rural', f, distance_km, app.hb, app.hm);
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
    % Calculates path loss using the COST-231 Hata model.

    % Correction factor for mobile antenna height a(hm)
    % This is for medium-sized city and suburban areas.
    ahm = (1.1 * log10(f) - 0.7) * hm - (1.56 * log10(f) - 0.8);
    
    % The Cm parameter depends on the terrain type.
    switch lower(terrain)
        case 'urban'
            Cm = 3; % For metropolitan/urban centers
        case 'suburban/rural'
            Cm = 0; % For suburban or rural areas
        otherwise
            error('Invalid terrain type specified.');
    end
    
    % COST-231 Hata formula
    PL = 46.3 + 33.9 * log10(f) - 13.82 * log10(hb) - ahm ...
         + (44.9 - 6.55 * log10(hb)) * log10(d) + Cm;
end
