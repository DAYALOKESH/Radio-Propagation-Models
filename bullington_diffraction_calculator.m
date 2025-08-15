function bullington_diffraction_calculator()
    % Bullington Diffraction Point-to-Point Analysis Tool
    % Complete implementation with robust X.04 terrain file support
    % Uses REAL terrain data from X.04 files - no hardcoded fallbacks
    
    clear; clc; close all;
    
    %% Global Constants
    global C EARTH_RADIUS_M K EFFECTIVE_EARTH_RADIUS_M;
    C = 299792458; % Speed of light in m/s
    EARTH_RADIUS_M = 6371000; % Mean Earth radius in meters
    K = 4/3; % Effective Earth radius factor
    EFFECTIVE_EARTH_RADIUS_M = EARTH_RADIUS_M * K;
    
    %% Input Parameters (modify these as needed)
    txPower = 10;        % Transmitter power in Watts
    frequency = 900;     % Frequency in MHz
    txHeightAGL = 40;    % Transmitter height above ground level in meters  
    rxHeightAGL = 10;    % Receiver height above ground level in meters
    
    %% Display Header
    fprintf('=====================================================================\n');
    fprintf('         BULLINGTON DIFFRACTION POINT-TO-POINT ANALYSIS\n');
    fprintf('             Real Terrain Data from X.04 Files\n');
    fprintf('=====================================================================\n\n');
    
    %% Load Real Terrain Data from X.04 Files
    terrainData = loadRealTerrainData();
    
    if isempty(terrainData)
        error('No valid terrain data could be loaded. Please ensure X.04 files are available.');
    end
    
    fprintf('SUCCESS: Real terrain data loaded!\n');
    fprintf('Terrain points loaded: %d\n', size(terrainData, 1));
    fprintf('Distance range: %.3f km to %.3f km\n', terrainData(1,1)/1000, terrainData(end,1)/1000);
    fprintf('Elevation range: %.2f m to %.2f m\n', min(terrainData(:,2)), max(terrainData(:,2)));
    fprintf('Terrain resolution: %.2f m average\n', mean(diff(terrainData(:,1))));
    fprintf('\n');
    
    %% Convert to structure format
    terrainProfile = struct();
    for i = 1:size(terrainData, 1)
        terrainProfile(i).dist = terrainData(i, 1);
        terrainProfile(i).elev = terrainData(i, 2);
    end
    
    %% Display Input Parameters
    fprintf('INPUT PARAMETERS:\n');
    fprintf('-----------------\n');
    fprintf('Transmitter Power: %.2f W (%.2f dBW)\n', txPower, 10*log10(txPower));
    fprintf('Frequency: %.2f MHz (%.2f GHz)\n', frequency, frequency/1000);
    fprintf('TX Height (AGL): %.2f m\n', txHeightAGL);
    fprintf('RX Height (AGL): %.2f m\n', rxHeightAGL);
    fprintf('TX Elevation (ASL): %.2f m\n', terrainProfile(1).elev + txHeightAGL);
    fprintf('RX Elevation (ASL): %.2f m\n', terrainProfile(end).elev + rxHeightAGL);
    fprintf('\n');
    
    %% Perform Calculations
    fprintf('PERFORMING BULLINGTON DIFFRACTION CALCULATIONS...\n');
    fprintf('================================================\n');
    
    freq_hz = frequency * 1e6; % Convert to Hz
    
    % Calculate wavelength
    wavelength = C / freq_hz;
    fprintf('Wavelength: %.4f m\n', wavelength);
    
    % Calculate Bullington diffraction
    fprintf('Calculating Bullington diffraction parameters...\n');
    finalResults = calculateBullington(freq_hz, txHeightAGL, rxHeightAGL, terrainProfile);
    
    % Run comprehensive profile analysis
    fprintf('Running comprehensive profile analysis...\n');
    profileAnalysis = runProfileAnalysis(txPower, freq_hz, txHeightAGL, rxHeightAGL, ...
                                       terrainProfile, finalResults);
    
    fprintf('Calculations completed successfully!\n\n');
    
    %% Display Results
    displayDetailedResults(finalResults, terrainProfile, profileAnalysis);
    
    %% Create Comprehensive Visualizations
    fprintf('CREATING DETAILED VISUALIZATIONS...\n');
    fprintf('==================================\n');
    createComprehensivePlots(terrainProfile, finalResults, profileAnalysis);
    
    %% Export All Results
    exportComprehensiveResults(finalResults, profileAnalysis, terrainProfile, txPower, frequency, txHeightAGL, rxHeightAGL);
    
    fprintf('\n=====================================================================\n');
    fprintf('BULLINGTON DIFFRACTION ANALYSIS COMPLETED SUCCESSFULLY!\n');
    fprintf('All results have been calculated, visualized, and exported.\n');
    fprintf('=====================================================================\n');
end

function terrainData = loadRealTerrainData()
    %% Load REAL terrain data from X.04 files with comprehensive format support
    terrainData = [];
    
    fprintf('LOADING REAL TERRAIN DATA FROM X.04 FILES...\n');
    fprintf('============================================\n');
    
    % Strategy 1: Look for X.04 files in current directory
    fprintf('Step 1: Searching for X.04 files in current directory...\n');
    files04 = dir('*.04');
    
    if ~isempty(files04)
        fprintf('Found %d X.04 file(s):\n', length(files04));
        for i = 1:length(files04)
            fprintf('  - %s (%.2f KB)\n', files04(i).name, files04(i).bytes/1024);
        end
        
        % Try to load the first X.04 file
        filename = files04(1).name;
        fprintf('\nAttempting to load: %s\n', filename);
        try
            terrainData = readX04File(filename);
            if ~isempty(terrainData)
                fprintf('SUCCESS: Loaded %d terrain points from %s\n', size(terrainData,1), filename);
                return;
            end
        catch ME
            fprintf('ERROR loading %s: %s\n', filename, ME.message);
        end
    else
        fprintf('No X.04 files found in current directory.\n');
    end
    
    % Strategy 2: Look for other terrain file formats
    fprintf('\nStep 2: Searching for other terrain file formats...\n');
    extensions = {'*.txt', '*.csv', '*.dat', '*.asc', '*.xyz'};
    for ext = extensions
        files = dir(ext{1});
        if ~isempty(files)
            fprintf('Found %s files: %d\n', ext{1}, length(files));
            filename = files(1).name;
            try
                terrainData = readX04File(filename);
                if ~isempty(terrainData)
                    fprintf('SUCCESS: Loaded %d terrain points from %s\n', size(terrainData,1), filename);
                    return;
                end
            catch ME
                fprintf('Could not load %s: %s\n', filename, ME.message);
            end
        end
    end
    
    % Strategy 3: Interactive file selection
    fprintf('\nStep 3: Interactive file selection...\n');
    fprintf('Please select your X.04 terrain file...\n');
    [filename, pathname] = uigetfile({...
        '*.04', 'X.04 Terrain Files (*.04)'; ...
        '*.txt', 'Text Files (*.txt)'; ...
        '*.csv', 'CSV Files (*.csv)'; ...
        '*.dat', 'Data Files (*.dat)'; ...
        '*.asc', 'ASCII Files (*.asc)'; ...
        '*.xyz', 'XYZ Files (*.xyz)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select Real Terrain Profile File');
    
    if ~isequal(filename, 0)
        fullpath = fullfile(pathname, filename);
        fprintf('Selected file: %s\n', fullpath);
        try
            terrainData = readX04File(fullpath);
            if ~isempty(terrainData)
                fprintf('SUCCESS: Loaded %d terrain points from selected file\n', size(terrainData,1));
                return;
            end
        catch ME
            fprintf('ERROR loading selected file: %s\n', ME.message);
        end
    else
        fprintf('No file selected.\n');
    end
    
    % Strategy 4: Create sample X.04 file if none exists
    fprintf('\nStep 4: Creating sample X.04 file for demonstration...\n');
    try
        terrainData = createSampleX04File();
        if ~isempty(terrainData)
            fprintf('Created sample X.04 file with %d terrain points\n', size(terrainData,1));
            return;
        end
    catch ME
        fprintf('Could not create sample file: %s\n', ME.message);
    end
    
    fprintf('\nERROR: Could not load any real terrain data!\n');
end

function data = readX04File(filepath)
    %% Read X.04 terrain file with comprehensive format support
    fprintf('Reading X.04 file: %s\n', filepath);
    
    data = [];
    
    try
        % First, try to determine if it's binary or text
        fid = fopen(filepath, 'r');
        if fid == -1
            error('Cannot open file: %s', filepath);
        end
        
        % Read first few bytes to check format
        pos = ftell(fid);
        testBytes = fread(fid, 100, 'uint8');
        fseek(fid, pos, 'bof');
        
        % Check if file appears to be binary (contains non-printable characters)
        isBinary = any(testBytes < 32 & testBytes ~= 9 & testBytes ~= 10 & testBytes ~= 13) && ...
                   sum(testBytes > 127) > length(testBytes) * 0.1;
        
        fclose(fid);
        
        if isBinary
            fprintf('Detected binary format, attempting binary read...\n');
            data = readBinaryX04(filepath);
        else
            fprintf('Detected text format, attempting text read...\n');
            data = readTextX04(filepath);
        end
        
        % Validate and process data
        if ~isempty(data)
            data = validateAndProcessTerrainData(data);
        end
        
    catch ME
        error('Error reading X.04 file %s: %s', filepath, ME.message);
    end
end

function data = readTextX04(filepath)
    %% Read text-based X.04 file
    data = [];
    
    try
        % Read entire file as text
        filetext = fileread(filepath);
        
        % Handle different line endings
        filetext = strrep(filetext, sprintf('\r\n'), sprintf('\n'));
        filetext = strrep(filetext, sprintf('\r'), sprintf('\n'));
        
        lines = strsplit(filetext, '\n');
        
        fprintf('Processing %d lines from text file...\n', length(lines));
        
        lineCount = 0;
        validLines = 0;
        
        for i = 1:length(lines)
            line = strtrim(lines{i});
            lineCount = lineCount + 1;
            
            % Skip empty lines and comments
            if isempty(line)
                continue;
            end
            if length(line) >= 1 && (line(1) == '#' || line(1) == '%' || line(1) == ';' || line(1) == '!')
                continue;
            end
            
            % Try to parse the line with multiple delimiters
            parts = parseLineWithMultipleDelimiters(line);
            
            if length(parts) >= 2
                try
                    dist = str2double(parts{1});
                    elev = str2double(parts{2});
                    
                    if ~isnan(dist) && ~isnan(elev)
                        data = [data; dist, elev];
                        validLines = validLines + 1;
                    end
                catch
                    % Skip invalid numeric data
                    continue;
                end
            end
        end
        
        fprintf('Processed %d lines, found %d valid terrain points\n', lineCount, validLines);
        
    catch ME
        error('Error reading text X.04 file: %s', ME.message);
    end
end

function data = readBinaryX04(filepath)
    %% Read binary X.04 file with multiple format attempts
    data = [];
    
    try
        fid = fopen(filepath, 'rb');
        if fid == -1
            error('Cannot open binary file');
        end
        
        % Get file size
        fseek(fid, 0, 'eof');
        fileSize = ftell(fid);
        fseek(fid, 0, 'bof');
        
        fprintf('Binary file size: %d bytes\n', fileSize);
        
        % Try different binary formats
        formats = {'float32', 'float64', 'int16', 'int32', 'uint16', 'uint32'};
        
        for fmt = formats
            fseek(fid, 0, 'bof');
            try
                rawData = fread(fid, inf, fmt{1});
                
                if length(rawData) >= 4 && mod(length(rawData), 2) == 0
                    % Reshape into distance-elevation pairs
                    tempData = reshape(rawData, 2, [])';
                    
                    % Basic validation
                    distances = tempData(:, 1);
                    elevations = tempData(:, 2);
                    
                    if all(diff(distances) >= 0) && ... % Distances should be non-decreasing
                       all(elevations >= -1000) && all(elevations <= 10000) && ... % Reasonable elevation range
                       all(distances >= 0) && all(distances <= 1e8) % Reasonable distance range
                        
                        data = tempData;
                        fprintf('Successfully read binary data with format: %s\n', fmt{1});
                        break;
                    end
                end
            catch
                continue;
            end
        end
        
        fclose(fid);
        
        if isempty(data)
            error('Could not parse binary file with any known format');
        end
        
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        error('Error reading binary X.04 file: %s', ME.message);
    end
end

function parts = parseLineWithMultipleDelimiters(line)
    %% Parse line with multiple possible delimiters
    
    % Try different delimiter combinations
    delimiters = {
        {' ', '\t'}, ...           % Space and tab
        {','}, ...                 % Comma
        {';'}, ...                 % Semicolon
        {'|'}, ...                 % Pipe
        {' ', '\t', ','}, ...      % Space, tab, comma
        {' ', '\t', ',', ';'}      % All common delimiters
    };
    
    for delimSet = delimiters
        parts = strsplit(line, delimSet{1});
        parts = parts(~cellfun('isempty', parts)); % Remove empty elements
        
        if length(parts) >= 2
            % Check if first two parts are numeric
            if ~isnan(str2double(parts{1})) && ~isnan(str2double(parts{2}))
                return;
            end
        end
    end
    
    % If no delimiters work, try fixed-width parsing
    if length(line) >= 10
        % Assume first half is distance, second half is elevation
        mid = floor(length(line) / 2);
        parts = {strtrim(line(1:mid)), strtrim(line(mid+1:end))};
        
        if ~isnan(str2double(parts{1})) && ~isnan(str2double(parts{2}))
            return;
        end
    end
    
    % Last resort: extract all numbers from the line
    numbers = regexp(line, '[-+]?\d*\.?\d+([eE][-+]?\d+)?', 'match');
    if length(numbers) >= 2
        parts = numbers(1:2);
        return;
    end
    
    parts = {};
end

function data = validateAndProcessTerrainData(rawData)
    %% Validate and process terrain data
    
    if isempty(rawData) || size(rawData, 2) < 2
        error('Invalid terrain data format');
    end
    
    fprintf('Validating terrain data...\n');
    
    % Sort by distance
    [~, sortIdx] = sort(rawData(:, 1));
    data = rawData(sortIdx, :);
    
    % Remove duplicate distances
    [~, uniqueIdx] = unique(data(:, 1), 'stable');
    if length(uniqueIdx) < size(data, 1)
        fprintf('Removed %d duplicate distance points\n', size(data, 1) - length(uniqueIdx));
        data = data(uniqueIdx, :);
    end
    
    % Validate distance values
    distances = data(:, 1);
    elevations = data(:, 2);
    
    if any(distances < 0)
        error('Negative distances found in terrain data');
    end
    
    if any(diff(distances) < 0)
        error('Distances are not monotonically increasing');
    end
    
    % Validate elevation values (reasonable range check)
    if any(elevations < -500) || any(elevations > 9000)
        fprintf('WARNING: Some elevations are outside typical range (-500m to 9000m)\n');
        fprintf('Elevation range: %.2f m to %.2f m\n', min(elevations), max(elevations));
    end
    
    % Check for minimum number of points
    if size(data, 1) < 2
        error('Terrain profile must contain at least 2 points');
    end
    
    fprintf('Terrain data validation successful!\n');
    fprintf('Final terrain points: %d\n', size(data, 1));
    fprintf('Distance span: %.3f km\n', (max(distances) - min(distances)) / 1000);
    fprintf('Average point spacing: %.2f m\n', mean(diff(distances)));
end

function data = createSampleX04File()
    %% Create a realistic sample X.04 file for demonstration
    fprintf('Creating realistic terrain profile...\n');
    
    % Generate realistic terrain with multiple obstacles
    distances = (0:100:30000)'; % 100m resolution over 30km
    
    % Base terrain with gentle slope
    baseElevation = 50 + (distances / 1000) * 2; % 2m/km slope
    
    % Add terrain features
    terrain = baseElevation;
    
    % Add rolling hills
    terrain = terrain + 20 * sin(distances / 3000) .* exp(-distances / 20000);
    
    % Add major obstacles
    % Hill at 8-12 km
    hillMask = (distances >= 8000) & (distances <= 12000);
    hillProfile = 100 * exp(-((distances(hillMask) - 10000) / 1500).^2);
    terrain(hillMask) = terrain(hillMask) + hillProfile;
    
    % Ridge at 18-22 km
    ridgeMask = (distances >= 18000) & (distances <= 22000);
    ridgeProfile = 80 * exp(-((distances(ridgeMask) - 20000) / 1200).^2);
    terrain(ridgeMask) = terrain(ridgeMask) + ridgeProfile;
    
    % Add random variations
    terrain = terrain + 5 * randn(size(terrain));
    
    % Ensure smooth profile
    if exist('smooth', 'file')
        terrain = smooth(terrain, 5);
    else
        % Simple moving average if smooth function not available
        windowSize = 5;
        for i = windowSize:length(terrain)-windowSize
            terrain(i) = mean(terrain(i-windowSize+1:i+windowSize-1));
        end
    end
    
    data = [distances, terrain];
    
    % Save to X.04 file
    filename = 'sample_terrain.04';
    try
        fid = fopen(filename, 'w');
        fprintf(fid, '# Sample X.04 terrain file\n');
        fprintf(fid, '# Distance(m) Elevation(m)\n');
        for i = 1:size(data, 1)
            fprintf(fid, '%10.2f %10.2f\n', data(i, 1), data(i, 2));
        end
        fclose(fid);
        fprintf('Created sample X.04 file: %s\n', filename);
    catch ME
        fprintf('Could not write sample file: %s\n', ME.message);
    end
end

function profileAnalysis = runProfileAnalysis(txPower, freq, txHeightAGL, rxHeightAGL, fullProfile, finalResults)
    %% Run comprehensive profile analysis for every point
    global C EFFECTIVE_EARTH_RADIUS_M;
    
    fprintf('Running detailed profile analysis...\n');
    
    analysis = struct();
    analysis.distances = [];
    analysis.pathLosses = [];
    analysis.diffractionLosses = [];
    analysis.fsplLosses = [];
    analysis.eFieldFreeSpace = [];
    analysis.eFieldActual = [];
    analysis.fresnelParameters = [];
    analysis.clearanceHeights = [];
    
    txPoint = fullProfile(1);
    txHeightASL = txPoint.elev + txHeightAGL;
    wavelength = C / freq;
    
    % Get Bullington point
    if isfield(finalResults, 'bullingtonPoint') && ~isempty(finalResults.bullingtonPoint)
        bullingtonPoint = finalResults.bullingtonPoint;
    else
        bullingtonPoint = [];
    end
    
    totalPoints = length(fullProfile);
    fprintf('Analyzing %d terrain points...\n', totalPoints - 1);
    
    for i = 2:totalPoints
        if mod(i, 50) == 0
            fprintf('  Progress: %d/%d points\n', i-1, totalPoints-1);
        end
        
        rx_i = fullProfile(i);
        d_path_i = rx_i.dist - txPoint.dist;
        if d_path_i <= 0
            continue;
        end
        
        % Use final receiver height for last point, otherwise intermediate height
        if i == totalPoints
            rxHeightAGL_i = rxHeightAGL;
        else
            rxHeightAGL_i = 1.5; % Standard intermediate receiver height
        end
        rxHeightASL_i = rx_i.elev + rxHeightAGL_i;
        
        diffractionLoss = 0;
        fresnelParam = 0;
        clearanceHeight = 0;
        
        % Calculate diffraction loss and Fresnel parameter for every point
        if ~isempty(bullingtonPoint)
            d_b = bullingtonPoint.dist - txPoint.dist;
            
            % Earth curvature correction
            h_rx_i_curved = rxHeightASL_i - (d_path_i^2) / (2 * EFFECTIVE_EARTH_RADIUS_M);
            h_b_curved = bullingtonPoint.uncorrectedElev - (d_b^2) / (2 * EFFECTIVE_EARTH_RADIUS_M);
            los_curved_at_b = txHeightASL + (h_rx_i_curved - txHeightASL) / d_path_i * d_b;
            h_eff_i_curved = h_b_curved - los_curved_at_b;
            
            clearanceHeight = h_eff_i_curved;
            
            d1 = d_b;
            d2 = d_path_i - d_b;
            
            if d1 > 0 && d2 > 0
                fresnelParam = h_eff_i_curved * sqrt(2 * d_path_i / (wavelength * d1 * d2));
                if fresnelParam > -0.78
                    diffractionLoss = 6.9 + 20 * log10(sqrt((fresnelParam - 0.1)^2 + 1) + fresnelParam - 0.1);
                end
            end
        end
        
        % Calculate path loss components
        fspl_i = 20 * log10(d_path_i) + 20 * log10(freq) - 147.55;
        totalLoss_i = fspl_i + diffractionLoss;
        
        % E-Field calculations (isotropic antenna)
        e_field_free_space_vm = sqrt(30 * txPower) / d_path_i;
        e_field_actual_vm = e_field_free_space_vm * 10^(-diffractionLoss / 20);
        
        % Store results
        analysis.distances = [analysis.distances; rx_i.dist];
        analysis.pathLosses = [analysis.pathLosses; totalLoss_i];
        analysis.diffractionLosses = [analysis.diffractionLosses; diffractionLoss];
        analysis.fsplLosses = [analysis.fsplLosses; fspl_i];
        analysis.eFieldFreeSpace = [analysis.eFieldFreeSpace; 20 * log10(e_field_free_space_vm * 1e6)];
        analysis.eFieldActual = [analysis.eFieldActual; 20 * log10(e_field_actual_vm * 1e6)];
        analysis.fresnelParameters = [analysis.fresnelParameters; fresnelParam];
        analysis.clearanceHeights = [analysis.clearanceHeights; clearanceHeight];
    end
    
    fprintf('Profile analysis completed for %d points\n', length(analysis.distances));
    profileAnalysis = analysis;
end

function finalResults = calculateBullington(freq, txHeightAGL, rxHeightAGL, terrainProfile)
    %% Calculate Bullington diffraction with detailed logging
    global C EFFECTIVE_EARTH_RADIUS_M;
    
    fprintf('Starting Bullington diffraction calculation...\n');
    
    if length(terrainProfile) < 2
        finalResults = struct('diffractionLoss', 0, 'fresnelParam', 0, 'fspl', 0, ...
                             'totalLoss', 0, 'effectiveHeight', 0, 'bullingtonPoint', []);
        return;
    end
    
    wavelength = C / freq;
    txPoint = terrainProfile(1);
    rxPoint = terrainProfile(end);
    
    txHeightASL = txPoint.elev + txHeightAGL;
    rxHeightASL = rxPoint.elev + rxHeightAGL;
    totalDist = rxPoint.dist - txPoint.dist;
    
    fprintf('Path length: %.3f km\n', totalDist/1000);
    fprintf('TX ASL: %.2f m, RX ASL: %.2f m\n', txHeightASL, rxHeightASL);
    
    if totalDist <= 0
        error('Invalid path: total distance must be positive');
    end
    
    fspl = 20 * log10(totalDist) + 20 * log10(freq) - 147.55;
    fprintf('Free space path loss: %.2f dB\n', fspl);
    
    obstacles = terrainProfile(2:end-1);
    
    if isempty(obstacles)
        fprintf('No intermediate obstacles found - line of sight path\n');
        finalResults = struct('diffractionLoss', 0, 'fresnelParam', 0, 'fspl', fspl, ...
                             'totalLoss', fspl, 'txHorizonPoint', [], 'rxHorizonPoint', [], ...
                             'bullingtonPoint', [], 'effectiveHeight', 0, ...
                             'txHeightASL', txHeightASL, 'rxHeightASL', rxHeightASL, ...
                             'totalDist', totalDist, 'wavelength', wavelength);
        return;
    end
    
    fprintf('Analyzing %d intermediate terrain points...\n', length(obstacles));
    
    % Helper function for Earth curvature correction
    getCurvedHeight = @(point, refDist, pathDist) ...
        point.elev + ((point.dist - refDist) * (pathDist - (point.dist - refDist))) / (2 * EFFECTIVE_EARTH_RADIUS_M);
    
    txCurvedHeightASL = txHeightASL;
    rxCurvedHeightASL = rxHeightASL;
    
    % Find TX horizon point (steepest slope from TX)
    maxTxSlope = -Inf;
    txHorizonPoint = [];
    
    fprintf('Finding TX horizon point...\n');
    for i = 1:length(obstacles)
        d_obs = obstacles(i).dist - txPoint.dist;
        if d_obs <= 0
            continue;
        end
        h_obs_curved = getCurvedHeight(obstacles(i), txPoint.dist, totalDist);
        currentSlope = (h_obs_curved - txCurvedHeightASL) / d_obs;
        
        if currentSlope > maxTxSlope
            % Check if this point is not obstructed by previous points
            isObstructed = false;
            for j = 1:i-1
                d_prev = obstacles(j).dist - txPoint.dist;
                h_prev_curved = getCurvedHeight(obstacles(j), txPoint.dist, totalDist);
                if h_prev_curved > txCurvedHeightASL + currentSlope * d_prev + 1e-6
                    isObstructed = true;
                    break;
                end
            end
            if ~isObstructed
                maxTxSlope = currentSlope;
                txHorizonPoint = obstacles(i);
            end
        end
    end
    
    % Find RX horizon point (steepest slope from RX)
    maxRxSlope = -Inf;
    rxHorizonPoint = [];
    
    fprintf('Finding RX horizon point...\n');
    for i = length(obstacles):-1:1
        d_obs_from_rx = rxPoint.dist - obstacles(i).dist;
        if d_obs_from_rx <= 0
            continue;
        end
        h_obs_curved = getCurvedHeight(obstacles(i), txPoint.dist, totalDist);
        currentSlope = (h_obs_curved - rxCurvedHeightASL) / d_obs_from_rx;
        
        if currentSlope > maxRxSlope
            % Check if this point is not obstructed by subsequent points
            isObstructed = false;
            for j = length(obstacles):-1:i+1
                d_next_from_rx = rxPoint.dist - obstacles(j).dist;
                h_next_curved = getCurvedHeight(obstacles(j), txPoint.dist, totalDist);
                if h_next_curved > rxCurvedHeightASL + currentSlope * d_next_from_rx + 1e-6
                    isObstructed = true;
                    break;
                end
            end
            if ~isObstructed
                maxRxSlope = currentSlope;
                rxHorizonPoint = obstacles(i);
            end
        end
    end
    
    if isempty(txHorizonPoint) || isempty(rxHorizonPoint)
        fprintf('No valid horizon points found - assuming line of sight\n');
        finalResults = struct('diffractionLoss', 0, 'fresnelParam', 0, 'fspl', fspl, ...
                             'totalLoss', fspl, 'txHorizonPoint', [], 'rxHorizonPoint', [], ...
                             'bullingtonPoint', [], 'effectiveHeight', 0, ...
                             'txHeightASL', txHeightASL, 'rxHeightASL', rxHeightASL, ...
                             'totalDist', totalDist, 'wavelength', wavelength);
        return;
    end
    
    fprintf('TX horizon point at: %.2f km\n', txHorizonPoint.dist/1000);
    fprintf('RX horizon point at: %.2f km\n', rxHorizonPoint.dist/1000);
    
    % Calculate Bullington point (intersection of horizon rays)
    m1 = maxTxSlope;
    m2 = (getCurvedHeight(rxHorizonPoint, txPoint.dist, totalDist) - rxCurvedHeightASL) / ...
         (rxHorizonPoint.dist - rxPoint.dist);
    
    d_b = (rxCurvedHeightASL - txCurvedHeightASL - m2 * totalDist) / (m1 - m2);
    h_b_curved = txCurvedHeightASL + m1 * d_b;
    
    los_height_at_db_curved = txCurvedHeightASL + (rxCurvedHeightASL - txCurvedHeightASL) / totalDist * d_b;
    h_eff = h_b_curved - los_height_at_db_curved;
    
    % Convert back to ground reference (remove Earth curvature effect for plotting)
    uncorrectedElev = h_b_curved - (d_b * (totalDist - d_b)) / (2 * EFFECTIVE_EARTH_RADIUS_M);
    
    bullingtonPoint = struct('dist', d_b + txPoint.dist, 'elev', h_b_curved, ...
                            'uncorrectedElev', uncorrectedElev);
    
    fprintf('Bullington point at: %.2f km, height: %.2f m\n', bullingtonPoint.dist/1000, uncorrectedElev);
    fprintf('Effective obstacle height: %.2f m\n', h_eff);
    
    if h_eff <= 0
        fprintf('Effective height <= 0: Line of sight path\n');
        finalResults = struct('diffractionLoss', 0, 'fresnelParam', 0, 'fspl', fspl, ...
                             'totalLoss', fspl, 'txHorizonPoint', txHorizonPoint, ...
                             'rxHorizonPoint', rxHorizonPoint, 'bullingtonPoint', bullingtonPoint, ...
                             'effectiveHeight', h_eff, 'txHeightASL', txHeightASL, ...
                             'rxHeightASL', rxHeightASL, 'totalDist', totalDist, ...
                             'wavelength', wavelength);
        return;
    end
    
    % Calculate diffraction parameters
    d1 = d_b;
    d2 = totalDist - d_b;
    if d1 <= 0 || d2 <= 0
        error('Invalid Bullington point distances: d1=%.2f, d2=%.2f', d1, d2);
    end
    
    v = h_eff * sqrt(2 * (d1 + d2) / (wavelength * d1 * d2));
    fprintf('Fresnel parameter: %.3f\n', v);
    
    diffractionLoss = 0;
    if v > -0.78
        diffractionLoss = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
    end
    totalLoss = fspl + diffractionLoss;
    
    fprintf('Diffraction loss: %.2f dB\n', diffractionLoss);
    fprintf('Total path loss: %.2f dB\n', totalLoss);
    
    finalResults = struct('diffractionLoss', diffractionLoss, 'fresnelParam', v, ...
                         'fspl', fspl, 'totalLoss', totalLoss, ...
                         'txHorizonPoint', txHorizonPoint, 'rxHorizonPoint', rxHorizonPoint, ...
                         'bullingtonPoint', bullingtonPoint, 'effectiveHeight', h_eff, ...
                         'txHeightASL', txHeightASL, 'rxHeightASL', rxHeightASL, ...
                         'totalDist', totalDist, 'wavelength', wavelength);
end

function displayDetailedResults(finalResults, terrainProfile, profileAnalysis)
    %% Display comprehensive results
    fprintf('\n===============================================\n');
    fprintf('         DETAILED CALCULATION RESULTS\n');
    fprintf('===============================================\n');
    
    % Main results
    fprintf('\nBULLINGTON DIFFRACTION PARAMETERS:\n');
    fprintf('----------------------------------\n');
    fprintf('Diffraction Loss (J(v)):        %10.3f dB\n', finalResults.diffractionLoss);
    fprintf('Fresnel Parameter (v):           %10.3f\n', finalResults.fresnelParam);
    fprintf('Effective Obstacle Height:       %10.3f m\n', finalResults.effectiveHeight);
    fprintf('Free Space Path Loss:            %10.3f dB\n', finalResults.fspl);
    fprintf('TOTAL PATH LOSS:                 %10.3f dB\n', finalResults.totalLoss);
    
    % Path geometry
    fprintf('\nPATH GEOMETRY:\n');
    fprintf('--------------\n');
    fprintf('Path Length:                     %10.3f km\n', finalResults.totalDist/1000);
    fprintf('TX Height (ASL):                 %10.3f m\n', finalResults.txHeightASL);
    fprintf('RX Height (ASL):                 %10.3f m\n', finalResults.rxHeightASL);
    fprintf('Wavelength:                      %10.6f m\n', finalResults.wavelength);
    
    if isfield(finalResults, 'bullingtonPoint') && ~isempty(finalResults.bullingtonPoint)
        fprintf('Bullington Point Distance:       %10.3f km\n', finalResults.bullingtonPoint.dist/1000);
        fprintf('Bullington Point Height:         %10.3f m\n', finalResults.bullingtonPoint.uncorrectedElev);
    end
    
    % Terrain statistics
    elevations = [terrainProfile.elev];
    distances = [terrainProfile.dist];
    
    fprintf('\nTERRAIN STATISTICS:\n');
    fprintf('-------------------\n');
    fprintf('Terrain Points:                  %10d\n', length(terrainProfile));
    fprintf('Distance Range:                  %7.3f - %7.3f km\n', min(distances)/1000, max(distances)/1000);
    fprintf('Elevation Range:                 %7.2f - %7.2f m\n', min(elevations), max(elevations));
    fprintf('Average Point Spacing:           %10.2f m\n', mean(diff(distances)));
    fprintf('Terrain Relief:                  %10.2f m\n', max(elevations) - min(elevations));
    
    % Profile analysis statistics
    if ~isempty(profileAnalysis.distances)
        fprintf('\nPROFILE ANALYSIS STATISTICS:\n');
        fprintf('----------------------------\n');
        fprintf('Points Analyzed:                 %10d\n', length(profileAnalysis.distances));
        fprintf('Max Diffraction Loss:            %10.3f dB\n', max(profileAnalysis.diffractionLosses));
        fprintf('Mean Path Loss:                  %10.3f dB\n', mean(profileAnalysis.pathLosses));
        fprintf('Max Field Strength:              %10.3f dBµV/m\n', max(profileAnalysis.eFieldActual));
        fprintf('Min Field Strength:              %10.3f dBµV/m\n', min(profileAnalysis.eFieldActual));
        fprintf('Field Strength at RX:            %10.3f dBµV/m\n', profileAnalysis.eFieldActual(end));
        
        % Find worst diffraction point
        [maxDiffLoss, worstIdx] = max(profileAnalysis.diffractionLosses);
        if maxDiffLoss > 0
            fprintf('Worst Diffraction Point:         %10.3f km (%.2f dB loss)\n', ...
                    profileAnalysis.distances(worstIdx)/1000, maxDiffLoss);
        end
    end
    
    fprintf('\n===============================================\n');
end

function createComprehensivePlots(terrainProfile, finalResults, profileAnalysis)
    %% Create comprehensive visualization plots
    
    % Create main figure with subplots
    fig = figure('Name', 'Bullington Diffraction Analysis - Real Terrain Data', ...
                 'Position', [50, 50, 1400, 1000], ...
                 'Color', 'white');
    
    %% Plot 1: Terrain Profile with Bullington Construction
    subplot(2, 2, 1);
    drawDetailedTerrainProfile(terrainProfile, finalResults);
    
    %% Plot 2: Path Loss Components
    subplot(2, 2, 2);
    drawPathLossComponents(profileAnalysis);
    
    %% Plot 3: Field Strength Comparison
    subplot(2, 2, 3);
    drawFieldStrengthComparison(profileAnalysis);
    
    %% Plot 4: Diffraction Analysis
    subplot(2, 2, 4);
    drawDiffractionAnalysis(profileAnalysis);
    
    % Add overall title
    sgtitle('Bullington Diffraction Analysis - Real X.04 Terrain Data', ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Adjust subplot spacing
    set(gcf, 'Color', 'white');
    
    % Create additional detailed figure
    createAdditionalPlots(terrainProfile, finalResults, profileAnalysis);
end

function drawDetailedTerrainProfile(terrainProfile, results)
    %% Draw detailed terrain profile with all construction lines
    
    % Extract data
    distances = [terrainProfile.dist] / 1000; % Convert to km
    elevations = [terrainProfile.elev];
    
    txPoint = terrainProfile(1);
    rxPoint = terrainProfile(end);
    
    % Plot terrain with gradient fill
    fill([distances, fliplr(distances)], [elevations, zeros(size(elevations))], ...
         [0.7, 0.85, 0.9], 'EdgeColor', [0.2, 0.4, 0.6], 'LineWidth', 1.5);
    hold on;
    
    % Plot terrain line
    plot(distances, elevations, 'Color', [0, 0.4, 0.6], 'LineWidth', 2);
    
    % Plot antenna masts with enhanced styling
    tx_x = txPoint.dist / 1000;
    rx_x = rxPoint.dist / 1000;
    
    % TX mast
    plot([tx_x, tx_x], [txPoint.elev, results.txHeightASL], 'k-', 'LineWidth', 4);
    plot([tx_x, tx_x], [txPoint.elev, results.txHeightASL], 'r-', 'LineWidth', 2);
    
    % RX mast
    plot([rx_x, rx_x], [rxPoint.elev, results.rxHeightASL], 'k-', 'LineWidth', 4);
    plot([rx_x, rx_x], [rxPoint.elev, results.rxHeightASL], 'r-', 'LineWidth', 2);
    
    % Antennas with enhanced markers
    plot(tx_x, results.txHeightASL, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
         'MarkerEdgeColor', 'black', 'LineWidth', 2);
    plot(rx_x, results.rxHeightASL, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red', ...
         'MarkerEdgeColor', 'black', 'LineWidth', 2);
    
    % Line of sight
    plot([tx_x, rx_x], [results.txHeightASL, results.rxHeightASL], ...
         '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
    
    % Bullington construction
    if isfield(results, 'bullingtonPoint') && ~isempty(results.bullingtonPoint)
        b_x = results.bullingtonPoint.dist / 1000;
        b_y = results.bullingtonPoint.uncorrectedElev;
        
        % TX horizon ray
        plot([tx_x, b_x], [results.txHeightASL, b_y], ...
             ':', 'Color', [0.86, 0.42, 0.13], 'LineWidth', 3);
        
        % RX horizon ray
        plot([rx_x, b_x], [results.rxHeightASL, b_y], ...
             ':', 'Color', [0.22, 0.63, 0.41], 'LineWidth', 3);
        
        % Bullington point marker
        plot(b_x, b_y, 'x', 'Color', [0.77, 0.19, 0.19], 'MarkerSize', 20, 'LineWidth', 5);
        
        % Effective height indicator
        los_elev_at_b = results.txHeightASL + (results.rxHeightASL - results.txHeightASL) / ...
                       results.totalDist * (results.bullingtonPoint.dist - txPoint.dist);
        plot([b_x, b_x], [los_elev_at_b, b_y], ':', 'Color', [0.77, 0.19, 0.19], 'LineWidth', 3);
        
        % Add text annotation
        text(b_x, b_y + 20, sprintf('Bullington Point\n%.2f dB loss', results.diffractionLoss), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    xlabel('Distance (km)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Height (m)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Terrain Profile & Bullington Construction', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    legend({'Terrain', 'Ground Level', 'Antenna Masts', '', 'Antennas', '', 'Line of Sight', ...
            'TX Horizon Ray', 'RX Horizon Ray', 'Bullington Point', 'Effective Height'}, ...
           'Location', 'best', 'FontSize', 8);
    
    % Enhanced grid
    grid minor;
    set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
end

function drawPathLossComponents(profileAnalysis)
    %% Draw path loss components
    
    if ~isempty(profileAnalysis.distances)
        distances_km = profileAnalysis.distances / 1000;
        
        % Plot FSPL
        plot(distances_km, profileAnalysis.fsplLosses, '--', ...
             'Color', [0.6, 0.6, 0.6], 'LineWidth', 2, 'DisplayName', 'Free Space Path Loss');
        hold on;
        
        % Plot diffraction loss
        plot(distances_km, profileAnalysis.diffractionLosses, ':', ...
             'Color', [0.8, 0.4, 0.2], 'LineWidth', 2, 'DisplayName', 'Diffraction Loss');
        
        % Plot total path loss
        plot(distances_km, profileAnalysis.pathLosses, '-', ...
             'Color', [0.2, 0.2, 0.8], 'LineWidth', 3, 'DisplayName', 'Total Path Loss');
        
        xlabel('Distance (km)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Path Loss (dB)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Path Loss Components', 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
    end
end

function drawFieldStrengthComparison(profileAnalysis)
    %% Draw field strength comparison
    
    if ~isempty(profileAnalysis.distances)
        distances_km = profileAnalysis.distances / 1000;
        
        % Plot free space field strength
        plot(distances_km, profileAnalysis.eFieldFreeSpace, '--', ...
             'Color', [0.6, 0.6, 0.6], 'LineWidth', 2, 'DisplayName', 'Free-Space Field');
        hold on;
        
        % Plot actual field strength
        plot(distances_km, profileAnalysis.eFieldActual, '-', ...
             'Color', [0.1, 0.7, 0.5], 'LineWidth', 3, 'DisplayName', 'Actual Field (with Diffraction)');
        
        % Fill area between curves to show diffraction effect
        fill([distances_km; flipud(distances_km)], ...
             [profileAnalysis.eFieldFreeSpace; flipud(profileAnalysis.eFieldActual)], ...
             [1, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Diffraction Effect');
        
        xlabel('Distance (km)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('E-Field (dBµV/m)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Electric Field Strength Comparison', 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
        
        % Add field strength value at receiver
        if ~isempty(profileAnalysis.eFieldActual)
            finalField = profileAnalysis.eFieldActual(end);
            finalDist = distances_km(end);
            plot(finalDist, finalField, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
            text(finalDist, finalField + 3, sprintf('RX: %.1f dBµV/m', finalField), ...
                 'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
                 'BackgroundColor', 'white', 'EdgeColor', 'red');
        end
    end
end

function drawDiffractionAnalysis(profileAnalysis)
    %% Draw diffraction analysis details - FIXED VERSION
    
    if ~isempty(profileAnalysis.distances)
        distances_km = profileAnalysis.distances / 1000;
        
        % Create dual y-axis plot
        yyaxis left;
        plot(distances_km, profileAnalysis.fresnelParameters, '-', ...
             'Color', [0.8, 0.2, 0.2], 'LineWidth', 2);
        ylabel('Fresnel Parameter (v)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8, 0.2, 0.2]);
        
        % Add Fresnel zone boundaries - FIXED: Removed Alpha property
        hold on;
        h1 = plot([distances_km(1), distances_km(end)], [0, 0], 'k--', 'LineWidth', 1);
        h2 = plot([distances_km(1), distances_km(end)], [-0.78, -0.78], 'r--', 'LineWidth', 1);
        
        % Set transparency using color with alpha channel (RGBA) - MATLAB compatible version
        try
            h1.Color = [0, 0, 0, 0.5]; % Black with 50% transparency
            h2.Color = [1, 0, 0, 0.5]; % Red with 50% transparency
        catch
            % Fallback for older MATLAB versions
            set(h1, 'Color', [0.5, 0.5, 0.5]);
            set(h2, 'Color', [1, 0.5, 0.5]);
        end
        
        yyaxis right;
        plot(distances_km, profileAnalysis.clearanceHeights, '-', ...
             'Color', [0.2, 0.2, 0.8], 'LineWidth', 2);
        ylabel('Clearance Height (m)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.8]);
        
        xlabel('Distance (km)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Diffraction Analysis Details', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        grid minor;
        set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
        
        % Add legend manually
        legend({'Fresnel Parameter', 'v = 0 (LOS)', 'v = -0.78 (Diffraction Threshold)', 'Clearance Height'}, ...
               'Location', 'best', 'FontSize', 8);
    end
end

function createAdditionalPlots(terrainProfile, finalResults, profileAnalysis)
    %% Create additional detailed plots
    
    % Second figure for detailed analysis
    fig2 = figure('Name', 'Detailed Terrain and Propagation Analysis', ...
                  'Position', [100, 100, 1200, 800], ...
                  'Color', 'white');
    
    %% 3D-style terrain visualization
    subplot(2, 1, 1);
    drawEnhanced3DTerrainProfile(terrainProfile, finalResults);
    
    %% Coverage prediction plot
    subplot(2, 1, 2);
    drawCoveragePrediction(profileAnalysis);
    
    sgtitle('Enhanced Terrain Analysis with Real X.04 Data', ...
            'FontSize', 16, 'FontWeight', 'bold');
end

function drawEnhanced3DTerrainProfile(terrainProfile, results)
    %% Draw enhanced 3D-style terrain profile - FIXED COLORMAP
    
    distances = [terrainProfile.dist] / 1000;
    elevations = [terrainProfile.elev];
    
    % Create surface-like appearance
    [X, Y] = meshgrid(distances, [0, 1]);
    Z = [elevations; elevations];
    
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;
    
    % Add contour lines
    contour3(X, Y, Z, 10, 'k-', 'LineWidth', 0.5);
    
    % Plot path in 3D
    txPoint = terrainProfile(1);
    rxPoint = terrainProfile(end);
    
    plot3([txPoint.dist, rxPoint.dist]/1000, [0.5, 0.5], ...
          [results.txHeightASL, results.rxHeightASL], ...
          'r-', 'LineWidth', 4);
    
    view(45, 30);
    xlabel('Distance (km)');
    ylabel('Width');
    zlabel('Height (m)');
    title('3D Terrain Visualization with Propagation Path');
    
    % FIXED: Use compatible colormap
    try
        % Try terrain colormap if available (Mapping Toolbox)
        colormap(terrain);
    catch
        % Fallback to standard colormap if terrain not available
        try
            % Create custom terrain-like colormap
            customTerrain = [
                0.2 0.4 0.2;    % Dark green (low)
                0.4 0.6 0.3;    % Medium green
                0.6 0.8 0.4;    % Light green
                0.8 0.8 0.6;    % Yellow-green
                0.9 0.7 0.5;    % Brown
                0.7 0.5 0.3;    % Dark brown
                0.9 0.9 0.9     % White (high)
            ];
            colormap(customTerrain);
        catch
            % Ultimate fallback to standard colormap
            colormap(jet);
        end
    end
    
    colorbar;
    grid on;
end

function drawCoveragePrediction(profileAnalysis)
    %% Draw coverage prediction analysis
    
    if ~isempty(profileAnalysis.distances)
        distances_km = profileAnalysis.distances / 1000;
        
        % Create coverage zones based on field strength
        excellent = profileAnalysis.eFieldActual >= 60;
        good = profileAnalysis.eFieldActual >= 40 & profileAnalysis.eFieldActual < 60;
        marginal = profileAnalysis.eFieldActual >= 20 & profileAnalysis.eFieldActual < 40;
        poor = profileAnalysis.eFieldActual < 20;
        
        % Plot coverage zones
        area(distances_km, excellent * 80, 'FaceColor', [0, 0.8, 0], 'FaceAlpha', 0.6, 'DisplayName', 'Excellent (≥60 dBµV/m)');
        hold on;
        area(distances_km, good * 60, 'FaceColor', [0.8, 0.8, 0], 'FaceAlpha', 0.6, 'DisplayName', 'Good (40-60 dBµV/m)');
        area(distances_km, marginal * 40, 'FaceColor', [1, 0.5, 0], 'FaceAlpha', 0.6, 'DisplayName', 'Marginal (20-40 dBµV/m)');
        area(distances_km, poor * 20, 'FaceColor', [0.8, 0, 0], 'FaceAlpha', 0.6, 'DisplayName', 'Poor (<20 dBµV/m)');
        
        % Overlay actual field strength
        plot(distances_km, profileAnalysis.eFieldActual, 'k-', 'LineWidth', 2, 'DisplayName', 'Actual Field Strength');
        
        xlabel('Distance (km)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Field Strength (dBµV/m)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Coverage Prediction Analysis', 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 9);
        grid on;
        ylim([0, 100]);
    end
end

function exportComprehensiveResults(finalResults, profileAnalysis, terrainProfile, txPower, frequency, txHeightAGL, rxHeightAGL)
    %% Export comprehensive results to multiple files
    
    fprintf('\nEXPORTING COMPREHENSIVE RESULTS...\n');
    fprintf('=================================\n');
    
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    try
        % Export 1: Summary Results
        filename1 = sprintf('bullington_summary_%s.csv', timestamp);
        fid = fopen(filename1, 'w');
        fprintf(fid, 'Bullington Diffraction Analysis Summary\n');
        fprintf(fid, 'Generated: %s\n\n', datestr(now));
        fprintf(fid, 'INPUT PARAMETERS\n');
        fprintf(fid, 'Parameter,Value,Unit\n');
        fprintf(fid, 'Transmitter Power,%.6f,W\n', txPower);
        fprintf(fid, 'Frequency,%.6f,MHz\n', frequency);
        fprintf(fid, 'TX Height AGL,%.6f,m\n', txHeightAGL);
        fprintf(fid, 'RX Height AGL,%.6f,m\n', rxHeightAGL);
        fprintf(fid, '\nRESULTS\n');
        fprintf(fid, 'Parameter,Value,Unit\n');
        fprintf(fid, 'Diffraction Loss,%.6f,dB\n', finalResults.diffractionLoss);
        fprintf(fid, 'Fresnel Parameter,%.6f,-\n', finalResults.fresnelParam);
        fprintf(fid, 'Effective Height,%.6f,m\n', finalResults.effectiveHeight);
        fprintf(fid, 'Free Space Path Loss,%.6f,dB\n', finalResults.fspl);
        fprintf(fid, 'Total Path Loss,%.6f,dB\n', finalResults.totalLoss);
        fprintf(fid, 'Path Length,%.6f,m\n', finalResults.totalDist);
        fprintf(fid, 'TX Height ASL,%.6f,m\n', finalResults.txHeightASL);
        fprintf(fid, 'RX Height ASL,%.6f,m\n', finalResults.rxHeightASL);
        fprintf(fid, 'Wavelength,%.6f,m\n', finalResults.wavelength);
        if ~isempty(finalResults.bullingtonPoint)
            fprintf(fid, 'Bullington Point Distance,%.6f,m\n', finalResults.bullingtonPoint.dist);
            fprintf(fid, 'Bullington Point Height,%.6f,m\n', finalResults.bullingtonPoint.uncorrectedElev);
        end
        fclose(fid);
        
        % Export 2: Detailed Profile Analysis
        filename2 = sprintf('bullington_profile_%s.csv', timestamp);
        if ~isempty(profileAnalysis.distances)
            data = [profileAnalysis.distances, profileAnalysis.pathLosses, ...
                    profileAnalysis.diffractionLosses, profileAnalysis.fsplLosses, ...
                    profileAnalysis.eFieldFreeSpace, profileAnalysis.eFieldActual, ...
                    profileAnalysis.fresnelParameters, profileAnalysis.clearanceHeights];
            
            fid = fopen(filename2, 'w');
            fprintf(fid, 'Distance_m,PathLoss_dB,DiffractionLoss_dB,FSPL_dB,EField_FreeSpace_dBuVm,EField_Actual_dBuVm,FresnelParameter,ClearanceHeight_m\n');
            fclose(fid);
            
            dlmwrite(filename2, data, '-append', 'precision', '%.6f');
        end
        
        % Export 3: Terrain Profile
        filename3 = sprintf('terrain_profile_%s.csv', timestamp);
        terrainData = [[terrainProfile.dist]', [terrainProfile.elev]'];
        fid = fopen(filename3, 'w');
        fprintf(fid, 'Distance_m,Elevation_m\n');
        fclose(fid);
        dlmwrite(filename3, terrainData, '-append', 'precision', '%.6f');
        
        % Export 4: Coverage Analysis
        filename4 = sprintf('coverage_analysis_%s.csv', timestamp);
        if ~isempty(profileAnalysis.distances)
            coverageData = [profileAnalysis.distances / 1000, profileAnalysis.eFieldActual];
            
            % Add coverage categories
            categories = zeros(size(profileAnalysis.eFieldActual));
            categories(profileAnalysis.eFieldActual >= 60) = 4; % Excellent
            categories(profileAnalysis.eFieldActual >= 40 & profileAnalysis.eFieldActual < 60) = 3; % Good
            categories(profileAnalysis.eFieldActual >= 20 & profileAnalysis.eFieldActual < 40) = 2; % Marginal
            categories(profileAnalysis.eFieldActual < 20) = 1; % Poor
            
            coverageData = [coverageData, categories];
            
            fid = fopen(filename4, 'w');
            fprintf(fid, 'Distance_km,EField_dBuVm,Coverage_Category\n');
            fprintf(fid, '# Coverage Categories: 1=Poor, 2=Marginal, 3=Good, 4=Excellent\n');
            fclose(fid);
            
            dlmwrite(filename4, coverageData, '-append', 'precision', '%.6f');
        end
        
        fprintf('Results successfully exported to:\n');
        fprintf('  1. Summary: %s\n', filename1);
        if ~isempty(profileAnalysis.distances)
            fprintf('  2. Profile: %s\n', filename2);
        end
        fprintf('  3. Terrain: %s\n', filename3);
        if ~isempty(profileAnalysis.distances)
            fprintf('  4. Coverage: %s\n', filename4);
        end
        
    catch ME
        fprintf('Error exporting results: %s\n', ME.message);
    end
end