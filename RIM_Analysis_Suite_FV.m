%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIM Analysis Suite
%
% A MATLAB toolbox for Real & Imaginary Mapping (RIM) analysis, designed to
% compute and visualize spatial impedance maps using synchronized optical
% and electrical signals.
%
% -------------------------------------------------------------------------------------
% MIT License
%
% Copyright (c) 2026 Giagu Gabriele / Electrochemistry of Molecular and Functional Materials group headed by Prof. Francesco Paolucci / Unibo, Department of Chemistry "G.Ciamician" Via Gobetti 85, Bologna (BO)
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sub license, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% -------------------------------------------------------------------------------------
%
% Version: 1.0.0
% System requirements: MATLAB R2019b or above, Image Processing Toolbox
%
% How to cite this work:
% If you use this toolbox in your research, please cite:
% [1] Giagu Gabriele, Unraveling the Electrode-Reaction Layer in Electrochemiluminescence by Reflective Impedance Microscopy, Chemical and Biomedical Imaging, 2026.
% Link: [DOI or URL of the paper once published]
%
% For support or bug reports, please contact: gabriele.giagu2@unibo.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===================================================================
% MAIN SCRIPT - RIM Analysis Suite
% Entry point. Mode selection only. Parameters moved to specific modes.
% ===================================================================

% Retain persistent variables across script executions
clearvars -except persistentScriptDefaultPath;
close all; clc;

%% --- 1. Initial Settings ---
disp('--- Starting RIM Analysis ---');
if exist('persistentScriptDefaultPath','var') == 0
    persistentScriptDefaultPath = pwd;
end

%% --- 2. Mode Selection ---
listOptions = {'Calibration', 'Measurement', 'Real & Imaginary Mapping'};
[choice, tf] = listdlg('ListString', listOptions, ...
    'SelectionMode', 'single', ...
    'PromptString', 'Select RIM Analysis Mode:', ...
    'Name', 'RIM Mode Selection', ...
    'ListSize', [250, 100]);

if tf == 0 || isempty(choice)
    disp('Operation cancelled.');
    return;
end

%% --- 3. Execution ---
try
    if choice == 1 % Calibration
        persistentScriptDefaultPath = RIM_Calibration(persistentScriptDefaultPath);
    elseif choice == 2 % Measurement
        persistentScriptDefaultPath = RIM_Measurement(persistentScriptDefaultPath);
    elseif choice == 3 % Real & Imaginary Mapping
        persistentScriptDefaultPath = RIM_RedoxAnalysis(persistentScriptDefaultPath);
    end
catch ME
    errordlg(sprintf('An error occurred: %s', ME.message), 'Critical Error');
    fprintf(2, 'Error in %s (line %d): %s\n', ME.stack(1).name, ME.stack(1).line, ME.message);
end

disp('--- Script Finished ---');


%% ===================================================================
%  MAIN MODULE FUNCTIONS
%  ===================================================================

function updatedPath = RIM_Calibration(startPath)
% Executes the single Beta calibration workflow.
updatedPath = startPath;
calibFile = 'rim_calib_params.mat';

disp('--- Starting Single Beta Calibration ---');

analysisChoice = questdlg('Select the Signal Extraction Method:', ...
    'Analysis Method', ...
    'Standard FFT (Linear Sine)', 'Faradaic Matched Filter (Non-Linear Q)', 'Standard FFT (Linear Sine)');

if isempty(analysisChoice)
    disp('Calibration cancelled.'); return;
end
useMatchedFilter = (string(analysisChoice) == "Faradaic Matched Filter (Non-Linear Q)");

[f_ref, fr_int, applyFilter, sigmaValue, ~, electrodeArea] = utils_promptAnalysisParams('Calibration');
if isempty(f_ref) || isnan(f_ref)
    disp('Calibration cancelled.');
    return;
end

[tbl, updatedPath] = createInteractiveTableFolders(updatedPath, useMatchedFilter);
if isempty(tbl)
    disp('Calibration cancelled.'); return;
end

numPoints = height(tbl);

applyCircuit = false;
circuitFunc = [];
circParamsToSave = [];

choice = questdlg('Apply Equivalent Circuit correction to Calibration Voltages?', 'Circuit Builder', 'Yes', 'No', 'No');
if string(choice) == "Yes"
    dummyV = tbl.ActualVoltage(1);
    if isnan(dummyV)
        dummyV = 0.4;
    end

    circParamsIn = [];
    loadChoice = questdlg('Would you like to load a previously used Equivalent Circuit from an existing Calibration file?', 'Load Circuit', 'Yes', 'No', 'No');

    if string(loadChoice) == "Yes"
        [mat_file, mat_path] = uigetfile('*.mat', 'Select Calibration MAT File', updatedPath);
        if ~isequal(mat_file, 0)
            try
                loadedData = load(fullfile(mat_path, mat_file));
                if isfield(loadedData, 'circParamsToSave') && ~isempty(loadedData.circParamsToSave)
                    circParamsIn = loadedData.circParamsToSave;
                    disp(['>> Successfully loaded previous Equivalent Circuit parameters from: ', mat_file]);
                else
                    warndlg('No circuit parameters found in the selected file. Opening default builder.', 'Missing Parameters');
                end
            catch ME
                warndlg(sprintf('Failed to load file:\n%s', ME.message), 'Load Error');
            end
        end
    end

    [~, circuitFunc, canc_circ, C_eff_func, circParamsToSave] = RIM_CircuitBuilder(dummyV, f_ref, circParamsIn);
    if canc_circ
        disp('Calibration cancelled during Circuit Builder.'); return;
    end
    applyCircuit = true;
else
    ansC = inputdlg({'Enter the single Reference Capacitance (F) used for all these measurements:'}, ...
        'Global Reference Capacitance', [1 60], {'1e-5'});
    if isempty(ansC), disp('Calibration cancelled.'); return; end
    global_C_val = str2double(ansC{1});
    if isnan(global_C_val), errordlg('Invalid capacitance value.', 'Input Error'); return; end
    C_eff_func = @(f) global_C_val;
end

V_export = NaN(numPoints, 1);
C_export = NaN(numPoints, 1);
dR_export = NaN(numPoints, 1);

roiMaskBulk = [];

for k = 1:numPoints
    fprintf('\n--- Processing Point %d/%d ---\n', k, numPoints);

    folderPath = char(tbl.Folder(k));
    V_act = tbl.ActualVoltage(k);
    lasv_file = char(tbl.LASVFile(k));
    lasv_path = char(tbl.LASVPath(k));
    f_ref_opt = f_ref;

    has_lasv_for_point = false;
    if ~strcmp(lasv_file, 'None') && ~isempty(lasv_file)
        try
            [t_el, V_el, I_el] = utils_loadECLabLASV(fullfile(lasv_path, lasv_file));
            dt_el = mean(diff(t_el));

            search_radius = max(0.5, 0.2 * f_ref);
            [V_amp, ~, f_peak_el] = math_findMaxInFreqInterval(V_el, dt_el, max(0.1, f_ref - search_radius), f_ref + search_radius);
            canc_el = false;

            if ~canc_el && ~isnan(V_amp)
                fprintf('LASV Reference Loaded: True V_amp = %.4e V at %.3f Hz (Auto-FFT, Replaces Nominal %.3f V)\n', V_amp, f_peak_el, V_act);
                V_act = V_amp;
                f_ref_opt = f_peak_el;
                has_lasv_for_point = true;

                if useMatchedFilter
                    Q_el = cumtrapz(t_el, I_el);
                    [h_amps, h_phases] = math_extractHarmonicSignature(Q_el, dt_el, f_peak_el, 5);
                    fprintf('Point %d Harmonic Signature:\n', k);
                    disp(h_amps);
                end
            end
        catch ME
            warning('RIM_Calibration:LASVError', 'Failed to load LASV file: %s. Falling back to nominal.', ME.message);
        end
    end

    lbl_C = C_eff_func(f_ref_opt);
    lbl = sprintf('%.3fV, C_eff=%.2eF', V_act, lbl_C);

    files = dir(fullfile(folderPath, '*.tif'));
    if isempty(files), files = dir(fullfile(folderPath, '*.tiff')); end
    if isempty(files), warning('No TIF files found in %s', folderPath); continue; end

    fileNames = utils_naturalSort({files.name});
    numFrames = length(fileNames);
    info = imfinfo(fullfile(folderPath, fileNames{1}));
    img = zeros(info.Height, info.Width, numFrames, 'single');

    hW = waitbar(0, sprintf('Loading %d TIFFs...', numFrames));
    try
        for fIdx = 1:numFrames
            if ~isvalid(hW), error('User Cancelled Loading'); end
            img(:,:,fIdx) = single(imread(fullfile(folderPath, fileNames{fIdx})));
            if mod(fIdx, 20) == 0 && isvalid(hW), waitbar(fIdx/numFrames, hW); end
        end
    catch ME
        if isvalid(hW), close(hW); end
        warning('RIM_Calibration:LoadError', 'Error loading TIFFs: %s', ME.message);
        continue;
    end
    if isvalid(hW), close(hW); end

    if applyFilter
        img = math_applySpatialFilter(img, sigmaValue);
    end

    if isempty(roiMaskBulk)
        [roiMaskBulk, ~] = utils_selectROI(img(:,:,1), ['Setup Bulk ROI - ' lbl]);
        if isempty(roiMaskBulk), disp('Cancelled.'); break; end
    end
    ts_full = math_calculateROIAverage(img, roiMaskBulk);
    [ts_win, ~, ~] = utils_selectWindowFromSpectrogram(ts_full, fr_int, f_ref_opt, lbl);
    if isempty(ts_win), continue; end

    if useMatchedFilter && has_lasv_for_point
        fprintf('>>> Point %d: Extracting signal using Custom Non-Linear Lock-In (Synthetic Matched Filter)...\n', k);
        [~, ~, f_peak_opt, canc] = utils_selectFFTPeak(ts_win, fr_int, f_ref_opt, [lbl ' (Optical Sync)']);
        if canc || isnan(f_peak_opt), continue; end

        numFrames_win = length(ts_win);
        [K_I, K_Q] = math_buildSyntheticKernel(numFrames_win, fr_int, f_peak_opt, h_amps, h_phases);
        [dR_abs, ~] = math_calculateMatchedFilter1D(ts_win, K_I, K_Q);
        freq_found = f_peak_opt;
    else
        fprintf('>>> Point %d: Extracting signal using Standard FFT...\n', k);
        [dR_abs, ~, freq_found, canc] = utils_selectFFTPeak(ts_win, fr_int, f_ref_opt, lbl);
        if canc || isnan(freq_found), continue; end
    end

    R_base = mean(ts_win, 'omitnan');
    if R_base == 0, R_base = 1; end
    dR = dR_abs / R_base;

    dR_export(k) = dR ;

    current_C_eff = C_eff_func(freq_found);

    if applyCircuit
        V_act = abs(circuitFunc(V_act, freq_found));
    end

    V_export(k) = V_act;
    C_export(k) = current_C_eff;
end

valid = ~isnan(V_export);
V_export = V_export(valid); C_export = C_export(valid);

C_export = C_export / electrodeArea;

if sum(valid) >= 1
    dR_export = dR_export(valid);
    [V_clean, dR_clean, mask] = utils_excludeOutliers(V_export, dR_export);
    if sum(mask) > 0
        C_clean = C_export(mask);
        [betaValue, qValue] = math_calculateBeta(V_clean, dR_clean, C_clean);
        calibType = 'Single';

        save(calibFile, 'calibType', 'betaValue', 'qValue', 'circParamsToSave');

        plotCalibrationFit(C_clean.*V_clean, dR_clean, betaValue, qValue, 'Bulk Fit');
        msgbox(sprintf('Calibration Complete.\nBeta = %.4e\nIntercept = %.4e', betaValue, qValue), 'Success');

        choice_export = questdlg('Export a portable copy of the Calibration data (.mat) for future measurements?', 'Export Calibration', 'Yes', 'No', 'Yes');
        if string(choice_export) == "Yes"
            [fName, pName] = uiputfile('*.mat', 'Save Calibration Data', fullfile(updatedPath, 'Calibration_Data.mat'));
            if ~isequal(fName, 0)
                try
                    save(fullfile(pName, fName), 'calibType', 'betaValue', 'qValue', 'circParamsToSave');
                    msgbox('Portable calibration data saved successfully.', 'Export Complete');
                catch ME
                    errordlg(sprintf('An error occurred while saving the file:\n%s', ME.message), 'Export Error');
                end
            end
        end
    else
        errordlg('All points excluded.', 'Error');
    end
else
    errordlg('Not enough valid points for calibration.', 'Error');
end
end

function plotCalibrationFit(X, Y, beta, q, tLabel)
% Plots the linear fit calibration data.
figure('Name', 'Calibration Fit Result', 'Position', [200, 200, 700, 500]);
ax = axes; hold on;
plot(ax, X, Y, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
x_line = linspace(0, max(X)*1.1, 100);
plot(ax, x_line, beta * x_line + q, 'r-', 'LineWidth', 2);
xlabel(ax, 'C_{specific} \times V (F/cm^2 \cdot V)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, '\Delta R / R', 'FontSize', 12, 'FontWeight', 'bold');
title(ax, sprintf('%s\n\\beta = %.4e, q = %.4e', tLabel, beta, q), 'FontSize', 14);
legend(ax, 'Data Points', 'Linear Fit (y = \beta x + q)', 'Location', 'best');
grid on; hold off;
end

function [tbl, updatedPath] = createInteractiveTableFolders(startPath, useMatchedFilter)
% Manages the interactive folder selection and matched LASV association for bulk calibration.
tbl = table();
updatedPath = startPath;
folders = {};

choice_method = questdlg('How would you like to select the calibration folders?', ...
    'Folder Selection Method', ...
    'Select Parent Folder (Auto-Detect)', 'Select Manually One-by-One', 'Select Parent Folder (Auto-Detect)');

if string(choice_method) == "Select Parent Folder (Auto-Detect)"
    parentDir = uigetdir(updatedPath, 'Select the PARENT folder containing all Calibration subfolders');
    if isequal(parentDir, 0), return; end
    updatedPath = parentDir;

    d = dir(parentDir);
    isub = [d(:).isdir];
    nameFolds = {d(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];

    validFolders = {};
    for i = 1:length(nameFolds)
        subPath = fullfile(parentDir, nameFolds{i});
        tifs = dir(fullfile(subPath, '*.tif'));
        if isempty(tifs), tifs = dir(fullfile(subPath, '*.tiff')); end
        if ~isempty(tifs)
            validFolders{end+1} = subPath; %#ok<AGROW>
        end
    end
    folders = utils_naturalSort(validFolders);

elseif string(choice_method) == "Select Manually One-by-One"
    while true
        choice = questdlg('Add a TIFF Folder to Calibration?', 'Add Folder', 'Yes', 'No (Done)', 'Yes');
        if string(choice) == "No (Done)" || isempty(choice)
            break;
        end
        pn = uigetdir(updatedPath, 'Select TIFF Folder for Calibration Point');
        if ~isequal(pn, 0)
            folders{end+1} = pn; %#ok<AGROW>
            updatedPath = pn;
        end
    end
else
    return;
end

N = length(folders);
if N == 0
    errordlg('No valid TIFF folders selected or found.', 'Folder Error');
    return;
end

if useMatchedFilter
    choice_lasv = "Yes";
else
    choice_lasv = questdlg('Load EC-Lab LASV files for ALL points to extract the TRUE Voltage Amplitude?', 'LASV Bulk Selection', 'Yes', 'No', 'Yes');
end

lasv_list = repmat({'None'}, N, 1);
lasv_dir = '';

if string(choice_lasv) == "Yes"
    [lasv_files, lasv_dir] = uigetfile({'*.mpt;*.txt;*.csv', 'EC-Lab Data Files'}, ...
        'Select ALL corresponding LASV Files at once', updatedPath, 'MultiSelect', 'on');

    if ~isequal(lasv_files, 0)
        if ischar(lasv_files)
            lasv_files = {lasv_files};
        end

        lasv_files = utils_naturalSort(lasv_files);

        if length(lasv_files) ~= N
            warndlg(sprintf('You selected %d LASV files, but there are %d folders. They will be matched sequentially from the top.', length(lasv_files), N), 'Mismatch Warning');
        end

        for i = 1:min(N, length(lasv_files))
            lasv_list{i} = lasv_files{i};
        end
        updatedPath = lasv_dir;
    end
end

d = cell(N, 3);
d(:,1) = folders';
d(:,2) = num2cell(NaN(N, 1));
d(:,3) = lasv_list;

fig = uifigure('Name', 'Verify Subfolders, Voltage & LASV Matches', 'Position', [100 100 1000 400], 'WindowStyle', 'modal');
gl = uigridlayout(fig, [2 1]);
gl.RowHeight = {'1x', 40};

t = uitable(gl, 'Data', d, 'ColumnName', {'Folder', 'Nominal Delta V (V)', 'Matched LASV File'}, ...
    'ColumnEditable', [false, true, false], 'ColumnWidth', {'2x', 'auto', '1x'});
t.Layout.Row = 1;
t.Layout.Column = 1;

btn = uibutton(gl, 'Text', 'CONFIRM MATCHES & CONTINUE', ...
    'FontWeight', 'bold', 'BackgroundColor', '#77AC30', 'FontColor', 'w', ...
    'ButtonPushedFcn', @(src, event) uiresume(fig));
btn.Layout.Row = 2;
btn.Layout.Column = 1;

uiwait(fig);
if ~isvalid(fig), return; end
finalData = t.Data; close(fig);

ActualVoltage = cell2mat(finalData(:,2));
LASVFile = finalData(:,3);
LASVPath = repmat({lasv_dir}, N, 1);

for i = 1:N
    if useMatchedFilter || isnan(ActualVoltage(i))
        if strcmp(LASVFile{i}, 'None') || isempty(LASVFile{i})
            if useMatchedFilter
                msg = sprintf('Matched Filter active. LASV file is REQUIRED for Folder:\n%s', char(finalData{i,1}));
            else
                msg = sprintf('Nominal Voltage is empty for Folder:\n%s\n\nPlease select an LASV file to extract the true voltage automatically.', char(finalData{i,1}));
            end
            uiwait(msgbox(msg, 'Missing LASV/Voltage Required', 'warn'));

            [fName, pName] = uigetfile({'*.mpt;*.txt;*.csv', 'EC-Lab Data Files'}, 'Select LASV File', lasv_dir);
            if isequal(fName, 0)
                errordlg('Operation cancelled. LASV file is required.', 'Input Error');
                tbl = table(); return;
            end
            LASVFile{i} = fName;
            LASVPath{i} = pName;
            lasv_dir = pName;
        end
    end
end

Folder = finalData(:,1);
tbl = table(Folder, ActualVoltage, LASVFile, LASVPath);
end

function [V_eff, evalFunc, cancelled, C_eff_func, circParamsOut] = RIM_CircuitBuilder(V_preview, f_preview, circParamsIn)
% RIM_CIRCUITBUILDER Interactive UI to calculate voltage drop across a specific
% component in an equivalent circuit using the voltage divider principle.
if nargin < 3, circParamsIn = []; end
if nargin < 2
    V_preview = NaN; f_preview = NaN;
end

cancelled = true;
V_eff = NaN;
evalFunc = [];
C_eff_func = @(f) NaN;
circParamsOut = [];

fig = []; uit = []; eqCell = []; lblResult = []; btnApply = []; targetCapVar = '';

if isempty(V_preview) || isempty(f_preview) || isnan(V_preview) || isnan(f_preview)
    warning('RIM_CircuitBuilder:InvalidInput', 'Invalid or empty preview settings provided.');
    return;
end

V_eff = V_preview;
evalFunc = @(V, f) V;

defaultStr = 'Rs + p(Q_dl, Rct + W)';
if ~isempty(circParamsIn) && isfield(circParamsIn, 'circuitStr')
    defaultStr = circParamsIn.circuitStr;
end

fig = figure('Name', 'Equivalent Circuit Builder', 'Position', [250, 200, 700, 580], ...
    'MenuBar', 'none', 'NumberTitle', 'off', 'WindowStyle', 'modal');

statusStr = sprintf('Preview settings: V_app = %.3f V, f = %.2f Hz', V_preview, f_preview);
if ~isempty(circParamsIn), statusStr = [statusStr '  |  (Loaded Calibration Defaults)']; end

uicontrol(fig, 'Style', 'text', 'String', statusStr, ...
    'Units', 'normalized', 'Position', [0.05, 0.92, 0.90, 0.05], ...
    'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center');

uicontrol(fig, 'Style', 'text', 'String', 'Write your Equivalent Circuit (use p(A,B) for parallel):', ...
    'Units', 'normalized', 'Position', [0.05, 0.85, 0.60, 0.05], ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left');

eqCell = uicontrol(fig, 'Style', 'edit', 'String', defaultStr, ...
    'Units', 'normalized', 'Position', [0.05, 0.78, 0.60, 0.06], ...
    'HorizontalAlignment', 'left', 'FontSize', 12);

uicontrol(fig, 'Style', 'pushbutton', 'String', 'PARSE CIRCUIT', ...
    'Units', 'normalized', 'Position', [0.68, 0.78, 0.27, 0.06], ...
    'BackgroundColor', [0 0.4470 0.7410], 'ForegroundColor', 'w', 'FontWeight', 'bold', ...
    'Callback', @(~,~) parseCircuit());

uit = uitable(fig, 'Data', cell(0,2), 'ColumnName', {'Parameter', 'Value'}, ...
    'ColumnEditable', [false, true], 'Units', 'normalized', ...
    'Position', [0.05, 0.35, 0.90, 0.40]);

uicontrol(fig, 'Style', 'pushbutton', 'String', 'VIEW INSTRUCTIONS & LEGEND', ...
    'Units', 'normalized', 'Position', [0.05, 0.25, 0.90, 0.06], ...
    'Callback', @(~,~) showLegend());

lblResult = uicontrol(fig, 'Style', 'text', 'String', 'Calculated V_eff: (Parse & Preview first)', ...
    'Units', 'normalized', 'Position', [0.05, 0.15, 0.90, 0.05], ...
    'FontWeight', 'bold', 'ForegroundColor', [0 0 1], 'FontSize', 12, 'HorizontalAlignment', 'center');

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Preview Calculation', ...
    'Units', 'normalized', 'Position', [0.05, 0.05, 0.42, 0.07], ...
    'Callback', @(~,~) updatePreview());

btnApply = uicontrol(fig, 'Style', 'pushbutton', 'String', 'CONFIRM & APPLY', ...
    'Units', 'normalized', 'Position', [0.53, 0.05, 0.42, 0.07], ...
    'BackgroundColor', [0.4660 0.6740 0.1880], 'ForegroundColor', 'w', 'FontWeight', 'bold', ...
    'Enable', 'off', 'Callback', @(~,~) applyAndClose());

parseCircuit();
uiwait(fig);
if isvalid(fig), close(fig); end

    function parseCircuit()
        str = strtrim(eqCell.String);
        if isempty(str), errordlg('Circuit string cannot be empty.', 'Input Error'); return; end

        words = regexp(str, '[A-Za-z]\w*', 'match');
        words = unique(words);
        reserved = {'p', 's', 'V', 'w', 'beta', 'n'};
        parseVars = setdiff(words, reserved);

        isCap = startsWith(string(parseVars), "C", "IgnoreCase", true) | ...
            startsWith(string(parseVars), "Q", "IgnoreCase", true);
        capCandidates = parseVars(isCap);

        if isempty(capCandidates)
            errordlg('No Capacitor (C) or CPE (Q) found in your circuit string! Please add one.', 'Missing Component');
            btnApply.Enable = 'off'; return;
        end

        mappingCapVar = char(capCandidates{1});
        if length(capCandidates) > 1
            if ~isempty(circParamsIn) && isfield(circParamsIn, 'targetCapVar') && ismember(circParamsIn.targetCapVar, capCandidates)
                mappingCapVar = circParamsIn.targetCapVar;
            else
                [idx, tf] = listdlg('ListString', capCandidates, 'SelectionMode', 'single', ...
                    'PromptString', 'Multiple capacitors found. Which is the target interface?', ...
                    'Name', 'Select Interface', 'ListSize', [250 100]);
                if tf == 0, return; end
                mappingCapVar = char(capCandidates{idx});
            end
        end
        targetCapVar = mappingCapVar;

        hasCPE = any(cellfun(@(x) upper(x(1)) == 'Q', parseVars));
        if hasCPE && ~ismember('n', parseVars), parseVars{end+1} = 'n'; end

        newData = cell(length(parseVars), 2);
        for i = 1:length(parseVars)
            vName = parseVars{i}; newData{i, 1} = vName; hasSavedValue = false;

            if ~isempty(circParamsIn) && isfield(circParamsIn, 'paramsData')
                savedNames = string(circParamsIn.paramsData(:, 1));
                matchIdx = find(savedNames == string(vName), 1);
                if ~isempty(matchIdx)
                    newData{i, 2} = circParamsIn.paramsData{matchIdx, 2};
                    hasSavedValue = true;
                end
            end

            if ~hasSavedValue
                if strcmp(vName, 'n'), newData{i, 2} = 0.9;
                elseif strcmp(vName, targetCapVar), newData{i, 2} = 1e-5;
                elseif upper(vName(1)) == 'R', newData{i, 2} = 100;
                elseif upper(vName(1)) == 'W', newData{i, 2} = 15;
                elseif upper(vName(1)) == 'L', newData{i, 2} = 1e-6;
                else, newData{i, 2} = 1.0;
                end
            end
        end
        uit.Data = newData; btnApply.Enable = 'on';
        lblResult.String = sprintf('Target Interface: %s (Click Preview Calculation)', targetCapVar);
    end

    function showLegend()
        msg = {'--- CIRCUIT PARSER INSTRUCTIONS & COMPONENT LIBRARY ---', '', ...
            '1. TOPOLOGY & AUTO-DETECTION', 'Simply write your circuit. The engine auto-finds the interface C/Q.', '', ...
            '2. COMPONENT LIBRARY', 'R (Resistor), C (Capacitor), W (Warburg), L (Inductor), Q (CPE)', '', ...
            '3. SYNTAX RULES', 'Series: + (e.g. Rs + Rct)', 'Parallel: p(A,B) (e.g. p(Q_dl, Rct))'};
        msgbox(msg, 'Circuit Builder Instructions', 'help');
    end

    function updatePreview()
        if isempty(targetCapVar), errordlg('Please Parse the circuit first.', 'Error'); return; end
        [v, err] = calculateVoltageDrop(uit.Data, eqCell.String, targetCapVar, V_preview, f_preview);
        if ~isempty(err), errordlg(['Math Error: ' err], 'Error'); else, lblResult.String = sprintf('Calculated V_eff: %.4e V (Phase: %+.2f\x00B0)', abs(v), rad2deg(angle(v))); end
    end

    function applyAndClose()
        if isempty(targetCapVar), errordlg('Please Parse the circuit first.', 'Error'); return; end
        [v, err] = calculateVoltageDrop(uit.Data, eqCell.String, targetCapVar, V_preview, f_preview);
        if ~isempty(err), errordlg(['Fix errors before applying: ' err], 'Error'); return; end

        finalData = uit.Data; finalStr = eqCell.String; finalCap = targetCapVar;
        V_eff = v; evalFunc = @(V, f) circuitMathEngine(finalData, finalStr, finalCap, V, f);
        cancelled = false;

        circParamsOut.circuitStr = finalStr; circParamsOut.paramsData = finalData; circParamsOut.targetCapVar = finalCap;

        if upper(finalCap(1)) == 'Q'
            Q_val = NaN; n_val = 0.9;
            for i = 1:size(finalData, 1)
                if strcmp(finalData{i, 1}, finalCap), Q_val = finalData{i, 2}; end
                if strcmp(finalData{i, 1}, 'n'), n_val = finalData{i, 2}; end
            end
            if ischar(Q_val) || isstring(Q_val), Q_val = str2double(Q_val); end
            if ischar(n_val) || isstring(n_val), n_val = str2double(n_val); end
            C_eff_func = @(f) Q_val .* ((2 * pi * f).^(n_val - 1)) .* sin(n_val * pi / 2);
        else
            C_val = NaN;
            for i = 1:size(finalData, 1)
                if strcmp(finalData{i, 1}, finalCap), C_val = finalData{i, 2}; end
            end
            if ischar(C_val) || isstring(C_val), C_val = str2double(C_val); end
            C_eff_func = @(f) C_val;
        end
        uiresume(fig);
    end
end

function v_out = circuitMathEngine(paramsData, circuitStr, targetCap, V_in, f_in)
% Evaluates the effective voltage based on given circuit parameters.
[v_out, err] = calculateVoltageDrop(paramsData, circuitStr, targetCap, V_in, f_in);
if ~isempty(err), error('CircuitMathEngine: %s', err); end
end

function [v_out, err] = calculateVoltageDrop(paramsData, circuitStr, targetCap, V_in, f_in)
% Computes the theoretical voltage drop across the targeted capacitor.
err = ''; v_out = NaN;
try
    w = 2 * pi * f_in;
    p = @(a,b) (a .* b) ./ (a + b); %#ok<NASGU>

    valStruct = struct();
    for i = 1:size(paramsData, 1)
        varName = strtrim(char(paramsData{i, 1}));
        val = paramsData{i, 2};
        if ischar(val) || isstring(val), val = str2double(val); end
        if isnan(val), error('Invalid numeric value for parameter "%s".', varName); end
        valStruct.(varName) = val;
    end

    vars = fieldnames(valStruct);
    n_val = 1; if isfield(valStruct, 'n'), n_val = valStruct.n; end
    parsedStr = circuitStr;

    for i = 1:length(vars)
        var = vars{i};
        if strcmp(var, 'n'), continue; end
        val = valStruct.(var); firstLetter = upper(var(1));
        if firstLetter == 'R', Z_val = val;
        elseif firstLetter == 'C', Z_val = 1 / (1j * w * val);
        elseif firstLetter == 'W', Z_val = val * (1 - 1j) / sqrt(w);
        elseif firstLetter == 'L', Z_val = 1j * w * val;
        elseif firstLetter == 'Q', Z_val = 1 / (val * (1j * w)^n_val);
        else, Z_val = val;
        end

        eval(sprintf('Z_%s = %.16e + %.16ei;', var, real(Z_val), imag(Z_val)));
        searchPattern = sprintf('\\<%s\\>', var);
        parsedStr = regexprep(parsedStr, searchPattern, sprintf('Z_%s', var));
    end

    Z_tot = eval(parsedStr);
    eval(sprintf('Z_%s = 1e-20;', targetCap));
    Z_sh = eval(parsedStr);
    Z_interface = Z_tot - Z_sh;

    if isempty(Z_tot) || isnan(Z_tot) || abs(Z_tot) < 1e-6
        error('Total cell impedance evaluates to zero or NaN.');
    end
    v_out = V_in * (Z_interface / Z_tot);

catch ME
    err = ME.message;
end
end


%% ===================================================================
%  MATH MODULE FUNCTIONS
%  ===================================================================

function [max_amp, phase_at_max, freq_at_max] = math_findMaxInFreqInterval(timeSeries, frameInterval, freqStart, freqEnd)
% Performs a 1D FFT to find the peak amplitude and phase within a frequency interval.
max_amp = NaN; phase_at_max = NaN; freq_at_max = NaN;
if isempty(timeSeries) || frameInterval <= 0, return; end

N = length(timeSeries);
if N < 5
    return;
end

Fs = 1 / frameInterval;
freqStart = max(0, freqStart);
freqEnd = min(Fs/2, freqEnd);

try
    sig_proc = math_removeNonLinearDrift(timeSeries(:), 2) .* hann(N);
catch
    sig_proc = timeSeries(:) - mean(timeSeries);
end

Y = fft(sig_proc);
P2 = abs(Y / N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
P1 = P1 * 2;
f_axis = Fs * (0:(floor(N/2))) / N;

len = min(length(f_axis), length(P1));
f_axis = f_axis(1:len);
P1 = P1(1:len);

idx_start = find(f_axis >= freqStart, 1, 'first');
idx_end = find(f_axis <= freqEnd, 1, 'last');

if isempty(idx_start) || isempty(idx_end) || idx_start > idx_end, return; end

[max_amp, idx_rel] = max(P1(idx_start:idx_end));
idx_global = idx_start + idx_rel - 1;

freq_at_max = f_axis(idx_global);
phase_at_max = angle(Y(idx_global));
end

function [deltaR_map, phase_map] = math_calculateMapVectorized(img_stack, fr_int, f_start, f_end)
% Calculates amplitude and phase maps using a vectorized standard FFT.
[h, w, N] = size(img_stack);
Fs = 1 / fr_int;
f_axis = Fs * (0:(floor(N/2))) / N;

idx_start = find(f_axis >= f_start, 1, 'first');
idx_end = find(f_axis <= f_end, 1, 'last');

sig_2D = reshape(img_stack, h*w, N)';
sig_2D = math_removeNonLinearDrift(sig_2D, 2) .* hann(N);

Y = fft(sig_2D);
P2 = abs(Y / N);

P1 = P2(1:floor(N/2)+1, :);
P1(2:end-1, :);
P1(2:end-1, :) = 2 * P1(2:end-1, :);
P1 = P1 * 2;

P1_int = P1(idx_start:idx_end, :);
Y_int = Y(idx_start:idx_end, :);

[deltaR_flat, max_idx] = max(P1_int, [], 1);
linear_idx = sub2ind(size(Y_int), max_idx, 1:(h*w));
complex_peaks = Y_int(linear_idx);
phase_flat = angle(complex_peaks);

deltaR_map = reshape(deltaR_flat, h, w);
phase_map = reshape(phase_flat, h, w);
end

function [harmonicAmps, harmonicPhases] = math_extractHarmonicSignature(timeSeries, dt, f_fund, numHarmonics)
% Extracts the harmonic signature (amplitude and phase) from a reference signal.
N = length(timeSeries);
Fs = 1 / dt;

sig = math_removeNonLinearDrift(timeSeries(:), 2) .* hann(N);
Y = fft(sig);
P1 = Y(1:floor(N/2)+1) * 2 / N;
f_axis = Fs * (0:(floor(N/2))) / N;

harmonicAmps = ones(1, numHarmonics);
harmonicPhases = zeros(1, numHarmonics);

[~, idx_fund] = min(abs(f_axis - f_fund));
search_radius = max(2, floor(0.15 * f_fund / (Fs/N)));

[~, local_idx] = max(abs(P1(max(1, idx_fund-search_radius) : min(length(P1), idx_fund+search_radius))));
idx_fund_actual = max(1, idx_fund-search_radius) + local_idx - 1;

fund_complex = P1(idx_fund_actual);

for k = 1:numHarmonics
    f_target = k * f_fund;
    [~, idx_harm] = min(abs(f_axis - f_target));

    [~, local_idx] = max(abs(P1(max(1, idx_harm-search_radius) : min(length(P1), idx_harm+search_radius))));
    idx_harm_actual = max(1, idx_harm-search_radius) + local_idx - 1;

    harm_complex = P1(idx_harm_actual);

    harmonicAmps(k) = abs(harm_complex) / abs(fund_complex);
    harmonicPhases(k) = angle(harm_complex * exp(-1j * k * angle(fund_complex)));
end
end

function [K_I, K_Q] = math_buildSyntheticKernel(N, fr_int, f_opt, harmonicAmps, harmonicPhases)
% Synthesizes an analytical matched filter kernel based on extracted harmonics.
t = (0:N-1)' * fr_int;
K_I = zeros(N, 1);
K_Q = zeros(N, 1);

for k = 1:length(harmonicAmps)
    K_I = K_I + harmonicAmps(k) * cos(2*pi * k * f_opt * t + harmonicPhases(k));
    K_Q = K_Q + harmonicAmps(k) * sin(2*pi * k * f_opt * t + harmonicPhases(k));
end

K_I = K_I / sqrt(mean(K_I.^2) * 2);
K_Q = K_Q / sqrt(mean(K_Q.^2) * 2);
end

function [amp, phase] = math_calculateMatchedFilter1D(timeSeries, K_I, K_Q)
% Applies a 1D matched filter to a time series.
N = length(timeSeries);
sig = math_removeNonLinearDrift(timeSeries(:), 2);

X = (2/N) * dot(K_I(:), sig);
Y = (2/N) * dot(K_Q(:), sig);

amp = sqrt(X^2 + Y^2);
phase = atan2(Y, X);
end

function [deltaR_map, phase_map] = math_calculateMatchedFilterMap(img_stack, K_I, K_Q)
% Applies a vectorized 2D matched filter to an image stack.
[h, w, N] = size(img_stack);

sig_2D = reshape(img_stack, h*w, N)';
sig_2D = math_removeNonLinearDrift(sig_2D, 2);

X = (2/N) * (K_I(:)' * sig_2D);
Y = (2/N) * (K_Q(:)' * sig_2D);

deltaR_flat = sqrt(X.^2 + Y.^2);
phase_flat = atan2(Y, X);

deltaR_map = reshape(deltaR_flat, h, w);
phase_map = reshape(phase_flat, h, w);
end

function avgTimeSeries = math_calculateROIAverage(imageData, roiMask)
% Calculates the mean time-series for a selected spatial mask.
[h, w, numFrames] = size(imageData);
reshapedData = reshape(imageData, h*w, numFrames);
roiData = reshapedData(roiMask(:), :);
avgTimeSeries = mean(roiData, 1, 'omitnan')';
end

function imageStack = math_applySpatialFilter(imageStack, sigma)
% Applies a spatial Gaussian filter to the image stack.
if sigma <= 0, return; end
[~, ~, numFrames] = size(imageStack);
hF = waitbar(0, sprintf('Applying Gaussian filter (sigma=%.1f)...', sigma), 'Name', 'Filtering');

for i = 1:numFrames
    try
        imageStack(:,:,i) = imgaussfilt(imageStack(:,:,i), sigma);
        if mod(i, 50) == 0 && isvalid(hF), waitbar(i/numFrames, hF); end
    catch ME
        warning('RIM_Math:FilterError', 'Gaussian filter failed on frame %d: %s', i, ME.message);
        break;
    end
end

if isvalid(hF), close(hF); end
end

function [beta, q] = math_calculateBeta(V_list, deltaR_list, C_list)
% Calculates the calibration factor (beta) using linear regression.
Y = deltaR_list(:);
X = C_list(:) .* V_list(:);

valid = isfinite(X) & isfinite(Y);
if sum(valid) < 1, error('No valid points for fit.'); end

beta = X(valid) \ Y(valid);
q = 0;

if beta < 0
    warning('RIM_Math:NegativeBeta', 'Calibration returned negative beta (%.3e). Check V/deltaR sign convention.', beta);
end
end

function Z_mag = math_calculateImpedance(V, deltaR, beta, freq)
% Calculates specific impedance magnitude.
omega = 2 * pi * freq;
Z_mag = abs( (beta .* V ) ./ (deltaR .* omega) );
end

function C = math_calculateCapacitance(V, deltaR, beta)
% Calculates specific capacitance.
C = abs( deltaR ./ (beta .* V ) );
end

function [Z_real, Z_imag] = math_calculateComplexMaps(Z_mag, Z_phase)
% Converts polar impedance maps to Cartesian coordinates.
Z_real = Z_mag .* cos(Z_phase);
Z_imag = Z_mag .* sin(Z_phase);
end


%% ===================================================================
%  MEASUREMENT & ANALYSIS MODULE FUNCTIONS
%  ===================================================================

function updatedPath = RIM_Measurement(startPath)
% Executes the vectorized mapping workflow.
updatedPath = startPath;
calibFile = 'rim_calib_params.mat';

choice_mat = questdlg('How would you like to load the Calibration data?', ...
    'Calibration Source', ...
    'Load Portable .mat', 'Use Last Calibration', 'Load Portable .mat');

if string(choice_mat) == "Load Portable .mat"
    [mat_file, mat_path] = uigetfile('*.mat', 'Select Calibration MAT File', updatedPath);
    if isequal(mat_file, 0)
        disp('Measurement cancelled.'); return;
    end
    calibFilePath = fullfile(mat_path, mat_file);
elseif string(choice_mat) == "Use Last Calibration"
    calibFilePath = calibFile;
else
    disp('Measurement cancelled.'); return;
end

if ~exist(calibFilePath, 'file')
    errordlg('Calibration file not found. Please run Calibration first or select a valid .mat file.', 'Missing File');
    return;
end

try
    loaded = load(calibFilePath);
    disp(['>> Loaded Calibration File: ', calibFilePath]);
catch ME
    errordlg(['Failed to read MAT file: ', ME.message], 'File Error');
    return;
end

circParamsIn = [];
if isfield(loaded, 'circParamsToSave') && ~isempty(loaded.circParamsToSave)
    circParamsIn = loaded.circParamsToSave;
    disp('>> Restored Equivalent Circuit parameters.');
end

analysisChoice = questdlg('Select the Signal Extraction Method:', ...
    'Analysis Method', ...
    'Standard FFT (Linear Sine)', 'Faradaic Matched Filter (Non-Linear Q)', 'Standard FFT (Linear Sine)');

if isempty(analysisChoice)
    disp('Measurement cancelled.'); return;
end
useMatchedFilter = (string(analysisChoice) == "Faradaic Matched Filter (Non-Linear Q)");

[f_meas, fr_int, applyFilter, sigmaValue, V_meas, ~] = utils_promptAnalysisParams('Measurement');
if isempty(f_meas) || isnan(f_meas)
    disp('Measurement cancelled.'); return;
end

if useMatchedFilter
    choice_lasv = 'Yes';
else
    choice_lasv = questdlg('Load EC-Lab LASV file for electrical reference and phase calibration?', 'LASV Phase Sync', 'Yes', 'No', 'Yes');
end

has_lasv = false;
phase_shift = 0;
Z_ref_phase = NaN;
V_phase = 0;
I_phase = 0;

if string(choice_lasv) == "Yes"
    [lasv_file, lasv_path] = uigetfile({'*.mpt;*.txt;*.csv', 'EC-Lab Data Files (*.mpt, *.txt, *.csv)'}, 'Select LASV File', updatedPath);
    if lasv_file ~= 0
        try
            [t_el, V_el, I_el] = utils_loadECLabLASV(fullfile(lasv_path, lasv_file));
            dt_el = mean(diff(t_el));

            V_win = V_el;
            e_idx1 = 1;
            e_idx2 = length(V_el);

            search_radius = max(0.5, 0.2 * f_meas);
            [V_amp, V_phase, f_peak_el] = math_findMaxInFreqInterval(V_win, dt_el, max(0.1, f_meas - search_radius), f_meas + search_radius);

            if isnan(V_amp)
                errordlg('Could not auto-detect a valid voltage peak in LASV data.', 'Analysis Error');
                return;
            end

            I_win = I_el(e_idx1:e_idx2);
            [I_amp, I_phase, ~] = math_findMaxInFreqInterval(I_win, dt_el, f_peak_el - 0.05, f_peak_el + 0.05);

            Z_ref_mag = V_amp / I_amp;
            Z_ref_phase = V_phase - I_phase;

            V_meas = V_amp;
            has_lasv = true;

            fprintf('\n--- EC-Lab LASV Reference ---\n');
            fprintf('Frequency: %.3f Hz\n', f_peak_el);
            fprintf('Measured V_amp: %.4e V\n', V_amp);
            fprintf('Measured I_amp: %.4e A\n', I_amp);
            fprintf('Reference |Z|: %.2f Ohm\n', Z_ref_mag);
            fprintf('Reference Phase: %.2f deg (%.2f rad)\n', rad2deg(Z_ref_phase), Z_ref_phase);

        catch ME
            errordlg(sprintf('Failed to read or process LASV file:\n%s\n\nFalling back to manual voltage.', ME.message), 'LASV Error');
            if useMatchedFilter
                warndlg('Cannot proceed with Matched Filter without LASV. Reverting to Standard FFT.', 'Method Changed');
                useMatchedFilter = false;
            end
        end
    else
        if useMatchedFilter
            warndlg('LASV file selection cancelled. Reverting to Standard FFT.', 'Method Changed');
            useMatchedFilter = false;
        end
    end
end

[img, ~, h, w, updatedPath, true_dt] = utils_loadTiffFolder(updatedPath, 'Measurement Data', applyFilter, sigmaValue);
if isempty(img), return; end
opticalImg = img(:,:,1);

if ~isnan(true_dt) && true_dt > 0
        % Calculate the difference between input and reality
        time_diff_pct = abs(true_dt - fr_int) / fr_int * 100;
        
        % If the difference is meaningful (e.g., greater than 0.1%)
        if time_diff_pct > 0.1 
            msg = sprintf('Timing Mismatch Detected!\n\nManual Input: %.5f s (%.1f fps)\nCamera Hardware: %.5f s (%.1f fps)\n\nWhich interval should the script use for phase extraction?', fr_int, 1/fr_int, true_dt, 1/true_dt);
            
            choice_time = questdlg(msg, 'Phase Synchronization Warning', ...
                'Use Hardware (Recommended)', 'Keep Manual Input', 'Use Hardware (Recommended)');
                
            if string(choice_time) == "Use Hardware (Recommended)"
                fprintf('>> Opted to use hardware true dt: %.5f s\n', true_dt);
                fr_int = true_dt; 
            else
                fprintf('>> Warning: Kept manual input dt: %.5f s. Phase maps may drift.\n', fr_int);
            end
        else
            % If they basically match, silently adopt the hardware precision
            fr_int = true_dt; 
            fprintf('>> Hardware timing matches input. Using precision dt: %.6f s\n', true_dt);
        end
    elseif isnan(true_dt)
        warning('Could not extract hardware timestamps. Proceeding with manual input.');
end

if isfield(loaded, 'betaValue')
    betaValue = loaded.betaValue;
    disp(['>> Applied Single Bulk Beta: ', num2str(betaValue)]);
else
    errordlg('Calibration file does not contain a Single Beta value. Run Calibration again.', 'Calibration Error');
    return;
end

if isfield(loaded, 'qValue')
    qValue = loaded.qValue;
    disp(['>> Applied Intercept (q): ', num2str(qValue)]);
else
    qValue = 0;
end

BetaMap = ones(h, w) * betaValue;

[roiMask, ~] = utils_selectROI(opticalImg, 'Select Setup ROI (High Signal)');
if isempty(roiMask), return; end

ts_full = math_calculateROIAverage(img, roiMask);
[ts_win, idx_start, idx_end] = utils_selectWindowFromSpectrogram(ts_full, fr_int, f_meas, 'Measurement Window Selection');
if isempty(ts_win), return; end

img_windowed = img(:, :, idx_start:idx_end);

[~, theta_opt_roi, f_peak_opt, canc] = utils_selectFFTPeak(ts_win, fr_int, f_meas, 'Optical Sync Check');
if canc || isempty(f_peak_opt) || isnan(f_peak_opt), return; end

if useMatchedFilter && has_lasv
    disp('--- Building Synthesized Faradaic Kernel ---');
    Q_el = cumtrapz(t_el, I_el);
    Q_win = Q_el(e_idx1:e_idx2);

    numHarmonics = 5;
    [h_amps, h_phases] = math_extractHarmonicSignature(Q_win, dt_el, f_peak_el, numHarmonics);

    fprintf('Extracted Harmonics Signature (Relative Amplitudes):\n');
    disp(h_amps);

    numFrames = size(img_windowed, 3);
    [K_I, K_Q] = math_buildSyntheticKernel(numFrames, fr_int, f_peak_opt, h_amps, h_phases);

    t_opt_win = (0:numFrames-1) * fr_int;
    figK = figure('Name', 'Synthetic Matched Filter Kernels', 'Position', [200, 200, 800, 400]);
    plot(t_opt_win, K_I, 'b', 'LineWidth', 1.5); hold on;
    plot(t_opt_win, K_Q, 'r', 'LineWidth', 1.5);
    title(sprintf('Synthetic Faradaic Kernels (Synchronized to %.3f Hz Camera Clock)', f_peak_opt));
    legend('In-Phase (Synthesized Charge Q)', 'Quadrature (90^\circ Shift)', 'Location', 'best');
    xlabel('Time (s)'); ylabel('Normalized Amplitude'); grid on;

    [~, theta_opt_roi] = math_calculateMatchedFilter1D(ts_win, K_I, K_Q);
end

if has_lasv
    target_opt_phase = I_phase + (pi / 2);
    phase_shift = target_opt_phase - theta_opt_roi;
    fprintf('Hardware Sync: Shifting optical phase maps by %.2f deg.\n', rad2deg(phase_shift));
else
    Z_ref_phase = angle(exp(1j * V_phase - (theta_opt_roi + pi/2)));
end

choice = questdlg('Apply Equivalent Circuit Voltage Drop Correction?', 'Circuit Builder', 'Yes', 'No', 'No');
if string(choice) == "Yes"
    [V_eff_complex, ~, canc_circ, ~, ~] = RIM_CircuitBuilder(V_meas, f_peak_opt, circParamsIn);
    if canc_circ
        disp('Measurement cancelled during Circuit Builder.');
        return;
    end
    fprintf('\n--- Circuit Correction Applied ---\n');
    fprintf('V_applied = %.3f V --> V_effective = %.4e V\n', V_meas, abs(V_eff_complex));
    fprintf('Phase Shift Introduced: %+.2f deg\n', rad2deg(angle(V_eff_complex)));
    V_meas = abs(V_eff_complex);
    V_phase = V_phase + angle(V_eff_complex);
end

fprintf('Calculating Map at %.3f Hz (Window indices: %d to %d)...\n', f_peak_opt, idx_start, idx_end);

if useMatchedFilter
    fprintf('>>> Extracting signal using Custom Non-Linear Lock-In (Synthetic Matched Filter)...\n');
    [dR_map_raw, theta_map] = math_calculateMatchedFilterMap(img_windowed, K_I, K_Q);
else
    fprintf('>>> Extracting signal using Standard FFT...\n');
    f_start = f_peak_opt - max(0.2, 0.1*f_peak_opt);
    f_end = f_peak_opt + max(0.2, 0.1*f_peak_opt);
    [dR_map_raw, theta_map] = math_calculateMapVectorized(img_windowed, fr_int, f_start, f_end);
end
mean_img = mean(img_windowed, 3);

mean_img(mean_img == 0) = 1;
dR_map = dR_map_raw ./ mean_img;

theta_map = theta_map + phase_shift;

ZmagMap = NaN(h, w);
ZphaseMap = ones(h, w) * Z_ref_phase;

effective_dR_map = dR_map - qValue;
invalid_pixels = effective_dR_map <= 1e-6;
validMask = ~isnan(dR_map) & (dR_map > 0) & (opticalImg > 0) & ~invalid_pixels;

ZmagMap(validMask) = math_calculateImpedance(V_meas, effective_dR_map(validMask), BetaMap(validMask), f_peak_opt);
ZphaseMap(validMask) = angle(exp(1j * (V_phase - theta_map(validMask) + pi/2)));

tempPhase = ZphaseMap(validMask);
tempPhase = atan2(sin(tempPhase), cos(tempPhase));
ZphaseMap(validMask) = tempPhase;

disp('Mapping complete.');

ZmagPlot = ZmagMap;
ZphasePlot = ZphaseMap;
titlePrefix = '';

figure('Name', 'RIM Spatial Maps', 'Position', [100, 100, 1600, 450]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
imagesc(ax1, opticalImg); axis(ax1, 'tight'); colormap(ax1, 'gray');
title(ax1, 'Optical Frame');

ax2 = nexttile;
hMag = imagesc(ax2, ZmagPlot); axis(ax2, 'tight');
set(hMag, 'AlphaData', ~isnan(ZmagPlot));
colormap(ax2, 'turbo'); colorbar(ax2);
set(ax2, 'Color', [0.9 0.9 0.9]);
title(ax2, [titlePrefix '|Z| (\Omega \cdot cm^2)']);

validMaskPlot = ~isnan(ZmagPlot) & (ZmagPlot > 0);
if any(validMaskPlot(:))
    valid_Z = ZmagPlot(validMaskPlot);
    valid_Z = valid_Z(isfinite(valid_Z) & valid_Z > 0);
    if ~isempty(valid_Z)
        med_Z = median(valid_Z);
        c_max = prctile(valid_Z, 95);
        if c_max > (med_Z * 10), c_max = med_Z * 10; end
        if c_max <= 0 || isnan(c_max), c_max = 1; end
        try clim(ax2, [0, c_max]); catch; end
    end
end

ax3 = nexttile;
imagesc(ax3, ZphasePlot); axis(ax3, 'tight');
colormap(ax3, 'turbo'); colorbar(ax3);
set(ax3, 'Color', [0.9 0.9 0.9]);
title(ax3, [titlePrefix 'Phase (rad)']);

if any(validMaskPlot(:))
    medPhase = median(ZphasePlot(validMaskPlot), 'omitnan');
    stdPhase = std(ZphasePlot(validMaskPlot), 'omitnan');
    cRange = min(2 * stdPhase, 0.25);
    if cRange < 0.01, cRange = 0.01; end
    try clim(ax3, [medPhase - cRange, medPhase + cRange]); catch; end
end

choice = questdlg('Export Impedance, Phase, and Optical Maps?', 'Export Maps', 'Yes', 'No', 'No');
if string(choice) == "Yes"
    defName = fullfile(updatedPath, 'ImpedanceMap');
    [fName, pName] = uiputfile('*.*', 'Save Output (Base Name)', defName);
    if ~isequal(fName, 0)
        [~, baseName, ~] = fileparts(fName);
        try
            writematrix(ZmagMap, fullfile(pName, [baseName '_Magnitude.csv']));
            writematrix(ZphaseMap, fullfile(pName, [baseName '_Phase.csv']));
            writematrix(opticalImg, fullfile(pName, [baseName '_Optical.csv']));

            ZmagExport = ZmagMap;
            ZmagExport(isinf(ZmagExport)) = NaN;
            utils_saveFloatTiff(ZmagExport, fullfile(pName, [baseName '_Magnitude.tif']));
            utils_saveFloatTiff(ZphaseMap, fullfile(pName, [baseName '_Phase.tif']));
            utils_saveFloatTiff(opticalImg, fullfile(pName, [baseName '_Optical.tif']));

            if exist('figK', 'var') && isgraphics(figK)
                try exportgraphics(figK, fullfile(pName, [baseName '_Synthetic_Kernels.png']), 'Resolution', 300); catch; end
            end

            msgbox('Maps and figures exported successfully.', 'Export Complete');
        catch ME
            errordlg(sprintf('An error occurred while saving the files:\n%s', ME.message), 'Export Error');
        end
    end
end
end

function updatedPath = RIM_RedoxAnalysis(startPath)
% Processes and visualizes the real and imaginary components of the impedance map.
updatedPath = startPath;
warning('off', 'MATLAB:imagesci:Tiff:libraryWarning');

[fNameM, pNameM] = uigetfile({'*.csv;*.tif;*.tiff', 'Map Files (*.csv, *.tif, *.tiff)'}, 'Select Z Magnitude Map', updatedPath);
if isequal(fNameM,0), return; end
[fNameP, pNameP] = uigetfile({'*.csv;*.tif;*.tiff', 'Map Files (*.csv, *.tif, *.tiff)'}, 'Select Z Phase Map', pNameM);
if isequal(fNameP,0), return; end

try
    Z_mag = loadMapData(fullfile(pNameM, fNameM));
    Z_phase = loadMapData(fullfile(pNameP, fNameP));
catch ME
    errordlg(['Error loading maps: ', ME.message], 'File Read Error');
    return;
end

[Z_real, Z_imag] = math_calculateComplexMaps(Z_mag, Z_phase);
Z_imag_plot = -Z_imag;

figName = "Impedance Components (Real & Imaginary)";
figPos = [100 100 1000 450];
layoutCols = 2;

figure("Name", figName, "Position", figPos);
tiledlayout(1, layoutCols, "Padding", "compact", "TileSpacing", "compact");

ax1 = nexttile;
imagesc(ax1, Z_real); axis(ax1, "tight");
colorbar(ax1); colormap(ax1, "turbo"); title(ax1, "Real Part (Z') [\Omega \cdot cm^2]");
applySafeClim(ax1, Z_real);

ax2 = nexttile;
imagesc(ax2, Z_imag_plot); axis(ax2, "tight");
colorbar(ax2); colormap(ax2, "turbo"); title(ax2, "Imaginary Part (-Z'') [\Omega \cdot cm^2]");
applySafeClim(ax2, Z_imag_plot);

choice_export = questdlg("Export the generated maps as TIFF and CSV?", "Export Maps", "Yes", "No", "Yes");
if string(choice_export) == "Yes"
    [fName, pName] = uiputfile("*.csv", "Save Base Name", fullfile(updatedPath, "Real_Imag_Maps"));
    if ~isequal(fName, 0)
        [~, baseName, ~] = fileparts(fName);

        utils_saveFloatTiff(Z_real, fullfile(pName, [baseName '_Real.tif']));
        utils_saveFloatTiff(Z_imag_plot, fullfile(pName, [baseName '_Nyquist_Imag.tif']));
        writematrix(Z_real, fullfile(pName, [baseName '_Real.csv']));
        writematrix(Z_imag_plot, fullfile(pName, [baseName '_Nyquist_Imag.csv']));

        msgbox("Export complete. 2 Maps saved as both 32-bit TIFF and CSV.", "Success");
    end
end
end

function data = loadMapData(filePath)
% Reads map data from either CSV or TIFF files.
[~, ~, ext] = fileparts(filePath);

if strcmpi(ext, '.csv')
    data = readmatrix(filePath);
elseif strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
    t = Tiff(filePath, 'r');
    data = double(t.read());
    close(t);
else
    error('Unsupported file format. Please use .csv or .tif/.tiff.');
end
end

function applySafeClim(ax, data)
% Applies a bounded color limit to handle outliers in map visualization.
validData = data(isfinite(data));
if isempty(validData)
    return;
end

c_min = prctile(validData, 2);
c_max = prctile(validData, 95);

if isnan(c_max) || isnan(c_min) || c_max <= c_min
    c_max = c_min + 1e-6;
end

try clim(ax, [c_min, c_max]); catch; end
end

function [detrendedData] = math_removeNonLinearDrift(data, polyOrder)

if nargin < 2
    polyOrder = 2; % Default to a 2nd-order (parabolic) curve
end
N = size(data, 1);
t = linspace(-1, 1, N)'; 
V = zeros(N, polyOrder + 1);
for i = 0:polyOrder
    V(:, i+1) = t.^(polyOrder - i);
end
trend = V * (V \ data); 
detrendedData = data - trend;
end


%% ===================================================================
%  UTILITIES MODULE FUNCTIONS
%  ===================================================================

function [f_target, fr_int, applyFilter, sigmaValue, dV, area] = utils_promptAnalysisParams(modeName)
% Opens a dialog box to collect mode-specific numerical parameters from the user.
f_target = []; fr_int = []; applyFilter = false; sigmaValue = 0; dV = []; area = 1.0;

if strcmp(modeName, 'Calibration')
    prompt = {
        'Target Frequency (Hz):', ...
        'Frame Interval (s) [e.g., 0.01 for 100fps]:', ...
        'Gaussian Filter Sigma (0 to disable):', ...
        'Electrode Active Area (cm^2):'
        };
    defaults = {'5.0', '0.01', '1.0', '1.0'};

elseif strcmp(modeName, 'Measurement')
    prompt = {
        'Target Frequency (Hz):', ...
        'Frame Interval (s) [e.g., 0.01 for 100fps]:', ...
        'Electrode Active Area (cm^2):', ...
        'Oscillation Amplitude Delta V (V, 0-peak):'
        };
    defaults = {'5.0', '0.01', '1.0', '0.4'};
else
    prompt = {'Target Frequency (Hz):', 'Frame Interval (s):', 'Electrode Active Area (cm^2):'};
    defaults = {'5.0', '0.01', '1.0'};
end

titleBar = [modeName ' Settings'];
answer = inputdlg(prompt, titleBar, [1 60], defaults);

if isempty(answer)
    return;
end

f_target = str2double(answer{1});
fr_int = str2double(answer{2});

if strcmp(modeName, 'Calibration')
    sigmaValue = str2double(answer{3});
    area = str2double(answer{4});
elseif strcmp(modeName, 'Measurement')
    sigmaValue = 0;
    area = str2double(answer{3});
    dV = str2double(answer{4});
else
    sigmaValue = 0;
    area = str2double(answer{3});
end

if isnan(f_target) || isnan(fr_int) || isnan(sigmaValue) || isnan(area) || (strcmp(modeName, 'Measurement') && isnan(dV))
    errordlg('Invalid numeric values entered. Please check your inputs.', 'Input Error');
    f_target = [];
    return;
end

applyFilter = (sigmaValue > 0);
if applyFilter && (~license('test', 'image_toolbox') || isempty(ver('images')))
    warndlg('Image Processing Toolbox not found. Filtering disabled.');
    applyFilter = false;
end
end

function [img, numFrames, h, w, updatedPath, true_dt] = utils_loadTiffFolder(currPath, label, applyFilter, sigmaValue)
% Reads a sequence of TIFF images from a selected folder into a 3D matrix.
img=[]; numFrames=0; h=0; w=0; updatedPath=currPath;

pn = uigetdir(currPath, ['Select Folder for: ' label]);
if isequal(pn,0), return; end

files = dir(fullfile(pn, '*.tif'));
if isempty(files), files = dir(fullfile(pn, '*.tiff')); end
if isempty(files), warning('No TIF files found in selected folder.'); return; end

fileNames = utils_naturalSort({files.name});
numFrames = length(fileNames);
info = imfinfo(fullfile(pn, fileNames{1}));
w = info.Width; h = info.Height;
img = zeros(h,w,numFrames, 'single');
timestamps = zeros(numFrames, 1);

hW = waitbar(0, sprintf('Loading %d TIFFs...', numFrames));
try
    for k=1:numFrames
        if ~isvalid(hW), error('User Cancelled Loading'); end

        t = Tiff(fullfile(pn, fileNames{k}), 'r');
        img(:,:,k) = single(t.read());
        close(t);

        % --- NEW: Extract Hardware Timestamp ---
        try
            info = imfinfo(fullfile(pn, fileNames{k}));
            timeToken = regexp(info.ImageDescription, 'bofTime=([0-9.]+)us', 'tokens');
            if ~isempty(timeToken)
                timestamps(k) = str2double(timeToken{1}{1}) / 1000000; % Convert to seconds
            else
                timestamps(k) = NaN;
            end
        catch
            timestamps(k) = NaN;
        end
        % ---------------------------------------

        if mod(k, 50) == 0 && isvalid(hW), waitbar(k/numFrames, hW); end
    end

    % Calculate true delta t
    dt_array = diff(timestamps);
    true_dt = mean(dt_array, 'omitnan');
catch ME
    if isvalid(hW), close(hW); end
    warning('RIM_Utils:LoadError', 'Error loading TIFFs: %s', ME.message);
    img = []; return;
end
if isvalid(hW), close(hW); end
updatedPath = pn;

if nargin > 2 && applyFilter
    img = math_applySpatialFilter(img, sigmaValue);
end
end

function [mask, pos] = utils_selectROI(img, titleStr)
% Provides an interactive tool to define a rectangular Region of Interest.
mask=[]; pos=[];
f = figure('Name', titleStr); imshow(img, []); title('Draw Rect, Double Click to Confirm');
try
    h = drawrectangle('Label', 'ROI');
    wait(h);
    if isvalid(h)
        mask = createMask(h);
        pos = h.Position;
    end
    if isvalid(f), close(f); end
catch ME
    warning('RIM_Utils:ROIError', 'ROI selection interrupted: %s', ME.message);
    if isvalid(f), close(f); end
end
end

function [ts_out, idx1, idx2] = utils_selectWindowFromSpectrogram(ts_in, fr_int, f_target, titleStr)
% Provides an interactive spectrogram tool to bound the temporal analysis window.
ts_out = []; idx1 = 1; idx2 = length(ts_in);
fs = 1 / fr_int;
N = length(ts_in);
t_axis = (0:N-1) * fr_int;

windowLength = max(4, floor(N / 8));
noverlap = floor(windowLength * 0.9);
nfft = max(256, 2^nextpow2(windowLength));

f = figure('Name', ['Spectrogram Selection - ' titleStr], 'Position', [100, 100, 900, 600]);
try
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile;
    plot(t_axis, ts_in);
    title('Raw Time Series'); xlabel('Time (s)'); ylabel('Amplitude');
    grid on;

    ax2 = nexttile;
    [~, F_spec, T_spec, P_spec] = spectrogram(math_removeNonLinearDrift(ts_in(:), 2), windowLength, noverlap, nfft, fs);

    imagesc(ax2, T_spec, F_spec, 10*log10(abs(P_spec)));
    axis xy;

    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram (Target Freq: %.2f Hz)', f_target));

    if ~isempty(f_target) && ~isnan(f_target)
        ylim(ax2, [0, min(fs/2, f_target * 2.5)]);
    end
    colormap(ax2, 'turbo');

    title(tl, 'Draw a rectangle over the steady-state region, then double-click to confirm', 'FontSize', 12, 'FontWeight', 'bold');

    roi = drawrectangle(ax2, 'Color', [0 0.4470 0.7410], 'Label', 'Window');
    wait(roi);

    if isvalid(roi)
        t_start = roi.Position(1);
        t_end = t_start + roi.Position(3);
        delete(roi);

        xregion(ax1, t_start, t_end, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.2);
        xregion(ax2, t_start, t_end, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.2);
        pause(0.5);

        idx1 = max(1, round(t_start/fr_int) + 1);
        idx2 = min(N, round(t_end/fr_int) + 1);
        ts_out = ts_in(idx1:idx2);
    end
catch ME
    warning('RIM_Utils:SpectrogramError', 'Spectrogram window selection interrupted: %s', ME.message);
end
if isvalid(f), close(f); end
end

function [amp, phase, freq, cancelled] = utils_selectFFTPeak(ts, fi, f_target, ~)
% Interactive tool to isolate a targeted frequency peak from an FFT spectrum.
amp=NaN; phase=NaN; freq=NaN; cancelled=true;
if isempty(ts), return; end

N = length(ts);
if N < 5, return; end

Y = fft(math_removeNonLinearDrift(ts(:), 2) .* hann(N));
P1 = abs(Y/N); P1 = P1(1:floor(N/2)+1); P1(2:end-1)=2*P1(2:end-1);
f = (0:(floor(N/2)))/N * (1/fi);

hF = figure('Name', 'FFT Peak Selection');
try
    tiledlayout(1, 1, 'Padding', 'compact');
    ax = nexttile;

    plot(ax, f, P1, 'LineWidth', 1.5);
    if isempty(f_target) || isnan(f_target)
        xlim(ax, [0 max(f)]);
    else
        xlim(ax, [0 max(f_target*2, 10)]);
    end
    title(ax, 'Draw a rectangle around the central peak, then double-click to confirm.', 'FontWeight', 'bold');
    xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Amplitude');

    roi = drawrectangle(ax, 'Color', [0.8500 0.3250 0.0980]);
    wait(roi);

    if isvalid(roi)
        f_start = roi.Position(1);
        f_end = f_start + roi.Position(3);
        delete(roi);

        xregion(ax, f_start, f_end, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2);
        pause(0.5);

        [amp, phase, freq] = math_findMaxInFreqInterval(ts, fi, f_start, f_end);
        cancelled = false;
    end
catch ME
    warning('RIM_Utils:FFTPeakError', 'FFT Peak selection interrupted: %s', ME.message);
end
if isvalid(hF), close(hF); end
end

function [V_out, dR_out, mask] = utils_excludeOutliers(V, dR)
% Interactive lasso tool for removing anomalous data points from a calibration dataset.
f = figure('Name', 'Outlier Removal', 'Position', [150 150 800 500]);

uicontrol('Style', 'pushbutton', 'String', 'CONFIRM & CONTINUE', ...
    'Units', 'normalized', 'Position', [0.35 0.02 0.3 0.08], ...
    'FontWeight', 'bold', 'BackgroundColor', '#77AC30', 'ForegroundColor', 'w', ...
    'Callback', @(~, ~) utils_finalizeOutlierSelection(f));

setappdata(f, 'done', false);

tl = tiledlayout(1, 1, 'Padding', 'compact');
tl.OuterPosition = [0 0.1 1 0.9];
ax = nexttile;

hValid = plot(ax, V, dR, 'o', 'MarkerSize', 8, 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Valid Data');
hold(ax, 'on');
hOut = plot(ax, NaN, NaN, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Excluded');
hold(ax, 'off');
legend(ax, 'show', 'Location', 'best');
grid(ax, 'on');

title(ax, 'LASSO TOOL: Freehand draw around outliers to remove them.', 'FontWeight', 'bold');
xlabel(ax, 'V'); ylabel(ax, 'dR');
keep = true(size(V));

while isvalid(f) && ~getappdata(f, 'done')
    try
        roi = drawfreehand(ax, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'LabelVisible', 'off');
        wait(roi);

        if ~isvalid(roi) || getappdata(f, 'done')
            break;
        end

        in = inpolygon(V, dR, roi.Position(:,1), roi.Position(:,2));
        keep(in) = false;

        hValid.XData = V(keep); hValid.YData = dR(keep);
        hOut.XData = V(~keep); hOut.YData = dR(~keep);

        delete(roi);
    catch
        break;
    end
end

if isvalid(f), close(f); end
V_out = V(keep); dR_out = dR(keep); mask = keep;
end

function utils_finalizeOutlierSelection(fig)
% Helper function to trigger the closure of the outlier selection loop.
if isvalid(fig)
    setappdata(fig, 'done', true);
    delete(findobj(fig, 'Type', 'images.roi.Freehand'));
    uiresume(fig);
end
end

function sortedList = utils_naturalSort(listIn)
% Sorts alphanumeric strings mathematically (e.g., handles '2' before '10').
if isstring(listIn), listIn = cellstr(listIn); end

expressions = regexp(listIn, '\d+', 'match');
if isempty(expressions) || all(cellfun(@isempty, expressions))
    sortedList = sort(listIn);
    return;
end

paddedList = regexprep(listIn, '\d+', '${sprintf(''%010d'', str2double($0))}');
[~, idx] = sort(paddedList);
sortedList = listIn(idx);
end

function [t, V, I] = utils_loadECLabLASV(filePath)
% Parses EC-Lab exported text files to extract time, voltage, and current arrays.
fid = fopen(filePath, 'rt');
if fid == -1, error('Could not open the selected LASV file.'); end

rawText = fileread(filePath);
rawText = strrep(rawText, ',', '.');

lines = regexp(rawText, '\r?\n', 'split');
headerIdx = 1;
for i = 1:min(500, length(lines))
    str = lower(lines{i});
    if contains(str, 'time') && (contains(str, 'ewe') || contains(str, 'e_v')) && contains(str, 'i')
        headerIdx = i; break;
    end
end

tempFile = [tempname, '.csv'];
fid = fopen(tempFile, 'wt');
fprintf(fid, '%s', rawText);
fclose(fid);

warnState1 = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
warnState2 = warning('off', 'MATLAB:table:ModifiedVarnames');
try
    opts = detectImportOptions(tempFile, 'NumHeaderLines', headerIdx-1);
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(tempFile, opts);
catch
    try
        tbl = readtable(tempFile, 'HeaderLines', headerIdx-1, 'VariableNamingRule', 'preserve');
    catch
        tbl = readtable(tempFile, 'HeaderLines', headerIdx-1);
    end
end
warning(warnState1);
warning(warnState2);
delete(tempFile);

colNames = tbl.Properties.VariableNames;
t_idx = find(contains(colNames, 'time', 'IgnoreCase', true), 1);
v_idx = find(contains(colNames, 'Ewe', 'IgnoreCase', true) | ...
    startsWith(colNames, 'E_V', 'IgnoreCase', true) | ...
    startsWith(colNames, 'E_', 'IgnoreCase', true) | ...
    startsWith(colNames, 'E/', 'IgnoreCase', true), 1);
i_idx = find(contains(colNames, 'I_', 'IgnoreCase', true) | ...
    contains(colNames, 'I/', 'IgnoreCase', true) | ...
    strcmp(colNames, 'I'), 1);

if isempty(t_idx) || isempty(v_idx) || isempty(i_idx)
    error('Could not auto-detect Time, Voltage, or Current columns in the file headers.');
end

t = tbl{:, t_idx}; V = tbl{:, v_idx}; I = tbl{:, i_idx};
if contains(colNames{i_idx}, 'mA', 'IgnoreCase', true), I = I * 1e-3; end

valid = isfinite(t) & isfinite(V) & isfinite(I);
t = t(valid); V = V(valid); I = I(valid);

if numel(t) < 10
    error('LASV file contains too few valid rows (%d). File may be corrupted or format not recognized.', numel(t));
end
end

function utils_saveFloatTiff(imgData, filePath)
% Saves a numeric 2D array as a standard 32-bit floating-point TIFF.
t = Tiff(filePath, 'w');
tags.ImageLength = size(imgData, 1);
tags.ImageWidth = size(imgData, 2);
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.BitsPerSample = 32;
tags.SamplesPerPixel = 1;
tags.SampleFormat = Tiff.SampleFormat.IEEEFP;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

t.setTag(tags);
t.write(single(imgData));
t.close();
end