%% Imaging re-analysis
% 250903 updated (v1.5)
% You can apply an average intensity threshold to filter out
% low-mean-intensity clusters in each channel.

clc; clear; close('all');

%% Set path & file search (align with App Designer naming)
addpath(cd)
fold_path = uigetdir;                % Select folder via UI
cd(fold_path)                        % Change working directory
warning('off')

% Target .mat files (e.g., saved from the app: result-*.mat) OR raw inputs
mat_info = dir('*.mat');             % Find MAT files
mat_names = {mat_info.name};         % Collect file names
fname_list = mat_names;

%% Parameter input (match app variable names/style)
background_setting        = 1;   % green channel (default: 1)
peak_background_fold      = 2;   % green channel (default: 2)
tolerance_parameter       = 1.5; % green channel (default: 1.5)

redch_background_setting  = 1;   % red channel (default: 1)
redch_peak_background_fold= 2;   % red channel (default: 2)
redch_tolerance_parameter = 1.5; % red channel (default: 1.5)

% Average intensity thresholds (new feature; follow app-like naming)
avg_int_threshold         = 50; % 0 = no threshold
redch_avg_int_threshold   = 50; % 0 = no threshold

%% Main loop
for kk = 1:length(fname_list)
    tmp_name = fname_list{kk};
    tmp_dat  = load(tmp_name);

    % Expect structure saved by the app (img_ch_* present)
    % If running on raw “dat.ch1/ch2”, adapt here accordingly.
    if isfield(tmp_dat, 'img_ch_green_data') && isfield(tmp_dat, 'img_ch_red_data')
        % Already processed image data exists — reuse images if present
        % (Fallback to fields inside if necessary)
        if isfield(tmp_dat.img_ch_green_data, 'cell_selc')
            img_ch_green = tmp_dat.img_ch_green_data.cell_selc;
        elseif isfield(tmp_dat, 'dat') && isfield(tmp_dat.dat, 'ch1')
            img_ch_green = tmp_dat.dat.ch1.Image;
        else
            error('Green channel image not found in %s', tmp_name);
        end

        if isfield(tmp_dat.img_ch_red_data, 'cell_selc')
            img_ch_red = tmp_dat.img_ch_red_data.cell_selc;
        elseif isfield(tmp_dat, 'dat') && isfield(tmp_dat.dat, 'ch2')
            img_ch_red = tmp_dat.dat.ch2.Image;
        else
            error('Red channel image not found in %s', tmp_name);
        end
    elseif isfield(tmp_dat, 'dat') && isfield(tmp_dat.dat, 'ch1') && isfield(tmp_dat.dat, 'ch2')
        % Raw data structure
        img_ch_green = tmp_dat.dat.ch1.Image;
        img_ch_red   = tmp_dat.dat.ch2.Image;
    else
        % Skip files that don’t match expected structure
        fprintf('Skip: %s (no recognizable image fields)\n', tmp_name);
        continue
    end

    % Run re-analysis with thresholds (function signature v42 assumes threshold at the end)
    fig_num = 0;
    [cond_img_green, img_ch_green_data_collection] = ...
        condensate_analyzer_confocal_FUNCTION_reanalysis( ...
            img_ch_green, background_setting, peak_background_fold, ...
            tolerance_parameter, fig_num, 'Green Channel', avg_int_threshold);

    [cond_img_red, img_ch_red_data_collection] = ...
        condensate_analyzer_confocal_FUNCTION_reanalysis( ...
            img_ch_red, redch_background_setting, redch_peak_background_fold, ...
            redch_tolerance_parameter, fig_num+1, 'Red Channel', redch_avg_int_threshold);

    % Colocalization v1: Pearson correlation (match app naming)
    cell_green = img_ch_green_data_collection.cell_selc;
    cell_red   = img_ch_red_data_collection.cell_selc;
    R = corrcoef(cell_green, cell_red);
    pearson_cor = R(2);

    % Colocalization v2: overlapped clusters fraction (match app naming)
    matrix1 = cond_img_green; matrix1(cond_img_green>0) = 1;
    matrix2 = cond_img_red;   matrix2(cond_img_red>0)   = 1;

    [labeledMatrix1, numClusters1] = bwlabel(matrix1);
    [labeledMatrix2, ~] = bwlabel(matrix2);

    overlapCount = 0;
    for i = 1:numClusters1
        currentCluster = (labeledMatrix1 == i);
        overlapWithMatrix2 = currentCluster & labeledMatrix2;
        if sum(overlapWithMatrix2(:)) / sum(currentCluster(:)) ~= 0
            overlapCount = overlapCount + 1;
        end
    end
    overlapped_colocal = overlapCount / max(numClusters1, 1);

    % Save out (align with app-style field names)
    time = datestr(now, 'yyyymmddHHMMSS');
    base_name = erase(tmp_name, '.mat');
    save_fname = ['reanalyzed-', base_name, '_', time, '.mat'];

    img_ch_green_data = img_ch_green_data_collection;
    img_ch_red_data   = img_ch_red_data_collection;

    save(save_fname, ...
        'img_ch_green_data','img_ch_red_data', ...
        'overlapped_colocal','pearson_cor', ...
        'background_setting','peak_background_fold','tolerance_parameter', ...
        'redch_background_setting','redch_peak_background_fold','redch_tolerance_parameter', ...
        'avg_int_threshold','redch_avg_int_threshold');

 %   close('all')
end

%% Summarize (match App Designer SummarizeButton output schema)
save_info  = dir('*.mat');
savenames  = {save_info.name};
sname_list = savenames;

final_num_cond        = [];
final_dilute_avg_int  = [];
final_avg_int_cond    = [];
final_sum_int_cond    = [];
final_size_cond       = [];
final_Kp              = [];
final_G_transfer      = [];

final_avg_int_cond_raw = [];
final_sum_int_cond_raw = [];
final_size_cond_raw    = [];
final_Kp_raw           = [];
final_G_transfer_raw   = [];

final_overlapped_colocal = [];
final_pearson_cor         = [];

set_info_background_setting      = [];
set_info_peak_background_fold    = [];
set_info_tolerance_parameter     = [];

set_info_redch_background_setting  = [];
set_info_redch_peak_background_fold= [];
set_info_redch_tolerance_parameter = [];

final_name_collect = {};

for k = 1:length(sname_list)
    tmp_name   = sname_list{k};
    tmp_dat_all= load(tmp_name);

    % only summarize files that contain expected processed fields
    if ~isfield(tmp_dat_all, 'img_ch_green_data') || ~isfield(tmp_dat_all, 'overlapped_colocal')
        continue
    end

    tmp_dat = tmp_dat_all.img_ch_green_data;

    % per-file aggregates (match app)
    final_num_cond(end+1)       = tmp_dat.num_cond;
    final_dilute_avg_int(end+1) = tmp_dat.dilute_average_int;
    final_avg_int_cond(end+1)   = mean(tmp_dat.average_int_cond(:));
    final_sum_int_cond(end+1)   = mean(tmp_dat.sum_int_cond(:));
    final_size_cond(end+1)      = mean(tmp_dat.size_cond(:));
    final_Kp(end+1)             = mean(tmp_dat.Kp(:));
    final_G_transfer(end+1)     = mean(tmp_dat.G_transfer(:));

    final_overlapped_colocal(end+1) = tmp_dat_all.overlapped_colocal;
    final_pearson_cor(end+1)        = tmp_dat_all.pearson_cor;

    % parameter capture (match app field names)
    if isfield(tmp_dat_all,'background_setting')
        set_info_background_setting(end+1)   = tmp_dat_all.background_setting;
        set_info_peak_background_fold(end+1) = tmp_dat_all.peak_background_fold;
        set_info_tolerance_parameter(end+1)  = tmp_dat_all.tolerance_parameter;
    else
        set_info_background_setting(end+1)   = NaN;
        set_info_peak_background_fold(end+1) = NaN;
        set_info_tolerance_parameter(end+1)  = NaN;
    end

    if isfield(tmp_dat_all,'redch_background_setting')
        set_info_redch_background_setting(end+1)   = tmp_dat_all.redch_background_setting;
        set_info_redch_peak_background_fold(end+1) = tmp_dat_all.redch_peak_background_fold;
        set_info_redch_tolerance_parameter(end+1)  = tmp_dat_all.redch_tolerance_parameter;
    else
        set_info_redch_background_setting(end+1)   = NaN;
        set_info_redch_peak_background_fold(end+1) = NaN;
        set_info_redch_tolerance_parameter(end+1)  = NaN;
    end

    % raw distributions
    final_avg_int_cond_raw = [final_avg_int_cond_raw, tmp_dat.average_int_cond(:).'];
    final_sum_int_cond_raw = [final_sum_int_cond_raw, tmp_dat.sum_int_cond(:).'];
    final_size_cond_raw    = [final_size_cond_raw,   tmp_dat.size_cond(:).'];
    final_Kp_raw           = [final_Kp_raw,          tmp_dat.Kp(:).'];
    final_G_transfer_raw   = [final_G_transfer_raw,  tmp_dat.G_transfer(:).'];

    % expand names per-condensate (for raw table)
    final_name_collect = [final_name_collect; repmat({tmp_name}, [tmp_dat.num_cond, 1])];
end

% Summary table (same columns/order as the app)
dat_collection = table( ...
    sname_list.', ...
    final_num_cond.', ...
    final_dilute_avg_int.', ...
    final_avg_int_cond.', ...
    final_sum_int_cond.', ...
    final_size_cond.', ...
    final_Kp.', ...
    final_G_transfer.', ...
    final_overlapped_colocal.', ...
    final_pearson_cor.', ...
    set_info_background_setting.', ...
    set_info_peak_background_fold.', ...
    set_info_tolerance_parameter.', ...
    set_info_redch_background_setting.', ...
    set_info_redch_peak_background_fold.', ...
    set_info_redch_tolerance_parameter.' );

dat_collection.Properties.VariableNames = { ...
    'Fname','numCond','AvgDilute_int','AvgInt','SumInt','Size','Kp','Gtransfer', ...
    'Overlapped_Colocal','Pearson_cor', ...
    'Greench_Background_setting','Greench_Peak_Background_fold','Greench_Tolerance_parameter', ...
    'Redch_Background_setting','Redch_Peak_Background_fold','Redch_Tolerance_parameter'};

% Raw table (per-condensate)
dat_collection_raw = table( ...
    final_name_collect, ...
    final_avg_int_cond_raw.', ...
    final_sum_int_cond_raw.', ...
    final_size_cond_raw.', ...
    final_Kp_raw.', ...
    final_G_transfer_raw.' );

dat_collection_raw.Properties.VariableNames = {'Fname','AvgInt','SumInt','Size','Kp','Gtransfer'};

% Write files with the same names the app uses
writetable(dat_collection,     'summarized.txt',     'Delimiter','tab');
writetable(dat_collection_raw, 'summarized_raw.txt', 'Delimiter','tab');