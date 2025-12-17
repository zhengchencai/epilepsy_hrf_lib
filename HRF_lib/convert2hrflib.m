clear all;

addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/GLMdenoise-master'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/GLMsingle'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/knkutils-master'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/nsddatapaper-main'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/TDM-master'));

load('/Volumes/django/EP_EEGfMRI_deconvolution/code/hrf_lib_pts.mat','hrf_lib_pts');

first_half = hrf_lib_pts(1:9, :);
second_half = hrf_lib_pts(10:18, :);

[~, idx_max_first] = max(first_half, [], 2);
[~, idx_max_second] = max(second_half, [], 2);

[~, order_first] = sort(idx_max_first, 'ascend');
[~, order_second] = sort(idx_max_second, 'ascend');

first_half_sorted = first_half(order_first, :);
second_half_sorted = second_half(order_second, :);

hrf_lib_pts = [first_half_sorted; second_half_sorted];


hrf_time = 0:0.25:25;

hrf_lib_pts = hrf_lib_pts ./ max(hrf_lib_pts, [], 2);

varNames = ['hrf_time', compose("hrf%d", 1:size(hrf_lib_pts,1))];
T = array2table([hrf_time', hrf_lib_pts'], 'VariableNames', varNames);
%writetable(T, ['/Volumes/django/EP_EEGfMRI_deconvolution/code/hrflib18.csv']);

for i = 1:18
    figure(i); plot(hrf_time, hrf_lib_pts', 'k');
    hold on;
    plot(hrf_time, hrf_lib_pts(i,:), 'r')
end


% Plot
figure('Color', 'w');

% Top subplot: HRFs 1–9 (red hues)
subplot(2, 1, 1); hold on;
for i = 1:9
    color = [1, (9 - i)/9, (9 - i)/9];  % red hues
    plot(hrf_time, hrf_lib_pts(i, :), 'Color', color, 'LineWidth', 2);
end
xline(0, '--k');
ylim([-0.4 1.1]);  % ensure 0 is visible in y-axis
yticks([0]);
yticklabels({'0'});
set(gca, 'FontSize', 14);
title('HRFs 1–9 (group1: low undershoot)', 'FontSize', 16);
grid on;

% Bottom subplot: HRFs 10–18 (blue hues, reversed)
subplot(2, 1, 2); hold on;
for i = 10:18
    color = [(18 - i)/9, (18 - i)/9, 1];  % blue hues
    plot(hrf_time, hrf_lib_pts(i, :), 'Color', color, 'LineWidth', 2);
end
xline(0, '--k');
ylim([-1.6 1.1]);  % ensure 0 is visible
yticks([0]);
yticklabels({'0'});
xlabel('Time (s)');
set(gca, 'FontSize', 14);
title('HRFs 10–18 (group2: high undershoot)', 'FontSize', 16);
grid on;



afni_dir = '/Volumes/django/EP_EEGfMRI_deconvolution/derivatives/afni';
jsonText = fileread('/Volumes/django/EP_EEGfMRI_deconvolution/code/participants_markers.json');
subj_info = jsondecode(jsonText);

subj = fieldnames(subj_info);

idx_2to10 = find((hrf_time > 2) & (hrf_time < 10));
idx_20to25 = find((hrf_time > 20) & (hrf_time < 25));
idx_10to25 = find((hrf_time > 10) & (hrf_time < 25));


baseline_thre = 0.25;
ceiling_thre = 10;


list_subj = 1:length(subj);

table_hrflib = table([], {}, {}, [], [], [], {}, [], [], [], [], [], [], [],...
    'VariableNames', {'subj', 'type', 'type_name', 'r', 'hrf_id', 'fs_label', 'voxel',...
                      'posPeakAmp', 'posPeakTime', 'posPeakFWHM', 'negPeakAmp', 'negPeakTime', 'negPeakFWHM', 'sign'});

negall = [];

for i = list_subj

    isub = strrep(subj{i}, '_', '-');
    i

    hrfs = dir(fullfile(afni_dir, isub, 'TENT_orig_dm/HRF', '*.mat'));
    hrfs = hrfs(~cellfun('isempty', regexp({hrfs.name}, '^type[0-9]\.mat$')));
    

    for j = 1:length(hrfs)
        matFilePath = fullfile(hrfs(j).folder, hrfs(j).name);

        hrf_mat = load(matFilePath);
        hrf_tmp = hrf_mat.hrf_mat_subset;
        hrf_label_tmp = hrf_mat.selected_voxel_label;
        hrf_voxel_tmp = hrf_mat.selected_voxel;


        if size(hrf_tmp, 2) ~= length(hrf_time)
            hrf_tmp = hrf_tmp';
        end

        % Condition 1: Remove rows where any value is larger than 10
        condition1 = any(abs(hrf_tmp) > ceiling_thre, 2);

        % Condition 2: peak at 2-10s; always increase/decrease within 2-10s
        dhrf = diff(hrf_tmp(:, idx_2to10), 1, 2);
        condition2 = all(dhrf > 0, 2) | all(dhrf < 0, 2);

        % Condition 3: Remove rows where the absolute value of the first
        % element is larger than 0.5 (start from close to 0)
        condition3 = abs(hrf_tmp(:, 1)) > baseline_thre;

        % Condition 4: back to baseline during 20  to 25s
        condition4 = all(abs(hrf_tmp(:, idx_20to25)) > baseline_thre, 2);

        % Combine both conditions (rows that meet either condition)
        rows_to_remove = condition1 | condition2 | condition3 | condition4;

        rows_to_remove = condition1;

        % Remove the rows from hrf_tmp
        hrf_tmp = hrf_tmp(~rows_to_remove, :);
        hrf_label_tmp = hrf_label_tmp(~rows_to_remove, :);
        hrf_voxel_tmp(rows_to_remove) = [];

        if size(hrf_tmp, 1) == 101 && size(hrf_tmp, 2) == 1
            hrf_tmp = hrf_tmp';
        end

        numRows = size(hrf_tmp, 1);
        %rng(19850913);
        %hrf_tmp = hrf_tmp(randperm(numRows, min(500, numRows)), :);
        
        n = size(hrf_tmp, 1);
        best_idx = zeros(n, 1);
        best_r = zeros(n, 1);
        posPeakAmp = zeros(n, 1);
        posPeakTime = zeros(n, 1);
        posPeakFWHM = zeros(n, 1);
        negPeakAmp = zeros(n, 1);
        negPeakTime = zeros(n, 1);
        negPeakFWHM = zeros(n, 1);
        sign = ones(n, 1);
        for k = 1:size(hrf_tmp, 1)

            hrf_tmp_one = hrf_tmp(k, :);
            [~, idx_tmp] = min(abs(diff(hrf_tmp_one(idx_2to10)))); 
            
            % if negative peak within 2 to 10s, but whole area under HRF is positive, it is positive HRF

            if hrf_tmp_one(idx_2to10(idx_tmp)) < 0 && sum(hrf_tmp_one) < 0
               negall = [negall; hrf_tmp_one];
               hrf_tmp_one = -hrf_tmp_one;
               sign(k) = -1;
            end    
           

            r = corr(hrf_tmp_one', hrf_lib_pts');           
            [best_r(k), best_idx(k)] = max(abs(r));   


            % extract HRF features
            [yPosPeaks, locsPos, widthsPos] = findpeaks(hrf_tmp_one, 'WidthReference', 'halfheight');
            [~, maxIdx] = max(yPosPeaks);
            if isempty(maxIdx)
                posPeakAmp(k) = -1000;
                posPeakTime(k) = -1000;
                posPeakFWHM(k) = -1000;
            else
                posPeakAmp(k) = yPosPeaks(maxIdx);
                posPeakTime(k) = hrf_time(locsPos(maxIdx));
                posPeakFWHM(k) = widthsPos(maxIdx) * (hrf_time(2) - hrf_time(1));

            end

            [yNegPeaks, locsNeg, widthsNeg] = findpeaks(-hrf_tmp_one, 'WidthReference', 'halfheight');
            [~, minIdx] = max(yNegPeaks); 
            if isempty(minIdx)
                negPeakAmp(k) = -1000;
                negPeakTime(k) = -1000;
                negPeakFWHM(k) = -1000;
            else
                negPeakAmp(k) = -yNegPeaks(minIdx);
                negPeakTime(k) = hrf_time(locsNeg(minIdx));
                negPeakFWHM(k) = widthsNeg(minIdx) * (hrf_time(2) - hrf_time(1));
            end
        end

        if ~strcmpi(strjoin(subj_info.(subj{i}).event_type{j}), 'art-dur') 
            num_vox = numel(best_r);
            table_hrflib = [table_hrflib; table( ...
                repmat(subj(i), num_vox, 1), ...
                repmat({hrfs(j).name(1:end-4)}, num_vox, 1), ...
                repmat({strjoin(subj_info.(subj{i}).event_type{j}', ' ')}, num_vox, 1), ...
                best_r, ...
                best_idx, ...
                hrf_label_tmp, ...
                hrf_voxel_tmp, ...
                posPeakAmp, ...
                posPeakTime, ...
                posPeakFWHM, ...
                negPeakAmp, ...
                negPeakTime, ...
                negPeakFWHM, ...
                sign, ...
                'VariableNames', {'subj', 'type', 'type_name', 'r', 'hrf_id', 'fs_label', 'voxel', ...
                'posPeakAmp', 'posPeakTime', 'posPeakFWHM', 'negPeakAmp', 'negPeakTime', 'negPeakFWHM', 'sign'})];

        end
    end
end

%writetable(table_hrflib, ['/Volumes/django/EP_EEGfMRI_deconvolution/code/table_hrflib', num2str(hrf_num), '_withsign.csv']);
%writetable(table_hrflib, '/Volumes/django/EP_EEGfMRI_deconvolution/code/table_hrflib.csv');