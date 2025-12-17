clear all;
clc;

afni_dir = '/Volumes/django/EP_EEGfMRI_deconvolution/derivatives/afni';
jsonText = fileread('/Volumes/django/EP_EEGfMRI_deconvolution/code/participants_markers.json'); 
subj_info = jsondecode(jsonText);

subj = fieldnames(subj_info);

hrf_time = 0:0.25:25;
idx_0to8 = find((hrf_time > 0) & (hrf_time < 8));
idx_2to10 = find((hrf_time > 2) & (hrf_time < 10));
idx_20to25 = find((hrf_time > 20) & (hrf_time < 25));

baseline_thre = 0.25;
ceiling_thre = 10;


table_type = table([], {}, {}, [], 'VariableNames', {'subj', 'type', 'type_name', 'sigvox'});

hrf_all = [];

list_subj = 1:length(subj);

for i = list_subj

    isub = strrep(subj{i}, '_', '-');

    hrfs = dir(fullfile(afni_dir, isub, 'TENT_orig_dm/HRF', '*.mat'));
    hrfs = hrfs(~cellfun('isempty', regexp({hrfs.name}, '^type[0-9]\.mat$')));

    for j = 1:length(hrfs)
        matFilePath = fullfile(hrfs(j).folder, hrfs(j).name);

        hrf_tmp = load(matFilePath);
        hrf_tmp = hrf_tmp.hrf_mat_subset;

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

        % Remove the rows from hrf_tmp
        hrf_tmp = hrf_tmp(~rows_to_remove, :);

        if size(hrf_tmp, 1) == 101 && size(hrf_tmp, 2) == 1
            hrf_tmp = hrf_tmp';
        end

        numRows = size(hrf_tmp, 1);
        %rng(19850913);
        %hrf_tmp = hrf_tmp(randperm(numRows, min(50, numRows)), :);

        if size(hrf_tmp, 1) >= 5 && ~strcmpi(strjoin(subj_info.(subj{i}).event_type{j}), 'art-dur') % more than 5 voxels
            hrf_all = [hrf_all; hrf_tmp];
            table_type = [table_type; table({subj{i}}, {hrfs(j).name(1:end-4)}, ...
                {strjoin(subj_info.(subj{i}).event_type{j}', ' ')}, size(hrf_tmp, 1), 'VariableNames', {'subj', 'type', 'type_name', 'sigvox'})];
        end
            
    end
end


% now, using all of the timecourses, perform PCA
[u,s,v] = svd(hrf_all,0);
figure;
plot(hrf_all');
% manually flip, the first PC to be positive, the other 2 are fine, this is
% just to have the same sign as the paper Allen et al., 
if sum(v(:,1)) < 0
    v(:,1) = -v(:,1);
end

if sum(v(idx_0to8,2)) < 0
    v(:,2) = -v(:,2);
end

if sum(v(idx_0to8,3)) < 0
    v(:,3) = -v(:,3);
end



% inspect canonical set of PCs
temp = diag(s);
temp = (temp.^2) / sum(temp.^2);

figure;
h1 = plot(hrf_time, v(:,1), 'r', 'LineWidth', 2); % PC1 in red
hold on;
h2 = plot(hrf_time, v(:,2), 'g', 'LineWidth', 2); % PC2 in green
h3 = plot(hrf_time, v(:,3), 'b', 'LineWidth', 2); % PC3 in blue
yline(0); 
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
subtitle(['Variance Explained: ', num2str(sum(temp(1:3)) * 100, '%.2f'), '%']);
legend([h1, h2, h3], {'PC1', 'PC2', 'PC3'});


%save('/Volumes/django/EP_EEGfMRI_deconvolution/code/hrfbasis.mat','s','v');
