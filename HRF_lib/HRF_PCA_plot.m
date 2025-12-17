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
ceiling_list = [10, 5, 3];

drawlist = [50, 100, 300, 500, 1000000];

table_type = table([], {}, {}, 'VariableNames', {'subj', 'type', 'type_name'});

pca_table = table([], [], [], [], [], [], [], 'VariableNames', {'pc1', 'pc2', 'pc3', 'ceil', 'totvar', 'draw', 'idraw'});

list_subj = 1:length(subj);

for ceiling_thre = ceiling_list
    ceiling_thre

    for mindraw = drawlist

        mindraw

        if mindraw == 1000000 
            draws = 1;
        else
            draws = 1:1:500;
        end

        for idraw = draws

            hrf_all = [];

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
                    if mindraw ~= 1000000
                        hrf_tmp = hrf_tmp(randperm(numRows, min(mindraw, numRows)), :);
                    end


                    if size(hrf_tmp, 1) >= 5 && ~strcmpi(strjoin(subj_info.(subj{i}).event_type{j}), 'art-dur') % more than 5 voxels
                        hrf_all = [hrf_all; hrf_tmp];
                    end

                end
            end

            % now, using all of the timecourses, perform PCA
            [u,s,v] = svd(hrf_all,0);

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

            var_exp = sum(temp(1:3)) * 100;

            pca_table = [pca_table; table(v(:,1), v(:,2), v(:,3),...
                repmat(ceiling_thre, [length(v(:,1)), 1]),...
                repmat(var_exp, [length(v(:,1)), 1]),...
                repmat(mindraw, [length(v(:,1)), 1]),...
                repmat(idraw, [length(v(:,1)), 1]), 'VariableNames', {'pc1', 'pc2', 'pc3', 'ceil', 'totvar', 'draw', 'idraw'})];

        end

    end

end


%writetable(pca_table, '/Volumes/django/EP_EEGfMRI_deconvolution/code/PCA_alldraws.csv');



