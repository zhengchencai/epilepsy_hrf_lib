clear all;

addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/GLMdenoise-master'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/GLMsingle'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/knkutils-master'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/nsddatapaper-main'));
addpath(genpath('/Volumes/django/EP_EEGfMRI_deconvolution/code/external_tool/TDM-master'));

load('/Volumes/django/EP_EEGfMRI_deconvolution/code/hrfbasis.mat','s','v');

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

normalize_by_vectorlength = 0;

opt.bins = -1.5:.05:1.5;
dif0 = diff(opt.bins(1:2));
opt.binsEDGES = linspacefixeddiff(opt.bins(1)-dif0/2,dif0,length(opt.bins)+1);
numbins = length(opt.bins);

hrf_all = [];

ns = [];
allpoints = {};

list_subj = 1:length(subj);

pos_num = 0;
neg_num = 0;
negall = [];

for i = list_subj

    isub = strrep(subj{i}, '_', '-');
    i

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
        rng(19850913);
        %hrf_tmp = hrf_tmp(randperm(numRows, min(300, numRows)), :);
        
        % flip 
        sum_within_idx = sum(hrf_tmp(:, idx_0to8), 2);

        negall = [negall; hrf_tmp(sum_within_idx < 0, :)];
        hrf_tmp(sum_within_idx < 0, :) = -hrf_tmp(sum_within_idx < 0, :);
        
        neg_num = neg_num + sum(sum_within_idx < 0);
        pos_num = pos_num + sum(sum_within_idx >= 0);

        % compute density
        if size(hrf_tmp, 1) >= 5 && ~strcmpi(strjoin(subj_info.(subj{i}).event_type{j}), 'art-dur')
            
            % intensity
            m0 = hrf_tmp * v(:, 1:3);
            m0l = vectorlength(m0,2);
            m0 = unitlength(m0,2);

            % get points
            XXX0 = m0(:,2);
            YYY0 = m0(:,3);
            ZZZ0 = m0(:,1);

            % calc
            [ndensity,xx,yy] = hist2d(XXX0,YYY0,opt.bins,opt.bins);

            % vector length
            nvectorlength = zeros(numbins,numbins);
            for cc=1:numbins
              ok = find(XXX0 > opt.binsEDGES(cc) & XXX0 <= opt.binsEDGES(cc+1));
              for rr=1:numbins
                ok2 = find(YYY0(ok) > opt.binsEDGES(rr) & YYY0(ok) <= opt.binsEDGES(rr+1));
                nvectorlength(rr,cc) = median(m0l(ok(ok2)));
              end
            end

            % normalize_by_vectorlength
            if normalize_by_vectorlength == 1
                nvectorlength = normalizerange(nvectorlength,0,1,0);
                nvectorlength(isnan(nvectorlength)) = 0;
                ndensity = ndensity.* nvectorlength;
            end

            allpoints{i} = [XXX0 YYY0 ZZZ0];
            ns(:,:,i) = ndensity/sum(ndensity(:));
            hrf_all = [hrf_all; hrf_tmp];

        end
    end
end


%% voxel wise
m0 = hrf_all * v(:, 1:3);
m0l = vectorlength(m0,2);
m0 = unitlength(m0,2);

% get points
XXX0 = m0(:,2);
YYY0 = m0(:,3);
ZZZ0 = m0(:,1);

% calc
[ndensity,xx,yy] = hist2d(XXX0,YYY0,opt.bins,opt.bins);

% vector length
nvectorlength = zeros(numbins,numbins);
for cc=1:numbins
  ok = find(XXX0 > opt.binsEDGES(cc) & XXX0 <= opt.binsEDGES(cc+1));
  for rr=1:numbins
    ok2 = find(YYY0(ok) > opt.binsEDGES(rr) & YYY0(ok) <= opt.binsEDGES(rr+1));
    nvectorlength(rr,cc) = median(m0l(ok(ok2)));
  end
end

% normalize_by_vectorlength
if normalize_by_vectorlength == 1
    nvectorlength = normalizerange(nvectorlength,0,1,0);
    nvectorlength(isnan(nvectorlength)) = 0;
    ndensity = ndensity.* nvectorlength;
end


% visualize
fig = figure; hold on;
imagesc(xx(1,:),yy(:,1),ndensity);
colormap(hot);
axis equal tight;
xlabel('loading on PC2');
ylabel('loading on PC3');
drawellipse(0,0,0,1,1,[],[],'w-');


% manually define high density points
% % figure; hold on;
% % imagesc(xx(1,:),yy(:,1),ndensity);
% % colormap(hot);
% % axis equal tight;
% % %      straightline(0,'h','w-');
% % %      straightline(0,'v','w-');
% % xlabel('loading on PC2');
% % ylabel('loading on PC3');
% % drawellipse(0,0,0,1,1,[],[],'w-');
% % [x1,y1] = ginput(9);
% % 
% % 
% % figure; hold on;
% % imagesc(xx(1,:),yy(:,1),ndensity);
% % colormap(hot);
% % axis equal tight;
% % %      straightline(0,'h','w-');
% % %      straightline(0,'v','w-');
% % xlabel('loading on PC2');
% % ylabel('loading on PC3');
% % drawellipse(0,0,0,1,1,[],[],'w-');
% % [x2,y2] = ginput(9);



%save('/Volumes/django/EP_EEGfMRI_deconvolution/code/xypoints.mat','x1', 'x2', 'y1', 'y2');
load('/Volumes/django/EP_EEGfMRI_deconvolution/code/xypoints.mat'); 

xypoint1 = [x1, y1];
xypoint2 = [x2, y2];

fig = figure('Units','inches', 'Position',[1 1 5.5 5.5]); hold on;
imagesc(xx(1,:),yy(:,1),ndensity);
colormap(hot);
hold on;
scatter(xypoint1(:,1), xypoint1(:,2), 200, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.6);
scatter(xypoint2(:,1), xypoint2(:,2), 200, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.6);

axis equal tight;
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
xlabel('loading on PC2');
ylabel('loading on PC3');
drawellipse(0,0,0,1,1,[],[],'w-');
exportgraphics(fig, fullfile('/Volumes/django/EP_EEGfMRI_deconvolution/code/tmp/hist2d', sprintf('hist2dall.png')), 'Resolution', 300);


% uniformly distribute N particles across the surface of a unit sphere.
% let's take V(:,3) to be the inward dimension
[V,Tri,~,Ue] = ParticleSampleSphere('N',100*2);
fv=struct('faces',Tri,'vertices',V); 
fv_new=SubdivideSphericalMesh(fv,2);   % subdivide to become denser
V = fv_new.vertices;
V = V(V(:,3)>0,:);  % T x 3

% how dense should we make it??
angularsep = 4;

thepath = [];
pts = [];

xypoint = [x1, y1];

% make nice unit-sphere points (1st column is inward, 2nd and 3rd columns are x- and y-)
pts_tmp = [sqrt(1-sum(xypoint.^2,2)) xypoint];  % clicked-points x 3

% construct path
thepath_tmp = constructpathonsphere(pts_tmp,angularsep);  % new-points x 3

hrf_lib1 = thepath_tmp*v(:, 1:3)';

thepath = [thepath; thepath_tmp];
pts = [pts; pts_tmp];


xypoint = [x2, y2];

% make nice unit-sphere points (1st column is inward, 2nd and 3rd columns are x- and y-)
pts_tmp = [sqrt(1-sum(xypoint.^2,2)) xypoint];  % clicked-points x 3

% construct path
thepath_tmp = constructpathonsphere(pts_tmp,angularsep);  % new-points x 3

thepath = [thepath; thepath_tmp];
pts = [pts; pts_tmp];


hrf_lib2 = thepath_tmp*v(:, 1:3)';

hrf_lib_pts = pts*v(:, 1:3)';

figure; 
plot(hrf_time, (hrf_lib1 ./ max(hrf_lib1, [], 2))', 'r');
xlabel("Time(s)")
hold on;
plot(hrf_time, (hrf_lib2 ./ max(hrf_lib2, [], 2))', 'b');
hold off;
ylim([-1.5, 1.5]);

figure; 
plot(hrf_time, (hrf_lib_pts ./ max(hrf_lib_pts, [], 2))', 'k');
ylim([-1.5, 1.5]);

% save('/Volumes/django/EP_EEGfMRI_deconvolution/code/hrf_lib_pts.mat','hrf_lib_pts');
