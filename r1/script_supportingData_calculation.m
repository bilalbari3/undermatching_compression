%% load data
clearvars -except fs_b data
if exist('fs_b','var') == 0
    load('C:\Users\bilal\Dropbox\Research\data\analysis\dynamicForaging\revision\fullStructures\qMods_b.mat')
    fprintf('Bari 2019 data loaded\n')
end
if exist('data_table', 'var') == 0
    load('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\matching_data.mat')
    fprintf('Matching data loaded \n')
end

%% ITI calculation
allITI = [];
mask_405 = contains(data_table.Task, '40/5');
mask_4010 = contains(data_table.Task, '40/10');
mask_task = mask_405 | mask_4010;

all_sessions = data_table(mask_task, :).Session_Name';
for sess = all_sessions
    sess = sess{:};
    a = sess(2:6);
    allITI = [allITI fs_b.(a).(sess).pd.fbITI];
end
allITI(allITI > 30) = [];

nanmean(allITI)

%% number of mice, etc
mask_405 = contains(data_table.Task, '40/5');
mask_4010 = contains(data_table.Task, '40/10');
mask_mult = contains(data_table.Task, 'mult');

mice_all = data_table(mask_405 | mask_4010 | mask_mult, :).Mouse;
mice_405 = data_table(mask_405, :).Mouse;
mice_4010 = data_table(mask_4010, :).Mouse;
mice_mult = data_table(mask_mult, :).Mouse;

fprintf(['Total: %i mice %i sessions\n', ...
         '4050: %i mice %i sessions\n', ...
         '4010: %i mice %i sessions\n', ...
         'Mult: %i mice %i sessions\n'], ...
    length(unique(mice_all)), length(mice_all), ...
    length(unique(mice_405)), length(mice_405), ...
    length(unique(mice_4010)), length(mice_4010), ...
    length(unique(mice_mult)), length(mice_mult))

allTrials = data_table(mask_405 | mask_4010 | mask_mult, :).Trials;
fprintf('Median trials: %i (min: %i, max: %i)\n', median(allTrials), min(allTrials), max(allTrials))

%%
group_titles = {'PD off Rx', 'PD on Rx'};
group_mat = zeros(length(data), 2);
rfrac_log2_all = [];
cfrac_log2_all = [];
for g_ind = 3:4
    if g_ind == 3
        s_ind = 1;
    end
    group_mat(s_ind:length(summ.(groups{g_ind}).cfrac) + s_ind - 1, g_ind - 2) = 1; % matrix of group identity [yc ec pdOFF pdON];
    s_ind = length(summ.(groups{g_ind}).cfrac) + s_ind;

    rfrac_log2_all = [rfrac_log2_all; summ.(groups{g_ind}).rfrac_log2']; % concatenate rfrac and cfrac
    cfrac_log2_all = [cfrac_log2_all; summ.(groups{g_ind}).cfrac_log2'];
end
cfrac_log2_all(cfrac_log2_all > 4) = NaN;
%%
mod = fitlm([rfrac_log2_all group_mat(:, 2).*rfrac_log2_all], cfrac_log2_all)
% mod_coef = mod.Coefficients.Estimate(2:end);
% mod_CI = coefCI(mod); mod_CI = mod_CI(2:end, :);