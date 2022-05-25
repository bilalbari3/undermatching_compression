%% load data
if exist('data','var') == 0 % not loaded
    load('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\rutledge09_data.mat')
    fprintf('Rutledge 2009 data loaded\n')
end
myColors = importColors_bb();
c_group.yc = myColors.darkGray;
c_group.ec = myColors.lightGray;
c_group.pdOFF = myColors.reddishPurple_dull;
c_group.pdON = myColors.reddishPurple_bright;
c_group.pd = myColors.reddishPurple;

saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\figures\rutledge09';
% N - # trials
% C - # options (always 2)
% pop - population (1 = young controls, 2 = elderly controls, 3 = PD off meds, 4 = PD on meds)
% block - block #
% trial - trial # within block
% action - 1 = left, 2 = right
% reward - 0 = not rewarded, 1 = rewarded
% RT - response time (msec)
% acc - 1 if subject chose highest probability arm, 0 otherwise
% state - hidden state of the task, correspondindag to different reward contingencies (specified in the "prob" field), which change across blocks (note that the block changes are unsignaled):
%               s1: 1/7, 6/7
%               s2: 1/4, 3/4
%               s3: 3/4, 1/4
%               s4: 6/7, 1/7

% trim data for one PD subject who did not perform as many trials PD_on
if length(data) == 104
    s_all = [65 91]; % this subject did not perform as many trials PD_on (ind 91)    data(s_all) = [];
    data(s_all) = [];
    fprintf('Dropped one pd subject\n')
end

%% munge the data

trial_start = 21; % mimic analysis of paper
trial_iter_mask = false; % if true, go to trial_iter + trial_start
trial_iter = 49;

groups = {'yc','ec','pdOFF','pdON'};
fields = {'cfrac','rfrac','bits','RT','RT_session','um_slope'};
clearvars summ
for g = groups
    for f = fields
        summ.(g{:}).(f{:}) = [];
    end
end
for s = 1:length(data)
    block = data(s).block'; % block data

    tmp_cfrac = [];
    tmp_rfrac = [];
    used_inds = []; % use to calculate mutual information
    tmp_RT = [];
    tmp_state = []; % for picking RTs in each state
    for curr_b = unique(block) % iterate through each block
        curr_b_mask = block == curr_b; % mask that block

        choice = data(s).action(curr_b_mask) == 1; % all actions to red option in that block
        reward = data(s).reward(curr_b_mask);
        RT = data(s).RT(curr_b_mask);
        
        tmp_state = [tmp_state; unique(data(s).state(curr_b_mask))];
        % remove burn-in trials
        if trial_iter_mask == true
            choice = choice(trial_start:trial_iter + trial_start);
            reward = reward(trial_start:trial_iter + trial_start);
        else
            choice = choice(trial_start:end);
            reward = reward(trial_start:end);
        end

        % calculate choice and reward fractions
        tmp_cfrac = [tmp_cfrac sum(choice == 1)/sum(choice == 0)];
        tmp_rfrac = [tmp_rfrac sum(reward(choice == 1))/sum(reward(choice == 0))];

        % get indices to calculate mutual information
        tmp_inds = find(curr_b_mask); % use to calculate mutual information
        if trial_iter_mask == true
            tmp_inds = tmp_inds(trial_start:trial_iter + trial_start);
        else
            tmp_inds = tmp_inds(trial_start:end);
        end
        used_inds = [used_inds tmp_inds];

        % get reaction times
        if trial_iter_mask == true
            tmp_RT = [tmp_RT median(RT(trial_start:trial_iter + trial_start))];
        else
            tmp_RT = [tmp_RT median(RT(trial_start:end))];
        end
    end
    % calculate mutual information
    choice = data(s).action(used_inds);
    state = data(s).state(used_inds);
    pi_base = [mean(choice == 1) mean(choice == 2)]; % marginal distribution of actions
    p_state = [mean(state == 1) mean(state == 2) mean(state == 3) mean(state == 4)]; % distribution of states
    pi = [mean(choice(state == 1) == 1) mean(choice(state == 1) == 2) % joint distribution of actions and states
          mean(choice(state == 2) == 1) mean(choice(state == 2) == 2)
          mean(choice(state == 3) == 1) mean(choice(state == 3) == 2)
          mean(choice(state == 4) == 1) mean(choice(state == 4) == 2)];
    tmp_bits = p_com(pi, pi_base, p_state); % mutual information

%     tmp_bits = mutual_information(data(s).state(used_inds), data(s).action(used_inds), 0.1);

    pop = data(s).pop;

    summ.(groups{pop}).cfrac = [summ.(groups{pop}).cfrac tmp_cfrac];
    summ.(groups{pop}).rfrac = [summ.(groups{pop}).rfrac tmp_rfrac];
    summ.(groups{pop}).bits = [summ.(groups{pop}).bits tmp_bits];

    % RT analyses
    summ.(groups{pop}).RT_session = [summ.(groups{pop}).RT_session median(data(s).RT(used_inds))];

    % count each state once
    state_tmp = data(s).state(used_inds);
    RT_tmp = data(s).RT(used_inds);
    RT_out = [];
    for s_ind = unique(state_tmp')
        RT_out = [RT_out median(RT_tmp(state_tmp == s_ind))];
    end
    summ.(groups{pop}).RT = [summ.(groups{pop}).RT RT_out];
    % within-session undermatching slope
    tmp_rfrac = log2(tmp_rfrac);
    tmp_rfrac(isinf(tmp_rfrac)) = NaN;
    tmp_cfrac = log2(tmp_cfrac);
    tmp_mod = fitlm(tmp_rfrac, tmp_cfrac);
    summ.(groups{pop}).um_slope = [summ.(groups{pop}).um_slope tmp_mod.Coefficients.Estimate(2)];
end

for g = groups
    g = g{:};
    % clean up reward fraction
    summ.(g).rfrac_log2 = log2(summ.(g).rfrac);
    summ.(g).rfrac_log2(isinf(summ.(g).rfrac_log2)) = NaN;

    summ.(g).cfrac_log2 = log2(summ.(g).cfrac);
end
%% fit linear regression
group_titles = {'Young controls', 'Elderly controls', 'PD off Rx', 'PD on Rx'};
group_mat = zeros(length(data), 4);
rfrac_log2_all = [];
cfrac_log2_all = [];
for g_ind = 1:4
    if g_ind == 1
        s_ind = 1;
    end
    group_mat(s_ind:length(summ.(groups{g_ind}).cfrac) + s_ind - 1, g_ind) = 1; % matrix of group identity [yc ec pdOFF pdON];
    s_ind = length(summ.(groups{g_ind}).cfrac) + s_ind;

    rfrac_log2_all = [rfrac_log2_all; summ.(groups{g_ind}).rfrac_log2']; % concatenate rfrac and cfrac
    cfrac_log2_all = [cfrac_log2_all; summ.(groups{g_ind}).cfrac_log2'];
end
mod = fitlm(group_mat.*rfrac_log2_all, cfrac_log2_all)
mod_coef = mod.Coefficients.Estimate(2:end);
mod_CI = coefCI(mod); mod_CI = mod_CI(2:end, :);

h_rwdSensitivity = figure('units','normalized','outerposition',[0 0.2 0.4 0.75]);
h_scatter_sp(1) = subplot(221); hold on
h_scatter_sp(2) = subplot(222); hold on
h_scatter_sp(3) = subplot(223); hold on
h_scatter_sp(4) = subplot(224); hold on

h_coef_polComp = figure('units','normalized','outerposition',[0 0.1 0.45 0.75]);
h_coef = subplot(221); hold on
h_polComp_hist = subplot(222); hold on
h_withinsubj = subplot(223); hold on
h_RT_hist = subplot(224); hold on
x_pred = [-5 5];
for g_ind = 1:4
    figure(h_rwdSensitivity)
    subplot(h_scatter_sp(g_ind)); axis(h_scatter_sp(g_ind), 'square'); xlim([-5 5]); ylim([-5 5])
    xlabel('$\textrm{log}_2(\frac{\textrm{R}_{\textrm{red}}}{\textrm{R}_{\textrm{green}}})$','Interpreter','Latex')
    ylabel('$\textrm{log}_2(\frac{\textrm{C}_{\textrm{red}}}{\textrm{C}_{\textrm{green}}})$','Interpreter','Latex')
    title(group_titles{g_ind})
    plot([-5 5],[-5 5],':','Color', myColors.darkGray, 'linewidth', 1.3)
    scatter(summ.(groups{g_ind}).rfrac_log2, summ.(groups{g_ind}).cfrac_log2, 5, 'filled', 'MarkerFaceColor', c_group.(groups{g_ind}))
    g_pred = [0 0 0 0]; g_pred(g_ind) = 1;
    plot([x_pred], mod.predict(x_pred'.*g_pred), 'k', 'linewidth', 1)

    figure(h_coef_polComp)
    subplot(h_coef)
    errorbar(g_ind, mod_coef(g_ind), mod_coef(g_ind) - mod_CI(g_ind, 1), mod_coef(g_ind) - mod_CI(g_ind, 2), 'o', 'linewidth', 2, ...
        'Color', c_group.(groups{g_ind}))
end
xlim([0.5 4.5])
ylim([0.15 0.5])
set(gca, 'xtick', 1:4, 'xticklabel', group_titles)
xlabel('Group')
ylabel(sprintf('Matching\n(slope of reward x choice)'))

% policy complexity histogram
subplot(h_polComp_hist)
plot([0 0],[0 1],'k:')
histogram(summ.pdON.bits - summ.pdOFF.bits, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', c_group.pd, 'FaceAlpha', 1)
xlim([-0.2 0.2])
ylim([0 0.6])
xlabel(sprintf('Change in policy complexity\nON minus OFF (bits)'))
ylabel('Probability')
set(gca, 'tickdir', 'out', 'ytick', 0:0.1:0.6)

% policy complexity scatter
subplot(h_withinsubj)
scatter(summ.pdON.bits - summ.pdOFF.bits, summ.pdON.um_slope - summ.pdOFF.um_slope, 10, 'filled', 'MarkerFaceColor', c_group.pd)
xlabel(sprintf('Change in policy complexity\nON minus OFF (bits)'))
ylabel(sprintf('Change in matching\nON minus OFF (arb)'))
xlim(h_withinsubj, [-0.15 0.15])
fitlm(summ.pdON.bits - summ.pdOFF.bits, summ.pdON.um_slope - summ.pdOFF.um_slope)
set([h_scatter_sp h_coef h_withinsubj],'tickdir','out')

% stats
[~,p] = lillietest(summ.pdON.bits - summ.pdOFF.bits)
[p,h,stats] = signrank(summ.pdON.bits, summ.pdOFF.bits)

% response time histogram
subplot(h_RT_hist)
plot([0 0],[0 1],'k:')
bins = -550:100:550;
histogram(summ.pdON.RT - summ.pdOFF.RT, bins, 'Normalization', 'Probability', ...
          'EdgeColor', 'none', 'FaceColor', c_group.pd, 'FaceAlpha', 1)
xlabel(sprintf('Change in response time\nON minus OFF (ms)'))
ylabel('Probability')
xlim([-550 550])
ylim([0 0.4])
set(gca, 'tickdir', 'out', 'ytick', 0:0.1:0.4, 'xtick', [-500:250:500])

% stats
[~,p] = lillietest(summ.pdON.RT - summ.pdOFF.RT)
[p,h,stats] = signrank(summ.pdON.RT, summ.pdOFF.RT)

saveFigureIteration(h_rwdSensitivity, saveFigLoc, 'rwdSensitivity')
saveFigureIteration(h_coef_polComp, saveFigLoc, 'policyComplexity')