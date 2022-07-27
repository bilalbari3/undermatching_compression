%% load data
clearvars -except fs_b
if exist('fs_b','var') == 0
    load('C:\Users\bilal\Dropbox\Research\data\analysis\dynamicForaging\revision\fullStructures\qMods_b.mat')
end
myColors = importColors_bb();

% find relevant sessions and do a whole lot of cleanup
sess_40_10 = {};
sess_40_7 = {};
sess_40_5 = {};
sess_50_5 = {};
sess_mult = {};
all_m = fields(fs_b)';
all_bp = {};
for m = all_m
    m = m{:};
    all_sess = fields(fs_b.(m))';
    for sess = all_sess
        sess = sess{:};
        if strcmp(sess, 'mBB115d20181230') % session that needs FIXING
            fprintf('Skipping %s\n', sess)
        else
            unique_bp = unique(fs_b.(m).(sess).pd.bp);
            if length(unique_bp) == 2 & any(contains(unique_bp, '40/10'))
                sess_40_10 = [sess_40_10 {sess}];
            end
            if length(unique_bp) == 2 & any(contains(unique_bp, '40/7'))
                sess_40_7 = [sess_40_7 {sess}];
            end
            if length(unique_bp) == 2 & any(contains(unique_bp, '40/5'))
                sess_40_5 = [sess_40_5 {sess}];
            end
            if length(unique_bp) == 2 & any(contains(unique_bp, '50/5'))
                sess_50_5 = [sess_50_5 {sess}];
            end
            if length(unique_bp) > 2 & ...
                    (any(contains(unique_bp, '40.00')) | ...
                     any(contains(unique_bp, '5.00')) | ...
                     any(contains(unique_bp, '38.75')) | ...
                     any(contains(unique_bp, '6.43')) | ...
                     any(contains(unique_bp, '6')) | ...
                     any(contains(unique_bp, '33.75')) | ...
                     any(contains(unique_bp, '11.25')) | ...
                     any(contains(unique_bp, '11')) | ...
                     any(contains(unique_bp, '22.50')) | ...
                     any(contains(unique_bp, '22'))) & ...
                    (any(contains(unique_bp, '30')) == 0)
                 sess_mult = [sess_mult {sess}];
                 all_bp = [all_bp {unique_bp}];
            end
        end
    end
end

%% munge all thems data
ind_offset = 21; % trials to reach steady state for each block

summ = struct();
summ.rwd = [];
summ.mouse_name = {};
summ.task_name = {};
summ.bits = [];
summ.cond_ent = [];
summ.session_actual = []; % actual session number
summ.session_reset = []; % only iterate when that session type reached again
summ.allRT = [];
summ.nTrials = [];
summ.false_alarm = {};

% for computing undermatching
summ.um_soltani = []; % by block
summ.um_slope_session = [];
summ.um_task_block = {};
summ.cfrac = [];
summ.rfrac = [];

summ.setsize = []; % for doing set size comparisons

all_m = fields(fs_b)';
for state = {'40/10','40/7','40/5'}
    state = state{:};
    switch state
        case '40/10'
            relev_sess = sess_40_10;
        case '40/7'
            relev_sess = sess_40_7;
        case '40/5'
            relev_sess = sess_40_5;
    end
    for m = all_m
        m = m{:};
        all_sess = fields(fs_b.(m))';

        sess_ind = 0; % initialize
        sess_reset = 0;
        for sess = all_sess
            sess = sess{:};
            sess_ind = sess_ind + 1;
            if any(strcmp(relev_sess, sess)) % within each session type (40/10, 40/7, 40/5)

                allC_R_hi = []; % choices to rightward side when it is the richer option
                allC_R_lo = []; % choices to the rightward side when it is the leaner option
                all_hi_state = []; % 1 if the leftward side is richer (for computing MI)

                allC_R = fs_b.(m).(sess).pd.allC_R; % all rightward choices
                allC_L = fs_b.(m).(sess).pd.allC_L; % all leftward choices
                allR = abs(fs_b.(m).(sess).pd.allR); % all rewards
                allR_R = fs_b.(m).(sess).pd.allR_R; % all rewards from right
                allR_L = fs_b.(m).(sess).pd.allR_L; % all rewards from left
                bs_all = fs_b.(m).(sess).pd.bs_corrected; bs_all = [bs_all length(allC_R)+1]; % all block switches
                bp = fs_b.(m).(sess).pd.bp; % all block probabilities
                cfrac_tmp = []; % for storing single session cfrac values
                rfrac_tmp = [];
                used_inds = [];
                for bs_ind = 1:length(bs_all) - 1 % iterate through each block
                    inds_of_int = bs_all(bs_ind):bs_all(bs_ind + 1) - 1; % all indices of that block
                    
                    % remove early trials to capture steady state
                    if length(inds_of_int) > ind_offset
                        inds_of_int = inds_of_int(ind_offset:end); % remove trials 1:ind_offset-1

                        if contains(bp{bs_ind}, state) % if left side is richer (since state is 40/x and not x/40)
                            allC_R_lo = [allC_R_lo allC_R(inds_of_int)]; % R choices towards low probability side
                            all_hi_state = [all_hi_state ones(1, length(inds_of_int))];
                        else
                            allC_R_hi = [allC_R_hi allC_R(inds_of_int)]; % R choices to high probability side
                            all_hi_state = [all_hi_state zeros(1, length(inds_of_int))];
                        end
                        used_inds = [used_inds inds_of_int]; % for Hutter estimator
                        summ.um_task_block = [summ.um_task_block {state}]; % task for each block
                        % compute soltani-style undermatching (not used)
                        choice_f = sum(allC_L(inds_of_int))/length(allC_L(inds_of_int));
                        reward_f = sum(fs_b.(m).(sess).pd.allR_L(inds_of_int))/ ...
                                   sum(allR(inds_of_int));

                        if sign(reward_f - 0.5) ~= 0
                            tmp_undermatching = (choice_f - reward_f)*sign(reward_f - 0.5);
                        else
                            tmp_undermatching = (choice_f - reward_f);
                        end
                        summ.um_soltani = [summ.um_soltani tmp_undermatching];

                        % store cfrac and rfrac values, within session, for within session undermatching slope
                        cfrac_tmp = [cfrac_tmp sum(allC_R(inds_of_int)) / sum(allC_L(inds_of_int))];
                        rfrac_tmp = [rfrac_tmp sum(allR_R(inds_of_int)) / sum(allR_L(inds_of_int))];
                    end
                end
                % compute generalized matching style undermatching (for regression, block-by-block values)
                summ.cfrac = [summ.cfrac cfrac_tmp];
                summ.rfrac = [summ.rfrac rfrac_tmp];

                % calculate session-wide undermatching slope
                rfrac_log2_tmp = log2(rfrac_tmp);
                rfrac_log2_tmp(isinf(rfrac_log2_tmp)) = NaN;
                cfrac_log2_tmp = log2(cfrac_tmp);
                um_mod = fitlm(rfrac_log2_tmp, cfrac_log2_tmp);
                summ.um_slope_session = [summ.um_slope_session um_mod.Coefficients.Estimate(2)];

                % calculate MI
                pi_base = [1-mean(allC_R(used_inds)) mean(allC_R(used_inds))];
                p_state = [1-mean(all_hi_state) mean(all_hi_state)];
                pi = [1-mean(allC_R_hi) mean(allC_R_hi);   % joint [state x action] distribution
                      1-mean(allC_R_lo) mean(allC_R_lo) ];

                [tmp_bits, ~, tmp_cond_ent] = p_com(pi, pi_base, p_state);
                summ.bits = [summ.bits tmp_bits];
                summ.cond_ent = [summ.cond_ent tmp_cond_ent];

%                 % calculate MI using Hutter estimator
%                 tmp_bits = mutual_information(all_hi_state, allC_R(used_inds), 0.1);
%                 summ.bits = [summ.bits tmp_bits];

                % reward
                summ.rwd = [summ.rwd mean(allR)];

                % session and mouse info
                summ.mouse_name = [summ.mouse_name {m}];
                summ.task_name = [summ.task_name {state}];
                summ.session_actual = [summ.session_actual sess_ind];

                sess_reset = sess_reset + 1;
                summ.session_reset = [summ.session_reset sess_reset];
                
                summ.allRT = [summ.allRT median(fs_b.(m).(sess).pd.rt)];

                summ.setsize = [summ.setsize 2]; % only two states

                summ.nTrials = [summ.nTrials length(allC_R)];

                % calculate false alarms
                summ.false_alarm = [summ.false_alarm {calculate_false_alarm_trials(fs_b.(m).(sess))}];
            end
        end
    end
end

% do sess_mult
s1 = '0540'; % all reward probability variants
s2 = '0638';
s3 = '1133';
s4 = '2222';
s5 = '3311';
s6 = '3806';
s7 = '4005';

relev_sess = sess_mult; state = 'mult';
for m = all_m
    m = m{:};
    all_sess = fields(fs_b.(m))';
    
    sess_ind = 0; % initialize
    sess_reset = 0;
    for sess = all_sess
        sess = sess{:};
        sess_ind = sess_ind + 1;
        if any(strcmp(relev_sess, sess)) % sess_mult sessions
            % rightward choices in each state
            allC_R_s1 = [];
            allC_R_s2 = [];
            allC_R_s3 = [];
            allC_R_s4 = [];
            allC_R_s5 = [];
            allC_R_s6 = [];
            allC_R_s7 = [];
            % marginal distribution of each state
            all_s1 = [];
            all_s2 = [];
            all_s3 = [];
            all_s4 = [];
            all_s5 = [];
            all_s6 = [];
            all_s7 = [];

            allC_R = fs_b.(m).(sess).pd.allC_R;
            allC_L = fs_b.(m).(sess).pd.allC_L;
            allR_R = fs_b.(m).(sess).pd.allR_R;
            allR_L = fs_b.(m).(sess).pd.allR_L;
            allR = abs(fs_b.(m).(sess).pd.allR);
            bs_all = fs_b.(m).(sess).pd.bs_corrected; bs_all = [bs_all length(allC_R)+1];
            bp = fs_b.(m).(sess).pd.bp;
            cfrac_tmp = []; % for storing single session cfrac values
            rfrac_tmp = [];
            used_inds = []; % used for MI calculation
            for bs_ind = 1:length(bs_all) - 1 % iterate through each block
                inds_of_int = bs_all(bs_ind):bs_all(bs_ind + 1) - 1;

                % remove early trials to capture steady state
                if length(inds_of_int) > ind_offset
                    inds_of_int = inds_of_int(ind_offset:end); % remove trials 1:ind_offset-1

                    tmp_bp = bp{bs_ind};
                    tmp_bp = regexp(tmp_bp, '/', 'split'); 
                    tmp_bp = tmp_bp{1}; % just first block
                    tmp_bp = floor(str2double(tmp_bp));
                    if tmp_bp == 5
                        allC_R_s1 = [allC_R_s1 allC_R(inds_of_int)];
                        all_s1 = [all_s1 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 6
                        allC_R_s2 = [allC_R_s2 allC_R(inds_of_int)];
                        all_s2 = [all_s2 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 11
                        allC_R_s3 = [allC_R_s3 allC_R(inds_of_int)];
                        all_s3 = [all_s3 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 22
                        allC_R_s4 = [allC_R_s4 allC_R(inds_of_int)];
                        all_s4 = [all_s4 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 33
                        allC_R_s5 = [allC_R_s5 allC_R(inds_of_int)];
                        all_s5 = [all_s5 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 38
                        allC_R_s6 = [allC_R_s6 allC_R(inds_of_int)];
                        all_s6 = [all_s6 ones(1, length(inds_of_int))];
                    elseif tmp_bp == 40
                        allC_R_s7 = [allC_R_s7 allC_R(inds_of_int)];
                        all_s7 = [all_s7 ones(1, length(inds_of_int))];
                    else 
                        error('BP not found\n')
                    end
                    used_inds = [used_inds inds_of_int];
                    summ.um_task_block = [summ.um_task_block {state}]; % task for each block

                    % compute soltani-style undermatching (not used)
                    choice_f = sum(allC_R(inds_of_int))/length(allC_R(inds_of_int));
                    reward_f = sum(fs_b.(m).(sess).pd.allR_R(inds_of_int))/ ...
                               sum(allR(inds_of_int));
                    if sign(reward_f - 0.5) ~= 0
                        tmp_undermatching = (choice_f - reward_f)*sign(reward_f - 0.5);
                    else
                        tmp_undermatching = (choice_f - reward_f);
                    end
                    summ.um_soltani = [summ.um_soltani tmp_undermatching];
    
                    % store cfrac and rfrac values, within session, for within session undermatching slope
                    cfrac_tmp = [cfrac_tmp sum(allC_R(inds_of_int)) / sum(allC_L(inds_of_int))];
                    rfrac_tmp = [rfrac_tmp sum(allR_R(inds_of_int)) / sum(allR_L(inds_of_int))];
                end
            end
            % compute generalized matching style undermatching (for regression, block-by-block values)
            summ.cfrac = [summ.cfrac cfrac_tmp];
            summ.rfrac = [summ.rfrac rfrac_tmp];

            % calculate session-wide undermatching slope
            rfrac_log2_tmp = log2(rfrac_tmp);
            rfrac_log2_tmp(isinf(rfrac_log2_tmp)) = NaN;
            cfrac_log2_tmp = log2(cfrac_tmp);
            um_mod = fitlm(rfrac_log2_tmp, cfrac_log2_tmp);
            summ.um_slope_session = [summ.um_slope_session um_mod.Coefficients.Estimate(2)];

            % calculate MI
            all_s1 = length(all_s1)/length(allR(used_inds)); 
            all_s2 = length(all_s2)/length(allR(used_inds));
            all_s3 = length(all_s3)/length(allR(used_inds));
            all_s4 = length(all_s4)/length(allR(used_inds));
            all_s5 = length(all_s5)/length(allR(used_inds));
            all_s6 = length(all_s6)/length(allR(used_inds));
            all_s7 = length(all_s7)/length(allR(used_inds));     

            pi_base = [1-mean(allC_R(used_inds)) mean(allC_R(used_inds))];
            p_state = [all_s1 all_s2 all_s3 all_s4 all_s5 all_s6 all_s7];
            pi = [1-mean(allC_R_s1) mean(allC_R_s1);
                  1-mean(allC_R_s2) mean(allC_R_s2);
                  1-mean(allC_R_s3) mean(allC_R_s3);
                  1-mean(allC_R_s4) mean(allC_R_s4);
                  1-mean(allC_R_s5) mean(allC_R_s5);
                  1-mean(allC_R_s6) mean(allC_R_s6);
                  1-mean(allC_R_s7) mean(allC_R_s7)];
            [tmp_bits, ~, tmp_cond_ent] = p_com(pi, pi_base, p_state);
            summ.bits = [summ.bits tmp_bits];
            summ.cond_ent = [summ.cond_ent tmp_cond_ent];
            summ.setsize = [summ.setsize numel(unique(bp))];
            summ.rwd = [summ.rwd mean(allR)];   

            summ.mouse_name = [summ.mouse_name {m}];
            summ.task_name = [summ.task_name {state}];
            summ.session_actual = [summ.session_actual sess_ind];

            sess_reset = sess_reset + 1;
            summ.session_reset = [summ.session_reset sess_reset];
            
            summ.allRT = [summ.allRT median(fs_b.(m).(sess).pd.rt)];

            summ.nTrials = [summ.nTrials length(allC_R)];

            % calculate false alarms
            summ.false_alarm = [summ.false_alarm {calculate_false_alarm_trials(fs_b.(m).(sess))}];
        end
    end
end
summ.rfrac_log2 = log2(summ.rfrac);
summ.rfrac_log2(isinf(summ.rfrac_log2)) = NaN;
summ.cfrac_log2 = log2(summ.cfrac);

data_table = table(summ.session_actual', summ.session_reset', summ.bits', summ.rwd', summ.cond_ent', ...
                   summ.allRT', summ.mouse_name', summ.task_name', summ.um_slope_session', summ.nTrials', ...
                   summ.false_alarm', ...
    'VariableNames', {'Session_Actual', 'Session_Reset', 'Bits', 'Reward', 'Conditional_Entropy', ...
                      'RT', 'Mouse', 'Task', 'UM_Slope_Session', 'Trials', 'False_Alarm'});
um_table = table(summ.um_soltani', summ.um_task_block', summ.cfrac', summ.rfrac', summ.cfrac_log2', summ.rfrac_log2', ...
    'VariableNames', {'UM_Soltani','UM_Task','CFrac','RFrac','CFrac_Log2','RFrac_Log2'});
save('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\matching_data.mat', 'data_table', 'um_table', ...
    'sess_40_10', 'sess_40_7', 'sess_40_5', 'sess_50_5', 'sess_mult')
fprintf('Data saved\n')
%% load data (start here to reproduce all analyses)

if exist('data_table', 'var') == 0
    load('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\matching_data.mat')
    saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\figures\bari2019';
    myColors = importColors_bb();
    fprintf('Matching data loaded \n')
end
c_4010 = myColors.blue_dull;
c_407 = myColors.darkGray;
c_405 = myColors.blue_bright;

%% complexity-reward frontier for one task variant
f_one = figure('units','normalized','position',[0.1 0.3 0.5 0.5]);
f_prob = subplot(121); hold on
f_comp = subplot(122); hold on

state = '40/10';

subplot(f_prob)
p_right = str2double(extractBefore(state, '/'))/100;
p_left = str2double(extractAfter(state, '/'))/100;
pi = linspace(0,1,1e3);
rhat_right = p_right./(pi + p_right*(1 - pi));
rhat_left = p_left./((1 - pi) + p_left.*pi);
rhat = rhat_right.*pi + rhat_left.*(1 - pi);

plot(pi, rhat_right, ':', 'linewidth', 2, 'Color', myColors.lightGray)
plot(pi, rhat_left, '--', 'Color', myColors.darkGray)
plot(pi, rhat, 'Linewidth', 1, 'Color', c_4010)
[max_rhat, max_rhat_ind] = max(rhat);
plot(pi(max_rhat_ind), max_rhat, 'o', 'linewidth', 2, 'Color', c_4010)
xlim([0.5 1])
ylabel('Average reward')
xlabel('P(rightward choice)')
set(gca,'tickdir','out')
legend({'$\hat{r}_{\textrm{right}}$','$\hat{r}_{\textrm{left}}$','$\hat{r}$','Matching'},'interpreter','latex')

% policy complexity plot
subplot(f_comp)
if strcmp(state, '40/10')
    r_i = 0.4;
    r_j = 0.1;
    c_of_int = c_4010;
elseif strcmp(state, '40/7')
    r_i = 0.4;
    r_j = 0.07;
    c_of_int = c_407;
elseif strcmp(state, '40/5')
    r_i = 0.4;
    r_j = 0.05;
    c_of_int = c_405;
else
    error('State not found\n')
end

%policy complexity
pi_base = [0.5 0.5];
p_state = [0.5 0.5];

total_n_bits = [];
total_r = [];
total_probs = linspace(0.5, 0, 1e2);
for prob_ind = total_probs
    pi = [prob_ind 1-prob_ind;
          1-prob_ind prob_ind];
    n_bits = p_com(pi, pi_base, p_state);

    total_n_bits = [total_n_bits n_bits];
    total_r = [total_r calc_rhat_matching(1-prob_ind, r_i, r_j)];
end
[max_r, max_bit_ind] = max(total_r);
switch state
    case '40/10'
        %match_bits_4010 = total_n_bits(max_bit_ind); % match_bits_4010 = 0.4121
        match_bits_4010 = 0.4121;
    case '40/5'
        %match_bits_405 = total_n_bits(max_bit_ind); % 0.6314
        match_bits_405 = 0.6314;
end

t(1) = plot(total_n_bits(max_bit_ind), max_r, 'o', 'linewidth', 2, 'MarkerEdgeColor', c_of_int);
% curve if MI is fixed
t(2) = plot(total_n_bits, total_r, '--', 'Color', c_of_int, 'linewidth', 1);
% curve if MI is not fixed
t(3) = plot(total_n_bits(1:max_bit_ind), total_r(1:max_bit_ind), 'Color', c_of_int, 'linewidth', 1);
plot(total_n_bits(max_bit_ind + 1:end), total_r(max_bit_ind)*ones(1, length(total_r) - max_bit_ind), ...
    'Color', c_of_int, 'linewidth', 1)
xlabel('Policy complexity (bits)')
ylabel('Average reward')
ylim([0.3 0.5])
legend(t,{'Matching','Fixed Constraint','Reward-Complexity Frontier'},'location','southeast')
set(gca,'tickdir','out')

saveFigureIteration(f_one, saveFigLoc, 'rwdComplexity_ideal')

%% ideal policies for each bit rate
f_policy = figure('units','normalized','position',[0 0 0.4 0.3]);
sp(1) = subplot(131); hold on
bar([0.5 0.5])
sp(2) = subplot(132); hold on
bar([total_probs(max_bit_ind) 1-total_probs(max_bit_ind)])
sp(3) = subplot(133); hold on
bar([0 1])
set(sp, 'tickdir', 'out', 'xtick', 1:2, 'xticklabel', {'Left','Right'}, 'ytick', [0:0.25:1])
ylim(sp, [0 1])
ylabel(sp, '$\pi(a | p_{L} = 0.1, p_{R} = 0.4)$','interpreter','latex')
xlabel(sp, 'Actions')

saveFigureIteration(f_policy, saveFigLoc, 'rwdComplexity_ideal_policy')


%% generate complexity-reward frontier w/ empiric data
f_all = figure('units','normalized','position',[0.1 0.1 0.5 0.4]); 
f_all_polComp = subplot(121); hold on
f_all_polComp_mean = subplot(122); hold on
all_states = {'40/10','40/5'};
clear legend_polcomp_emp
for state = all_states
    state = state{:};

    if strcmp(state, '40/10')
        r_i = 0.4;
        r_j = 0.1;
        c_of_int = c_4010;
        relev_sess = sess_40_10;
    elseif strcmp(state, '40/7')
        r_i = 0.4;
        r_j = 0.07;
        c_of_int = c_407;
        relev_sess = sess_40_7;
    elseif strcmp(state, '40/5')
        r_i = 0.4;
        r_j = 0.05;
        c_of_int = c_405;
        relev_sess = sess_40_5;
        eb_ind = 2;
    elseif strcmp(state, '50/5')
        r_i = 0.5;
        r_j = 0.05;
        c_of_int = myColors.gray;
        relev_sess = c_505;
    else
        error('State not found\n')
    end

    %policy complexity
    pi_base = [0.5 0.5];
    p_state = [0.5 0.5];

    total_n_bits = [];
    total_r = [];
    total_probs = linspace(0.5, 0, 1e2);
    for prob_ind = total_probs
        pi = [prob_ind 1-prob_ind;
              1-prob_ind prob_ind];
        n_bits = p_com(pi, pi_base, p_state);

        total_n_bits = [total_n_bits n_bits];
        total_r = [total_r calc_rhat_matching(1-prob_ind, r_i, r_j)];
    end
    [max_r, max_bit_ind] = max(total_r);
    max_prob = 1-total_probs(max_bit_ind);
    max_bits = total_n_bits(max_bit_ind);
    switch state
        case '40/10'
            match_bits_4010 = max_bits;
        case '40/7'
            match_bits_407 = max_bits;
        case '40/5'
            match_bits_405 = max_bits;
    end

    subplot(f_all_polComp)
    plot(total_n_bits(max_bit_ind), max_r, 'o', 'MarkerEdgeColor', c_of_int, 'linewidth', 2);
    tmp_legend_polcomp_emp = plot(total_n_bits(1:max_bit_ind), total_r(1:max_bit_ind), 'Color', c_of_int, 'linewidth', 1);
    plot(total_n_bits(max_bit_ind:end), total_r(max_bit_ind)*ones(1, length(total_r) - max_bit_ind + 1), ...
    'Color', c_of_int, 'linewidth', 1)
    if exist('legend_polcomp_emp','var') == 0
        legend_polcomp_emp = tmp_legend_polcomp_emp;
    else
        legend_polcomp_emp(end + 1) = tmp_legend_polcomp_emp;
    end
    
    for m = unique(data_table.Mouse)'
        m = m{:};
        m_mask = contains(data_table.Mouse, m);
        task_mask = contains(data_table.Task, state);
        scatter(mean(data_table(m_mask & task_mask, :).Bits), ...
            mean(data_table(m_mask & task_mask, :).Reward), 25, 'Filled', 'MarkerFaceColor', c_of_int);
    end

    subplot(f_all_polComp_mean)
    [~, eb_ind] = intersect(all_states, state);
    bit_data = data_table(task_mask, :).Bits;
    bit_data_median = median(bit_data);
    bit_data_ci = prctile(bootstrp(1e3, @median, bit_data), [2.5 97.5]);
    errorbar(eb_ind, bit_data_median, bit_data_median-bit_data_ci(1), bit_data_median-bit_data_ci(2), 'o', ...
             'linewidth', 2, 'Color', c_of_int', 'MarkerFaceColor', c_of_int)
end
set([f_all_polComp f_all_polComp_mean], 'tickdir', 'out')
xlabel(f_all_polComp, 'Policy complexity (bits)')
ylabel(f_all_polComp, 'Average reward')
legend(legend_polcomp_emp, all_states)

xlabel(f_all_polComp_mean, 'Group')
ylabel(f_all_polComp_mean, 'Policy complexity (bits)')
xlim(f_all_polComp_mean, [0.5 length(all_states)+0.5])
ylim(f_all_polComp_mean, [0 0.2])
set(f_all_polComp_mean,'xtick',1:length(all_states),'xticklabel',all_states)

saveFigureIteration(f_all, saveFigLoc, 'rwdComplexity_empiric')

% stats
state_mask = contains(data_table.Task, '40/10');
b4010 = data_table(state_mask, :).Bits;
state_mask = contains(data_table.Task, '40/5');
b405 = data_table(state_mask, :).Bits;
[~,p4010_norm] = lillietest(b4010);
[~,p405_norm] = lillietest(b405);
fprintf('pval 4010 normality: %0.4f\n', p4010_norm);
fprintf('pval 405 normality: %0.4f\n', p405_norm);
[p,h,stats] = ranksum(b4010, b405)
    

%% cfrac vs rfrac data
% simulate fixed policy
f_um = figure('units','normalized','position',[0.1 0.1 0.37 0.7]);
scale_var = 5; % stretch the lines out
pi_fixed = 0.75;
cfrac_log2 = log2(pi_fixed/(1 - pi_fixed));
all_states = {'40/10','40/5'};
h_sim_um_slope = subplot(221); hold on; title('Simulation')
h_sim_um_coef = subplot(222); hold on
clear legend_sim
for state = all_states
    state = state{:};
    switch state
        case '40/10'
            p_set = [0.4 0.1];
            c_of_int = c_4010;
            match_bits = 0.4121;
        case '40/5'
            p_set = [0.4 0.05];
            c_of_int = c_405;
            match_bits = 0.6314;
    end
    rfrac = (p_set(1)/(pi_fixed + p_set(1)*(1 - pi_fixed))*pi_fixed) / ...
            (p_set(2)/((1 - pi_fixed) + p_set(2)*pi_fixed)*(1 - pi_fixed));
    rfrac_log2 = log2(rfrac);

    subplot(h_sim_um_slope)
    plot([-5 5],[-5 5],':','Color', myColors.darkGray, 'linewidth', 1.3)
    plot(scale_var*[-rfrac_log2 rfrac_log2], scale_var*[-cfrac_log2 cfrac_log2], 'linewidth', 2, 'Color', c_of_int);

    subplot(h_sim_um_coef)
    tmp_legend = plot(match_bits, cfrac_log2/rfrac_log2, 'o', 'Color', c_of_int, 'MarkerFaceColor', c_of_int);
    if exist('legend_sim','var') == 0
        legend_sim = tmp_legend;
    else
        legend_sim(end + 1) = tmp_legend;
    end
end
set([h_sim_um_slope h_sim_um_coef],'tickdir','out')
xlim(h_sim_um_slope, [-5 5]); ylim(h_sim_um_slope, [-5 5])
axis(h_sim_um_slope, 'square')
xlabel(h_sim_um_slope, '$\textrm{log}_2(\frac{\textrm{R}_{\textrm{right}}}{\textrm{R}_{\textrm{left}}})$','Interpreter','Latex')
ylabel(h_sim_um_slope, '$\textrm{log}_2(\frac{\textrm{C}_{\textrm{right}}}{\textrm{C}_{\textrm{left}}})$','Interpreter','Latex')

legend(legend_sim, all_states, 'location', 'northeast')
ylim(h_sim_um_coef, [0.5 0.8])
xlim(h_sim_um_coef, [0.2 0.8])
xlabel(h_sim_um_coef, 'Matching policy (Bits)')
ylabel(h_sim_um_coef, sprintf('Matching\n(slope of reward x choice)'))
axis(h_sim_um_coef, 'square')

% behavioral data
all_states = {'40/10','40/5'};
h_um_coef = subplot(224); hold on
h_um_sp = subplot(223); hold on; title('Behavior')
clear legend_coef
for state = all_states
    state = state{:};
    switch state
        case '40/10'
            c_of_int = c_4010;
            match_bits = 0.4121; % hard-coded, can get this from above cell
        case '40/5'
            c_of_int = c_405;
            match_bits = 0.6314; % hard-coded, can get this from above cell
    end
    state_mask = contains(um_table.UM_Task, state);
    plot([-5 5],[-5 5],':','Color', myColors.darkGray, 'linewidth', 1.3)

    subplot(h_um_coef)
    m_um = fitlm(um_table(state_mask, :).RFrac_Log2, um_table(state_mask, :).CFrac_Log2);
    mu_coef = m_um.Coefficients.Estimate(2);
    ci_coef = coefCI(m_um); ci_coef = ci_coef(2, :);
    tmp_legend_coef = errorbar(match_bits, mu_coef , mu_coef-ci_coef(1), mu_coef-ci_coef(2), ...
            'o', 'linewidth', 2, 'Color', c_of_int, 'MarkerFaceColor', c_of_int);
    if exist('legend_coef','var') == 0
        legend_coef = tmp_legend_coef;
    else
        legend_coef(end + 1) = tmp_legend_coef;
    end

    x_pred = [-5 5]';
    y_pred = m_um.predict(x_pred);

    subplot(h_um_sp)
    plot(x_pred, y_pred, 'Color', c_of_int, 'Linewidth', 2)
end
set([h_um_sp h_um_coef], 'tickdir', 'out')
xlim(h_um_sp, [-5 5])
ylim(h_um_sp, [-5 5])
axis(h_um_sp, 'square')
xlabel(h_um_sp, '$\textrm{log}_2(\frac{\textrm{R}_{\textrm{right}}}{\textrm{R}_{\textrm{left}}})$','Interpreter','Latex')
ylabel(h_um_sp, '$\textrm{log}_2(\frac{\textrm{C}_{\textrm{right}}}{\textrm{C}_{\textrm{left}}})$','Interpreter','Latex')

ylim(h_um_coef, [0.5 0.8])
xlim(h_um_coef, [0.2 0.8])
legend(legend_coef, all_states)
xlabel(h_um_coef, 'Matching policy (Bits)')
ylabel(h_um_coef, sprintf('Matching\n(slope of reward x choice)'))
axis(h_um_coef, 'square')

set([h_sim_um_coef h_um_coef], 'ytick', 0.5:0.1:0.8)

saveFigureIteration(f_um, saveFigLoc, 'undermatching_slopes')