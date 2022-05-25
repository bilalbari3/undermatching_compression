clear; close all; clc
saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\figures\williams1985';
myColors = importColors_bb();
% 1: 90, 0.5
% 2: 90, 0.15
% 3: 90, 0.08
% 4: 30, 0.15
% 5: 90, 0.15
% 6: 270, 0.15

p_vi_set = 6./[90 90 90 30 90 270];
p_vr_set = [0.5 0.15 0.08 0.15 0.15 0.15];
animal_odds = [8.18 1.27 0.49 0.23 0.94 3.04]; % expressed as odds ratio
animal_prob = animal_odds./(animal_odds + 1); % convert to prob

figure
p_choice_vr = linspace(0,1,1e3);
for ind = 1:length(p_vi_set)
    if ind == 1
        f(1) = subplot(121); hold on
        ylim([0 0.7])
    elseif ind == 4
        f(2) = subplot(122); hold on
        ylim([0 0.35])
    end
    p_vi = p_vi_set(ind); % current vi probablity
    p_vr = p_vr_set(ind); % current vr probability
    

    r_vi = p_vi ./ ((1 - p_choice_vr) + p_vi.*p_choice_vr); % mean reward from vi option
    r_vr = p_vr; % mean reward from vr option
    r_tot = r_vi.*(1 - p_choice_vr) + r_vr.*p_choice_vr; % total reward from both options
    
    plot(p_choice_vr, r_tot, 'k', 'linewidth', 1); % plot r_tot line

    % get matching index
    [~, match_ind] = min(abs(r_vi - r_vr)); 
    t(1) = plot(p_choice_vr(match_ind), r_tot(match_ind), 'ko', 'Linewidth', 2);
    % get maximizing index
    [~, max_ind] = max(r_tot);
    t(2) = plot(p_choice_vr(max_ind), r_tot(max_ind), 'kx', 'MarkerSize', 10, 'Linewidth', 2, 'Color', [0 0 0]);
    % get animal index
    ap = animal_prob(ind);
    t(3) = plot(ap, (p_vi/((1-ap)+p_vi*ap)).*(1-ap) + r_vr*ap, 'o', 'MarkerSize', 8, 'LineWidth', 2, ...
                'Color', 'none', 'MarkerFaceColor', myColors.vermillion);
end
legend(t, {'Matching','Optimal','Behavior'})
set(f,'tickdir','out')
for cp = f
    subplot(cp)
    xlabel('P(choice to variable-ratio option)')
    ylabel('Average reward')
end

saveFigureIteration(gcf, saveFigLoc, 'williams1985');