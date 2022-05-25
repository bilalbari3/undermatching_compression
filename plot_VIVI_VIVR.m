clear; close all; clc
saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\figures\general';
myColor = importColors_bb();

ca = myColor.gray;
cb = myColor.gray;
cr = myColor.gray;

pi = [0.05 0.1 0.4 0.7];
t = (0:10)';
Pi =  1 - (1 - pi).^(t + 1);

h_VIVI_VIVR = figure('units','normalized','position',[0 0.25 0.8 0.5]);
h_baiting = subplot(131); hold on
plot(t, Pi, 'k', 'linewidth', 2, 'Color', myColor.gray)
set(h_baiting,'xtick', 0:10)
xlabel(sprintf('Number of consecutive\nchoices of opposite arm'))
ylabel('P(reward)')

h_VIVI = subplot(132); hold on
clear t
pa = 0.4;
pb = 0.1;
pi = linspace(0,1,1e3);
ra = pa./(pi + pa*(1 - pi));
rb = pb./((1 - pi) + pb*pi);
rhat = ra.*pi + rb.*(1 - pi);
t(1) = plot(pi, ra, '--', 'Color', ca, 'linewidth', 2);
t(2) = plot(pi, rb, ':', 'Color', cb, 'linewidth', 2);
t(3) = plot(pi, rhat, 'Color', cr, 'linewidth', 3);
[rhat_max, rhat_max_i] = max(rhat);
t(4) = plot(pi(rhat_max_i), rhat_max, 'ko', 'linewidth', 2);
legend(t,{'$r_a (p_a = 0.4)$', ...
          '$r_b (p_b = 0.1)$', ...
          '$V^{\pi}$', ...
          'Matching / Optimal'},'interpreter','latex', ...
          'location', 'best')

h_VIVR = subplot(133); hold on
clear t
pa = 0.3;
pb = 0.24;
pi = linspace(0,1,1e3);
ra = pa.*ones(1, length(pi));
rb = pb./((1 - pi) + pb*pi);
rhat = ra.*pi + rb.*(1 - pi);
t(1) = plot(pi, ra, '--', 'Color', ca, 'linewidth', 2);
t(2) = plot(pi, rb, ':', 'Color', cb, 'linewidth', 2);
t(3) = plot(pi, rhat, 'Color', cr, 'linewidth', 3);
[rhat_max, rhat_max_i] = max(rhat);
[~, rhat_match_i] = min(abs(rb-pa));
t(4) = plot(pi(rhat_match_i), rhat(rhat_match_i), 'ko', 'Linewidth', 2);
t(5) = plot(pi(rhat_max_i), rhat_max, 'kx', 'MarkerSize', 10, 'linewidth', 2);
legend(t,{'$r_a (p_a = 0.3)$', ...
          '$r_b (p_b = 0.24)$', ...
          '$V^{\pi}$', ...
          'Matching', ...
          'Optimal'},'interpreter','latex', ...
          'location', 'best')

set([h_baiting h_VIVI h_VIVR], 'tickdir','out')
ylim([h_VIVI h_VIVR], [0 1])
xlabel([h_VIVI h_VIVR],'P(choice to option a)')
ylabel([h_VIVI h_VIVR], 'Average reward')

saveFigureIteration(h_VIVI_VIVR, saveFigLoc, 'VIVI_VIVR');