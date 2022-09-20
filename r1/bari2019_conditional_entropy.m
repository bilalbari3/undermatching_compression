% run this after policyComplexity_bari2019 line 364
myColors = importColors_bb();
saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\r1\figures';

if exist('data_table', 'var') == 0
    load('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\matching_data.mat')
    fprintf('Matching data loaded \n')
end
%%

mask_405 = contains(data_table.Task, '40/5');
mask_4010 = contains(data_table.Task, '40/10');
mask_mult = contains(data_table.Task, 'mult');

bits_405 = data_table(mask_405, :).Bits;
bits_4010 = data_table(mask_4010, :).Bits;
bits_mult = data_table(mask_mult, :).Bits;

condent_405 = data_table(mask_405, :).Conditional_Entropy;
condent_4010 = data_table(mask_4010, :).Conditional_Entropy;
condent_2state = [condent_405; condent_4010];
condent_7state = data_table(mask_mult, :).Conditional_Entropy;

ci_2state = bootci(1e3, @median, condent_2state);
ci_7state = bootci(1e3, @median, condent_7state);

%%
f_condEnt = figure('units','normalized','position',[0.1 0.1 0.17 0.4]);
hold on
median_2state = median(condent_2state);
errorbar(1, median_2state, median_2state - ci_2state(1), median_2state - ci_2state(2), 'o', ...
    'linewidth', 2, 'Color', myColors.blue, 'MarkerFaceColor', myColors.blue)

median_7state = median(condent_7state);
errorbar(2, median_7state, median_7state - ci_7state(1), median_7state - ci_7state(2), 'o', ...
    'linewidth', 2, 'Color', myColors.vermillion, 'MarkerFaceColor', myColors.vermillion)

xlim([0.5 2.5])
ylim([0.7 1.0])
xlabel('Task')
ylabel('H(A|S)')
set(gca,'xtick',[1 2],'xticklabel',{'2 states', '7 states'},'tickdir','out')

saveFigureIteration(f_condEnt, saveFigLoc, 'condEnt')
%%
fprintf('2 state (median [95%% CI]): %0.3f [%0.3f - %0.3f]\n', median_2state, ci_2state(1), ci_2state(2))
fprintf('7 state (median [95%% CI]): %0.3f [%0.3f - %0.3f]\n', median_7state, ci_7state(1), ci_7state(2))