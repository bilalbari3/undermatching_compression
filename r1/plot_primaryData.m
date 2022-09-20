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
if exist('data','var') == 0 % not loaded
    load('C:\Users\bilal\Documents\gitRepositories\undermatching_compression\data\rutledge09_data.mat')
    fprintf('Rutledge 2009 data loaded\n')
end
myColors = importColors_bb();
saveFigLoc = 'C:\Users\bilal\Documents\gitRepositories\undermatching_compression\r1\figures';

%% plot the data
% Bari 2019
f_raw = figure('units','normalized','position',[0.1 0.1 0.7 0.6]);
sp(1) = subplot(211); hold on
rawSession = 'mBB041d20161028';
a = rawSession(2:6);
helper_plotRawData(fs_b.(a).(rawSession), sp(1))

% Rutledge 2009
sp(2) = subplot(212); hold on
os = writeRutledgeData(data(6));
helper_plotRawData(os, sp(2))

subplot(sp(1))
xlim([0 800])
title('Example mouse data')
ylabel(sprintf('Fraction choices\nto rightward target'))
subplot(sp(2))
title('Example human data')
ylabel(sprintf('Fraction choices\nto red target'))

saveFigureIteration(f_raw, saveFigLoc, 'rawData')


function helper_plotRawData(os, figureLegend)

myColors = importColors_bb();

subplot(figureLegend)

choiceTickOffset = 0.01;
choiceTickLength = 0.2; % length of choice raster ticks
textOffset = 0.1;
rwd_color = [0 0 0];
noRwd_color = [0 0 0];
for i = 1:length(os.pd.allC)
    if os.pd.allC(i) == 1 % right choice
        if os.pd.allR(i) == 1 % reward
            hTicks = plot([i i], [1 1 + choiceTickLength] + choiceTickOffset, 'Color', rwd_color);
        else
            hTicks = plot([i i], [1 1 + choiceTickLength/2] + choiceTickOffset, 'Color', noRwd_color);
        end
    elseif os.pd.allC(i) == -1 % left choice
        if os.pd.allR(i) == -1 % reward
            hTicks = plot([i i], [0 0 - choiceTickLength] - choiceTickOffset, 'Color', rwd_color);
        else
            hTicks = plot([i i], [0 0 - choiceTickLength/2] - choiceTickOffset, 'Color', noRwd_color);
        end
    end
end

% block overlay
blocks = [os.pd.bs_corrected length(os.pd.allC)];
bp = os.pd.bp;
for i = 1:length(blocks) - 1
    text(mean(blocks(i:i + 1)), 1 + choiceTickLength + textOffset + choiceTickOffset, bp{i}, 'HorizontalAlignment', 'center')
    if i ~= 1
        plot([blocks(i) blocks(i)] - 0.5, [0 1], '--', 'Color', myColors.gray)
    end
    bp_num = cellfun(@str2num, strsplit(bp{i}, '/'));
    prob_match = 1 - (bp_num(1)/sum(bp_num));
    hMatch = plot(blocks(i:i + 1), [prob_match prob_match], 'linewidth', 2, 'Color', myColors.yellow);
end

% smooth choice
smoothChoice = smooth(os.pd.allC_R, 11);
hSmooth = plot(smoothChoice, 'k', 'linewidth', 2);

% clean it up
xlabel('Trial')
ylabel('Fraction rightward choices')
set(figureLegend, 'XTick', [0 100:100:ceil(length(os.pd.allC)/100)*100], 'YTick', 0:0.5:1, 'tickdir', 'out')
xlim([0 length(os.pd.allC)])
ylim([-choiceTickLength-0.1 1+choiceTickLength+0.1])
% predicted choice
legend([hTicks hSmooth hMatch], 'Raw Choice', 'Choice (smoothed)', 'Perfect Matching')

end

function os = writeRutledgeData(os_orig)
    block_change = diff(os_orig.block); % find block changes
    block_start = find([1; block_change])'; % assign first trial of block change to an index
    allC_R = os_orig.action - 1;
    allC = allC_R; allC(allC == 0) = -1;
    reward = os_orig.reward';
    reward(allC' == -1 & reward == 1) = -1; % recode rewarded left choices as -1

    % assign block probabilities
    prob = os_orig.prob;
    prob = round(prob * 100)/100;
    prob(prob == 0.14) = 0.043*100;
    prob(prob == 0.25) = 0.075*100;
    prob(prob == 0.75) = 0.225*100;
    prob(prob == 0.86) = 0.257*100;
    bp_num = prob(block_start, :)'; % subselect probabilities corresponding to block start
    for i = 1:length(bp_num)
        bp{i} = [num2str(bp_num(1,i)) '/' num2str(bp_num(2,i))];
    end

    os.pd.bs_corrected = block_start;
    os.pd.allC = allC';
    os.pd.allC_R = allC_R';
    os.pd.allR = reward;
    os.pd.bp = bp;
end