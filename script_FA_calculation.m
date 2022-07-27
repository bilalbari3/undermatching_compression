state_mask_4010 = contains(data_table.Task, '40/10');
state_mask_405 = contains(data_table.Task, '40/5');
bits_4010 = data_table(state_mask_4010, :).Bits;
bits_405 = data_table(state_mask_405, :).Bits;
FA_4010 = data_table(state_mask_4010, :).False_Alarm;
FA_405 = data_table(state_mask_405, :).False_Alarm;

task = [zeros(length(bits_4010), 1); ones(length(bits_405), 1)];
all_bits = [bits_4010; bits_405];
all_FA = [FA_4010; FA_405];
% all_FA = [cellfun(@mean, FA_4010); cellfun(@mean, FA_405)];
% fitlm([bits_4010; bits_405], [cellfun(@mean, FA_4010); cellfun(@mean, FA_405)])

%%
bit_expanded = [];
FA_expanded = [];
for t = 1:length(all_bits)
    bit_expanded = [bit_expanded; all_bits(t)*ones(length(all_FA{t}), 1)];
    FA_expanded = [FA_expanded; all_FA{t}'];
end

my_mod = fitglm(bit_expanded, FA_expanded, 'Distribution', 'binomial', 'Link', 'logit');