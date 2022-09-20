mask_PDoff = [data.pop] == 3;
mask_PDon = [data.pop] == 4;

rwd_vec = [];
choice_vec = [];
group = [];
subject = [];
NaN_pad = 50;

data_PDoff = data(mask_PDoff);
for s = 1:length(data_PDoff)
    allC_R = data_PDoff(s).action - 1;
    allC = allC_R; allC(allC == 0) = -1;
    reward = data_PDoff(s).reward;
    reward(allC == -1 & reward == 1) = -1; % recode rewarded left choices as -1

    rwd_vec = [rwd_vec; NaN(NaN_pad, 1); reward];
    choice_vec = [choice_vec; NaN(NaN_pad, 1); allC_R];
    group = [group; NaN(NaN_pad, 1); zeros(length(allC_R), 1)];
end

data_PDon = data(mask_PDon);
for s = 1:length(data_PDon)
    allC_R = data_PDon(s).action - 1;
    allC = allC_R; allC(allC == 0) = -1;
    reward = data_PDon(s).reward;
    reward(allC == -1 & reward == 1) = -1; % recode rewarded left choices as -1

    rwd_vec = [rwd_vec; NaN(NaN_pad, 1); reward];
    choice_vec = [choice_vec; NaN(NaN_pad, 1); allC_R];
    group = [group; NaN(NaN_pad, 1); ones(length(allC_R), 1)];
end

subject = zeros(length(group), length(data_PDoff));
ind_start = 21251;
for s = 1:length(data_PDoff)
    subj_input = [NaN(NaN_pad, 1); ones(800, 1)];
    subject(ind_start:ind_start + length(subj_input) - 1, s) = subj_input;
    ind_start = ind_start + length(subj_input);
end

rwdHx = generateHistoryMatrix(rwd_vec, 5);
choiceHx = generateHistoryMatrix(choice_vec, 1);

mod_rwdHx = fitglm([rwdHx choiceHx rwdHx.*group choiceHx.*group], choice_vec, 'Distribution', 'binomial', 'Link', 'logit');
