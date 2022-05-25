function r_hat = calc_rhat_matchMax(pi_i, r_i, r_j)
%calculate reward from matching
% pi_i: probability of choosing option i (baited option)
% r_i: base reward probability of option i (baited option)
% r_j: base reward probability of option j (standard option)

pi_j = 1 - pi_i;
r_hat = r_i./(pi_i + r_i.*(1 - pi_i)).*pi_i + ...
        r_j.*pi_j;