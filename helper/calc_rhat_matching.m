function r_hat = calc_rhat_matching(pi_i, r_i, r_j)
%calculate reward from matching
% pi_i: probability of choosing option i
% r_i: base reward probability of option i
% r_j: base reward probability of option j

pi_j = 1 - pi_i;
r_hat = r_i/(pi_i + r_i*(1 - pi_i))*pi_i + ...
        r_j/(pi_j + r_j*(1 - pi_j))*pi_j;