function [beta_1, beta_2, beta_3, beta_4, beta_5, s_1, s_2, s_3, s_4, s_5] = attach_parameters(H_s, H_u)

sigma_L = 1/6;      % deatachment rate for larvae
sigma_N = 1/11;     % deatachment rate for nymph 
sigma_A = 1/17;     % deatachment rate for adult

beta_1 = 1/(20*H_s);    % attachment coefficient for larvae on small mammals
beta_2 = 1/(20*H_u);    % attachment coefficient for larvae on ungulates     
beta_3 = 1/(15*H_s);    % attachment coefficient for nymph on small mammals
beta_4 = 1/(15*H_u);    % attachment coefficient for nymph on ungulates     
beta_5 = 1/(20*H_u);    % attachment coefficient for adults on ungulates 

n_max_s = 900;      % max ticks on a small mammal host
n_max_u = 300;      % max ticks on an ungulate host

% Strength of density-dependence in attachment
s_1 = beta_1/(sigma_L*n_max_s);     % larvae on small mammals
s_2 = beta_2/(sigma_L*n_max_u);     % larvae on ungulates
s_3 = beta_3/(sigma_N*n_max_s);     % nymph on small mammals
s_4 = beta_4/(sigma_N*n_max_u);     % nymph on ungulates
s_5 = beta_5/(sigma_A*n_max_u);     % adults on ungulates
end