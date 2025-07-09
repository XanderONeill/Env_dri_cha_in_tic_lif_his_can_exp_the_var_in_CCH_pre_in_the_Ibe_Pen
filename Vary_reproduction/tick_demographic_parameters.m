function [preovi, ovi, nu, r_T, m_L, m_N, sigma_L, sigma_N, sigma_A, b_L, b_N, b_A, b_oA] = tick_demographic_parameters

preovi = 2/55;  % Rate of pre-oviposition 
ovi = 2/41;     % Rate of oviposition
nu = 1/31;      % Rate of hatching
r_T = 2000;     % Number of eggs laid on average per engorged female

m_L = 1/14;     % Moulting rate for larvae
m_N = 2/135;     % Moulting + Resting Period for Nymphs

sigma_L = 1/6;      % Detachment rate for larvae
sigma_N = 1/11;     % Detachment rate for nymphs
sigma_A = 1/17;     % Detachment rate for adults

b_L = -log(0.17)/20;    % Natural death rate of questing larvae
b_N = -log(0.5)/15;    % Natural death rate of questing nymph
b_A = -log(0.4481)/30;     % Natural death rate of questing adults

% b_oL = -log(0.8)/212;   % Overwinter death rate of questing larvae 
% b_oN = -log(0.8)/197;   % Overwinter death rate of questing nymph
b_oA = -log(0.8)/120;   % Overwinter death rate of questing adults

end