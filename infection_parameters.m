function [rho, p1, p2, p3, p4, p5, q1, q2, q3, q4, q5, cS, cU, gS, gU] = infection_parameters

rho = 0.17;         % Vertical transmission

p1 = .1;          % small mammal to larvae transmission 
p2 = .1;          % ungulate to larvae transmission 
p3 = 11/6*p1;      % small mammal to nymph transmission    
p4 = 11/6*p2;      % ungulate to nymph transmission 
p5 = 17/6*p2;      % ungulate to adult transmission 

q1 = .0003;       % larvae to small mammal transmission
q2 = .00042;       % larvae to ungulate transmission
q3 = 11/6*q1;      % nymph to small mammal transmission
q4 = 11/6*q2;      % nymph to ungulate transmission
q5 = 17/6*q2;      % adult to ungulate transmission

cS = 0.00000007*0.85;       % tick to tick transmission (small mammals)     
cU = 0.00000084*0.85;         % tick to tick transmission (ungulates)

% c_S = 0.000000059 - limit of where virus persists in ticks, in absence of
% other transmission routes.
% c_U = 0.00000063 - limit of where virus persist in ticks, in absence of other
% transmission routes.
% c_S, c_U combined need to be less than 80% of the above limits for the
% virus to not persist in the absence of other transmission routes.

gS = 1/14;          % small mammal recovery rate
gU = 1/7;           % ungulate recovery rate

end