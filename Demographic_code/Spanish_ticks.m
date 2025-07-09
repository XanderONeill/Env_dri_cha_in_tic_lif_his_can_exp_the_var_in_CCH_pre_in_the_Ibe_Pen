function dydt = Spanish_ticks(t, y)
dydt = zeros(15,1);

%% Tick and Host Densities
L_q = y(1);
%L_fs = y(2); L_fu = y(3);
N_q = y(4);
%N_fs = y(5); N_fu = y(6);
A_q = y(7); 
%A_f = y(8); 
%E_f = y(9) and y(10); 

%Host densities
H_s = 100;      H_u = 20;

%% Demographic parameters

[~, ~, ~, r_T, ~, ~, ~, ~, ~, ~, ~, ~, ~] = demographic_parameters;
% [beta_1, beta_2, beta_3, beta_4, beta_5, s_1, s_2, s_3, s_4, s_5] = attach_parameters(H_s, H_u);

%% Seasonality of terms in Equations

terms = seasonality_of_terms(t, H_s, H_u, L_q, N_q, A_q);

%% Tick Equations

% Larvae
dydt(1) =  (terms(1)*y(11) + terms(2)*y(12)) - (terms(3) + terms(4))*y(1)    - terms(5)*y(1);
dydt(2) =  terms(3)*y(1)        - terms(6)*y(2);                        
dydt(3) =  terms(4)*y(1)        - terms(6)*y(3);                         

% Nymphs
dydt(4) = terms(6)*(y(2)+y(3))  - (terms(7) + terms(8))*y(4)    - terms(9)*y(4);     
dydt(5) = terms(7)*y(4)         - terms(10)*y(5);                                     
dydt(6) = terms(8)*y(4)         - terms(10)*y(6);                      

% Adults
dydt(7) = terms(10)*(y(5) + y(6))   - terms(11)*y(7)       - terms(12)*y(7);    
dydt(8) = terms(11)*y(7)            - terms(13)*y(8);            

% Engorged Females
dydt(9) = terms(13)*y(8)          - terms(14)*y(9)        - terms(12)*y(9) ;  
dydt(10) =                      - terms(15)*y(10)       - terms(12)*y(10);

% Eggs
dydt(11) = r_T*terms(14)*y(9)   - terms(1)*y(11);
dydt(12) = r_T*terms(15)*y(10)  - terms(2)*y(12);

dydt(13) = (terms(1)*y(11) + terms(2)*y(12));                      % cumulative questing larvae density
dydt(14) = terms(6)*(y(2)+y(3));                       % cumulative questing nymphal density
dydt(15) = terms(10)*(y(5) + y(6));     % cumulative questing adult density
end