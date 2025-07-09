function dydt = Spanish_ticks_INF(t, y)
dydt = zeros(33,1);

%% Tick and Host Densities

H_s = y(25) + y(26) + y(27);            % Total Small Mammal Density
H_u = y(28) + y(29) + y(30);            % Total Ungulate Density

L_q = y(1) + y(13);         
N_q = y(4) + y(16);        
A_q = y(7) + y(19); 

I_TS = y(14) + y(17);                   % Total infected ticks feeding on small mammals
I_TU = y(15) + y(18) + y(20);           % Total infected ticks feeding on ungulates

%% Demographic parameters

[~, ~, ~, r_T, ~, ~, ~, ~, ~, ~, ~, ~, ~] = tick_demographic_parameters;

[a_S, a_U] = host_demographic_parameters;

%% Infection Parameters

[rho, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, c_S, c_U, gamma_S, gamma_U] = infection_parameters;

%% Seasonality of terms in Equations

terms = seasonality_of_terms(t, H_s, H_u, L_q, N_q, A_q);

infect_terms = seasonality_of_infection(t, y);

%% Susceptible Tick Equations

% Larvae
dydt(1) =   terms(1)*y(11) + terms(2)*y(12)  - (terms(3) + terms(4))*y(1)    - terms(5)*y(1);
dydt(2) =   infect_terms(1)*y(1)             - terms(6)*y(2)                 - c_S*y(2)*I_TS;                        
dydt(3) =   infect_terms(3)*y(1)             - terms(6)*y(3)                 - c_U*y(3)*I_TU;                         

% Nymphs
dydt(4) =   terms(6)*(y(2)+y(3))              - (terms(7) + terms(8))*y(4)    - terms(9)*y(4);     
dydt(5) =   infect_terms(5)*y(4)              - terms(10)*y(5)                - c_S*y(5)*I_TS;                                     
dydt(6) =   infect_terms(7)*y(4)              - terms(10)*y(6)                - c_U*y(6)*I_TU;                      

% Adults
dydt(7) =   terms(10)*(y(5) + y(6))           - terms(11)*y(7)                - terms(12)*y(7);    
dydt(8) =   infect_terms(9)*y(7)              - terms(13)*y(8)                  - c_U*y(8)*I_TU;            

% Engorged Females
dydt(9) =   terms(13)*y(8)                    - terms(14)*y(9)                - terms(12)*y(9);  
dydt(10) =                                    - terms(15)*y(10)               - terms(12)*y(10);

% Eggs
dydt(11) =  r_T*terms(14)*(y(9) + (1-rho)*y(21))           - terms(1)*y(11);
dydt(12) =  r_T*terms(15)*(y(10) + (1-rho)*y(22))          - terms(2)*y(12);

%% Infected Tick Equations

% Larvae
dydt(13) =  terms(1)*y(23) + terms(2)*y(24)         - (terms(3) + terms(4))*y(13)   - terms(5)*y(13);
dydt(14) =  terms(3)*y(13) + infect_terms(2)*y(1)   - terms(6)*y(14)                + c_S*y(2)*I_TS;  
dydt(15) =  terms(4)*y(13) + infect_terms(4)*y(1)   - terms(6)*y(15)                + c_U*y(3)*I_TU;

% Nymphs
dydt(16) = terms(6)*(y(14)+y(15))                   - (terms(7) + terms(8))*y(16)   - terms(9)*y(16);
dydt(17) = terms(7)*y(16) + infect_terms(6)*y(4)    - terms(10)*y(17)               + c_S*y(5)*I_TS;
dydt(18) = terms(8)*y(16) + infect_terms(8)*y(4)    - terms(10)*y(18)               + c_U*y(6)*I_TU;

% Adults
dydt(19) = terms(10)*(y(17) + y(18))                - terms(11)*y(19)               - terms(12)*y(19);
dydt(20) = terms(11)*y(19) + infect_terms(10)*y(7)  - terms(13)*y(20)                 + c_U*y(8)*I_TU;

% Engorged Females
dydt(21) = terms(13)*y(20)                - terms(14)*y(21)               - terms(12)*y(21);  
dydt(22) =                              - terms(15)*y(22)               - terms(12)*y(22);

% Eggs
dydt(23) = r_T*terms(14)*y(21)*rho          - terms(1)*y(23);
dydt(24) = r_T*terms(15)*y(22)*rho          - terms(2)*y(24);


%% Host Equations

dydt(25) = a_S*H_s    - (infect_terms(11) + infect_terms(13))*y(25)   - a_S*y(25);
dydt(26) =              (infect_terms(11) + infect_terms(13))*y(25)   - gamma_S*y(26)     - a_S*y(26);
dydt(27) = gamma_S*y(26)        - a_S*y(27);

dydt(28) = a_U*H_u    - (infect_terms(12) + infect_terms(14) + infect_terms(15))*y(28)   - a_U*y(28);
dydt(29) =              (infect_terms(12) + infect_terms(14) + infect_terms(15))*y(28)   - gamma_U*y(29)     - a_U*y(29);
dydt(30) = gamma_U*y(29)        - a_U*y(30);

%% Cumulative Densities

% Cumulative infected questing tick densities
dydt(31) = (terms(1)*y(23) + terms(2)*y(24));       % larvae
dydt(32) = terms(6)*(y(14)+y(15));                  % nymph
dydt(33) = terms(10)*(y(17) + y(18));               % adult 

% Cumulative total questing tick densities
% dydt(31) = terms(1)*(y(11)+y(23)) + terms(2)*(y(12)+y(24));  
% dydt(32) = terms(6)*(y(2)+y(3)+y(14)+y(15));  
% dydt(33) = terms(10)*(y(5)+y(6)+y(17)+y(18)); 
end