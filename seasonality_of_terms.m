function terms = seasonality_of_terms(t, H_s, H_u, L_q, N_q, A_q)
terms = zeros(1, 15);                    

[preovi, ovi, nu, ~, m_L, m_N, sigma_L, sigma_N, sigma_A, b_L, b_N, b_A, b_oA] = tick_demographic_parameters;
[beta_1, beta_2, beta_3, beta_4, beta_5, s_1, s_2, s_3, s_4, s_5] = attach_parameters(H_s, H_u);

%% Timings 
timings = [90, 120, 135, 151, 196, 237, 258, 273, 304];

% t = 90,   Start of Apr    START:  Adults are active, proportion of females lay eggs, other females start preoviposition 
% t = 120,  Start of May    START:  Initial egg hatching, initial larvae are active, initial nymphs are active
% t = 135,  Mid-May         START:  Oviposition of other females
% t = 151,  Start of Jun    END OF: initial egg hatching and larvae
% t = 196,  Mid-Jul         START:  Other egg hatching, larvae active  
%                           END OF: initial nymphs 
% t = 237,  Late-Aug        START:  Nymphs active  
%                           END OF: Spring adults
% t = 258,  Mid-Sept            START:  Adults are active, Engorged females are active
%                           END OF: Engorged females from spring adults
% t = 273,  Start of Oct    END OF: Egg laying, Egg hatching, larvae, nymphs
% t = 304,  End of Oct      END OF: Adults, engorged females (everything goes into over winter))

% terms(1) - hatching of eggs laid by late females
% terms(2) - hatching of eggs laid by females that manage to complete
%               preoviposition prior to overwintering

% terms(3) - attachment of larvae to small mammals
% terms(4) - attachment of larvae to ungulates
% terms(5) - death rate of questing larvae
% terms(6) - detachment and moulting rate of larvae

% terms(7) - attachment of nymphs to small mammals
% terms(8) - attachment of nymphs to ungulates
% terms(9) - death rate of questing nymphs
% terms(10) - detachment and moulting rate of nymphs

% terms(11) - attachment of adults
% terms(12) - death rate of adults and engorged females
% terms(13) - preoviposition and oviposition of females that did not
%               complete oviposition prior to overwintering
% terms(14) - rate of oviposition for females that completed preoviposition
%               before overwintering.

%%

% --------NOTE ------- Density-dependent attachment terms need to be 
% dependent on when ticks are actually active.

if rem(t, 365)>=timings(1) && rem(t, 365)<timings(2)        % APR            
    att_U = 1/(1 + s_5*A_q);

    terms(1) = 0;
    terms(2) = 0;

    terms(3) = 0;
    terms(4) = 0;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = 0;
    terms(8) = 0;
    terms(9) = b_N;
    terms(10) = m_N*sigma_N/(m_N + sigma_N);

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = 0;
    terms(15) = ovi;

elseif rem(t, 365)>=timings(2) && rem(t, 365)<timings(3)    % MAY
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);   

    terms(1) = 0;
    terms(2) = nu;

    terms(3) = beta_1*att_S*H_s;
    terms(4) = beta_2*att_U*H_u;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = beta_3*att_S*H_s;
    terms(8) = beta_4*att_U*H_u;
    terms(9) = b_N;
    terms(10) = 0;

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = 0;
    terms(15) = ovi;

elseif rem(t, 365)>=timings(3) && rem(t, 365)<timings(4)    % MID MAY
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);   

    terms(1) = 0;
    terms(2) = nu;

    terms(3) = beta_1*att_S*H_s;
    terms(4) = beta_2*att_U*H_u;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = beta_3*att_S*H_s;
    terms(8) = beta_4*att_U*H_u;
    terms(9) = b_N;
    terms(10) = 0;

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = preovi*ovi/(preovi + ovi);
    terms(15) = ovi;

elseif rem(t, 365)>=timings(4) && rem(t, 365)<timings(5)    % JUNE
    att_S = 1/(1 + s_3*N_q);             
    att_U = 1/(1 + s_4*N_q + s_5*A_q);   

    terms(1) = 0;
    terms(2) = nu;

    terms(3) = 0;   %let's turn this on/off here and see what happens
    terms(4) = 0;   %let's turn this on/off here and see what happens
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = beta_3*att_S*H_s;
    terms(8) = beta_4*att_U*H_u;
    terms(9) = b_N;
    terms(10) = 0;

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = preovi*ovi/(preovi + ovi);
    terms(15) = 0;

elseif rem(t, 365)>=timings(5) && rem(t, 365)<timings(6)    % JULY
    att_S = 1/(1 + s_1*L_q);             
    att_U = 1/(1 + s_2*L_q);   

    terms(1) = nu;
    terms(2) = nu;

    terms(3) = beta_1*att_S*H_s;
    terms(4) = beta_2*att_U*H_u;
    terms(5) = b_L;
    terms(6) = 0;

    terms(7) = 0;
    terms(8) = 0;
    terms(9) = b_N;
    terms(10) = 0;

    terms(11) = 0;
    terms(12) = b_oA;

    terms(13) = sigma_A;
    terms(14) = preovi*ovi/(preovi + ovi);
    terms(15) = 0;

elseif rem(t, 365)>=timings(6) && rem(t, 365)<timings(7)    % AUG
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q);   

    terms(1) = nu;
    terms(2) = nu;

    terms(3) = beta_1*att_S*H_s;
    terms(4) = beta_2*att_U*H_u;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = beta_3*att_S*H_s;
    terms(8) = beta_4*att_U*H_u;
    terms(9) = b_N;
    terms(10) = 0;

    terms(11) = 0;
    terms(12) = b_oA;

    terms(13) = sigma_A;
    terms(14) = preovi*ovi/(preovi + ovi);
    terms(15) = 0;
    
elseif rem(t, 365)>=timings(7) && rem(t, 365)<timings(8)    % SEPT
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);   

    terms(1) = nu;
    terms(2) = nu;

    terms(3) = beta_1*att_S*H_s;
    terms(4) = beta_2*att_U*H_u;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = beta_3*att_S*H_s;
    terms(8) = beta_4*att_U*H_u;
    terms(9) = b_N;
    terms(10) = m_N*sigma_N/(m_N + sigma_N);

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = preovi*ovi/(preovi + ovi);
    terms(15) = 0;

elseif rem(t, 365)>=timings(8) && rem(t, 365)<timings(9)    % OCT           
    att_U = 1/(1 + s_5*A_q);   

    terms(1) = 0;
    terms(2) = 0;

    terms(3) = 0;
    terms(4) = 0;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = 0;
    terms(8) = 0;
    terms(9) = b_N;
    terms(10) = m_N*sigma_N/(m_N + sigma_N);

    terms(11) = beta_5*att_U*H_u;
    terms(12) = b_A;

    terms(13) = sigma_A;
    terms(14) = 0;
    terms(15) = 0;

else                                                        % NOV-APR          
    % att_U = 1/(1 + s_5*A_q);   

    terms(1) = 0;
    terms(2) = 0;

    terms(3) = 0;
    terms(4) = 0;
    terms(5) = b_L;
    terms(6) = m_L*sigma_L/(m_L + sigma_L);

    terms(7) = 0;
    terms(8) = 0;
    terms(9) = b_N;
    terms(10) = m_N*sigma_N/(m_N + sigma_N);

    terms(11) = 0;%beta_5*att_U*H_u/2;
    terms(12) = b_oA;

    terms(13) = 0;
    terms(14) = 0;
    terms(15) = 0;
end

end