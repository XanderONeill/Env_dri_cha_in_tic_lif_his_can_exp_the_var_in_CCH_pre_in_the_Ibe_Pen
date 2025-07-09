function terms = seasonality_of_infection(t, y)
terms = zeros(1, 15);

H_s = y(25) + y(26) + y(27);            % Total Small Mammal Density
H_u = y(28) + y(29) + y(30);            % Total Ungulate Density

L_q = y(1) + y(13);         N_q = y(4) + y(16);        A_q = y(7) + y(19); 

[beta_1, beta_2, beta_3, beta_4, beta_5, s_1, s_2, s_3, s_4, s_5] = attach_parameters(H_s, H_u);
[~, p_1, p_2, p_3, p_4, p_5, q_1, q_2, q_3, q_4, q_5, ~, ~, ~, ~] = infection_parameters;

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

if rem(t, 365)>=timings(1) && rem(t, 365)<timings(2)            % APR
    att_U = 1/(1 + s_5*A_q);

    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    % Adult to ungulate transmission
    terms(1, 15) = beta_5*att_U*q_5*y(19);

    terms(1, [1:8, 11:14]) = 0;

elseif rem(t, 365)>=timings(2) && rem(t, 365)<timings(3)        % MAY
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);   

    % Small mammal to Larvae transmission
    terms(1, 1) = beta_1*att_S*(y(25) + (1-p_1)*y(26) + y(27));
    terms(1, 2) = beta_1*att_S*p_1*y(26);
    
    % Ungulate host to larvae transmission
    terms(1, 3) = beta_2*att_U*(y(28) + (1-p_2)*y(29) + y(30));
    terms(1, 4) = beta_2*att_U*p_2*y(29);
    
    % Small mammal to nymph transmission
    terms(1, 5) = beta_3*att_S*(y(25) + (1-p_3)*y(26) + y(27));
    terms(1, 6) = beta_3*att_S*p_3*y(26);
    
    % Ungulate to nymph transmission
    terms(1, 7) = beta_4*att_U*(y(28) + (1-p_4)*y(29) + y(30));
    terms(1, 8) = beta_4*att_U*p_4*y(29);
    
    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    
    % Tick to Host transmission
    terms(1, 11) = beta_1*att_S*q_1*y(13); 
    terms(1, 12) = beta_2*att_U*q_2*y(13); 
    terms(1, 13) = beta_3*att_S*q_3*y(16);
    terms(1, 14) = beta_4*att_U*q_4*y(16);
    terms(1, 15) = beta_5*att_U*q_5*y(19);

elseif rem(t, 365)>=timings(3) && rem(t, 365)<timings(4)        % MID-MAY
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);  

    % Small mammal to Larvae transmission
    terms(1, 1) = beta_1*att_S*(y(25) + (1-p_1)*y(26) + y(27));
    terms(1, 2) = beta_1*att_S*p_1*y(26);
    
    % Ungulate host to larvae transmission
    terms(1, 3) = beta_2*att_U*(y(28) + (1-p_2)*y(29) + y(30));
    terms(1, 4) = beta_2*att_U*p_2*y(29);
    
    % Small mammal to nymph transmission
    terms(1, 5) = beta_3*att_S*(y(25) + (1-p_3)*y(26) + y(27));
    terms(1, 6) = beta_3*att_S*p_3*y(26);
    
    % Ungulate to nymph transmission
    terms(1, 7) = beta_4*att_U*(y(28) + (1-p_4)*y(29) + y(30));
    terms(1, 8) = beta_4*att_U*p_4*y(29);
    
    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    
    % Tick to Host transmission
    terms(1, 11) = beta_1*att_S*q_1*y(13); 
    terms(1, 12) = beta_2*att_U*q_2*y(13); 
    terms(1, 13) = beta_3*att_S*q_3*y(16);
    terms(1, 14) = beta_4*att_U*q_4*y(16);
    terms(1, 15) = beta_5*att_U*q_5*y(19);

elseif rem(t, 365)>=timings(4) && rem(t, 365)<timings(5)        %JUN
    att_S = 1/(1 + s_3*N_q);             
    att_U = 1/(1 + s_4*N_q + s_5*A_q);     

    % Small mammal to nymph transmission
    terms(1, 5) = beta_3*att_S*(y(25) + (1-p_3)*y(26) + y(27));
    terms(1, 6) = beta_3*att_S*p_3*y(26);
    
    % Ungulate to nymph transmission
    terms(1, 7) = beta_4*att_U*(y(28) + (1-p_4)*y(29) + y(30));
    terms(1, 8) = beta_4*att_U*p_4*y(29);
    
    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    
    % Nymph to small mammal transmission
    terms(1, 13) = beta_3*att_S*q_3*y(16);
    % Nymph to ungulate transmission
    terms(1, 14) = beta_4*att_U*q_4*y(16);
    % Adult to ungulate transmission
    terms(1, 15) = beta_5*att_U*q_5*y(19);

    terms(1, [1:4, 11, 12]) = 0;

elseif rem(t, 365)>=timings(5) && rem(t, 365)<timings(6)        % JUL
    att_S = 1/(1 + s_1*L_q);             
    att_U = 1/(1 + s_2*L_q);

    % Small mammal to Larvae transmission
    terms(1, 1) = beta_1*att_S*(y(25) + (1-p_1)*y(26) + y(27));
    terms(1, 2) = beta_1*att_S*p_1*y(26);
    
    % Ungulate host to larvae transmission
    terms(1, 3) = beta_2*att_U*(y(28) + (1-p_2)*y(29) + y(30));
    terms(1, 4) = beta_2*att_U*p_2*y(29);
    
    % Larvae to Host transmission
    terms(1, 11) = beta_1*att_S*q_1*y(13); 
    terms(1, 12) = beta_2*att_U*q_2*y(13); 

    terms(1, [5:10, 13:15]) = 0;

 elseif rem(t, 365)>=timings(6) && rem(t, 365)<timings(7)       % AUG
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q);  

    % Small mammal to Larvae transmission
    terms(1, 1) = beta_1*att_S*(y(25) + (1-p_1)*y(26) + y(27));
    terms(1, 2) = beta_1*att_S*p_1*y(26);
    
    % Ungulate host to larvae transmission
    terms(1, 3) = beta_2*att_U*(y(28) + (1-p_2)*y(29) + y(30));
    terms(1, 4) = beta_2*att_U*p_2*y(29);
    
    % Small mammal to nymph transmission
    terms(1, 5) = beta_3*att_S*(y(25) + (1-p_3)*y(26) + y(27));
    terms(1, 6) = beta_3*att_S*p_3*y(26);
    
    % Ungulate to nymph transmission
    terms(1, 7) = beta_4*att_U*(y(28) + (1-p_4)*y(29) + y(30));
    terms(1, 8) = beta_4*att_U*p_4*y(29);

    % Tick to Host transmission
    terms(1, 11) = beta_1*att_S*q_1*y(13); 
    terms(1, 12) = beta_2*att_U*q_2*y(13); 
    terms(1, 13) = beta_3*att_S*q_3*y(16);
    terms(1, 14) = beta_4*att_U*q_4*y(16);

    terms(1, [9, 10, 15]) = 0;

elseif rem(t, 365)>=timings(7) && rem(t, 365)<timings(8)        % SEPT     
    att_S = 1/(1 + s_1*L_q + s_3*N_q);             
    att_U = 1/(1 + s_2*L_q + s_4*N_q + s_5*A_q);

    % Small mammal to Larvae transmission
    terms(1, 1) = beta_1*att_S*(y(25) + (1-p_1)*y(26) + y(27));
    terms(1, 2) = beta_1*att_S*p_1*y(26);
    
    % Ungulate host to larvae transmission
    terms(1, 3) = beta_2*att_U*(y(28) + (1-p_2)*y(29) + y(30));
    terms(1, 4) = beta_2*att_U*p_2*y(29);
    
    % Small mammal to nymph transmission
    terms(1, 5) = beta_3*att_S*(y(25) + (1-p_3)*y(26) + y(27));
    terms(1, 6) = beta_3*att_S*p_3*y(26);
    
    % Ungulate to nymph transmission
    terms(1, 7) = beta_4*att_U*(y(28) + (1-p_4)*y(29) + y(30));
    terms(1, 8) = beta_4*att_U*p_4*y(29);
    
    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    
    % Tick to Host transmission
    terms(1, 11) = beta_1*att_S*q_1*y(13); 
    terms(1, 12) = beta_2*att_U*q_2*y(13); 
    terms(1, 13) = beta_3*att_S*q_3*y(16);
    terms(1, 14) = beta_4*att_U*q_4*y(16);
    terms(1, 15) = beta_5*att_U*q_5*y(19);

elseif rem(t, 365)>=timings(8) && rem(t, 365)<timings(9)        % OCT
    att_U = 1/(1 + s_5*A_q);

    % Ungulate host to adult tick transmission
    terms(1, 9) = beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30));
    terms(1, 10) = beta_5*att_U*p_5*y(29);
    % Adult to ungulate transmission
    terms(1, 15) = beta_5*att_U*q_5*y(19);

    terms(1, [1:8, 11:14]) = 0;
else
    %att_U = 1/(1 + s_5*A_q); 

    % Ungulate host to adult tick transmission
    terms(1, 9) = 0;%beta_5*att_U*(y(28) + (1-p_5)*y(29) + y(30))/2;
    terms(1, 10) = 0;%beta_5*att_U*p_5*y(29)/2;
    % Adult to ungulate transmission
    terms(1, 15) = 0;%beta_5*att_U*q_5*y(19)/2;

    terms(1, [1:8, 11:14]) = 0;
end

end