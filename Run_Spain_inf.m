%%% Script to run the Spanish_v3 CCHFv infection model

close all
clearvars

%% ODE Solver
options = odeset('Refine', 1, 'NonNegative',1:33, 'RelTol', 1e-07, 'AbsTol', 1e-07);

tspan = 0:1:365;
%            Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec, Jan
month_vec = [1,   32,  60,  91,  121, 152, 182, 213, 244, 274, 305, 335, 366];

% y0 = [100*ones(1,12), zeros(1,12), 100, 0, 0, 20, 0, 0, ones(1,3)];
y0 = [10000*ones(1,24), 80, 20, 0, 10, 10, 0, zeros(1,3)];
for i = 1:50
    [t,y] = ode15s(@(t,y) Spanish_ticks_INF(t,y), tspan, y0, options);
    y0 = y(end,:);
    y0(10) = 0.05*y0(9);
    y0(22) = 0.05*y0(21);
    y0([11,12,23,24]) = 0;
    y0([31, 32, 33])=0;
end

%% Data Processing - Total Questing Densities

L_q = zeros(1,12);
N_q = zeros(1,12);
A_q = zeros(1,12);

for i = 1:12
    % L_q(i) = y(month_vec(i), 1) + y(month_vec(i), 13) + (y(month_vec(i+1),31)-y(month_vec(i),31));
    % N_q(i) = y(month_vec(i), 4) + y(month_vec(i), 16) + (y(month_vec(i+1),32)-y(month_vec(i),32));
    % A_q(i) = y(month_vec(i), 7) + y(month_vec(i), 19) + (y(month_vec(i+1),33)-y(month_vec(i),33));
    L_q(i) = y(month_vec(i), 13) + (y(month_vec(i+1),31)-y(month_vec(i),31));
    N_q(i) = y(month_vec(i), 16) + (y(month_vec(i+1),32)-y(month_vec(i),32));
    A_q(i) = y(month_vec(i), 19) + (y(month_vec(i+1),33)-y(month_vec(i),33));
end

L_q(1, 1:4) = 0;
L_q(1, 10:12) = 0;
N_q(1, 1:4) = 0;
N_q(1, 10:12) = 0;
A_q(1, 1:3) = 0;
A_q(1, 11:12) = 0;

%% Data Processing - Infected Densities

H_Ss = y(:,25); H_Si = y(:,26); H_Sr = y(:,27);
H_Us = y(:,28); H_Ui = y(:,29); H_Ur = y(:,30);

L_qs = y(:,1);  L_qi = y(:,13);     N_qs = y(:,4);  N_qi = y(:,16);
A_qs = y(:,7);  A_qi = y(:,19);

Prev_L = 100*L_qi./(L_qs + L_qi);       Prev_N = 100*N_qi./(N_qs + N_qi);
Prev_A = 100*A_qi./(A_qs + A_qi);

Prev_S = 100*(H_Si+H_Sr)./(H_Ss + H_Si + H_Sr);
Prev_U = 100*(H_Ui+H_Ur)./(H_Us + H_Ui + H_Ur);

[mean(Prev_A), mean(Prev_S), mean(Prev_U)]
[max(L_q), max(N_q), max(A_q)]

%% Data Figures
% 
figure
yyaxis left
plot(month_vec(1:end-1), A_q, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
ylabel({'Density of questing adult','\it{Hyalomma lusitanicum}, (per sq km)'})
% ylim([0 1e+03])
yyaxis right
plot(month_vec(1:end-1), L_q, ':', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
hold on
plot(month_vec(1:end-1), N_q*10, '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
ylabel({'Density of questing immature', '\it{Hyalomma lusitanicum}, (per sq km)'})
% ylim([0 4e+04])
xlim([0 365])
xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec', 'Jan'})
legend('A_q','L_q', 'N_q (x10)')
grid on
ax = gca;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
set(gca,'box','off')
set(gca, 'position', [0.12 0.2 0.75 0.5])
ax.FontSize = 11;

% ------------------------
figure
subplot(3,1,1)
plot(t, Prev_L)
hold on
plot(t, Prev_N)
plot(t, Prev_A)
xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
xticklabels('')
xlim([0 365])
%ylim([0 10])
%yticks([0 5 10])
set(gca,'box','off')
set(gca, 'position', [0.11 0.68 0.4 0.23])
legend('L_q','N_q','A_q' , 'Position', [0.55 0.68 0.1 0.1])
ax = gca;
ax.FontSize = 11;

subplot(3,1,2)
plot(t, Prev_S)
xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
xticklabels('')
xlim([0 365])
%ylim([4 10])
%yticks([4 6 8 10])
ylabel('Sero/Prevalences (%)')
set(gca,'box','off')
set(gca, 'position', [0.11 0.41 0.4 0.23])
ax = gca;
ax.FontSize = 11;

subplot(3,1,3)
plot(t, Prev_U)
xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec'})
xlim([0 365])
xlabel('Time (months)')
%ylim([38 44])
%yticks([38 40 42 44])
set(gca,'box','off')
set(gca, 'position', [0.11 0.14 0.4 0.23])
ax = gca;
ax.FontSize = 11;

