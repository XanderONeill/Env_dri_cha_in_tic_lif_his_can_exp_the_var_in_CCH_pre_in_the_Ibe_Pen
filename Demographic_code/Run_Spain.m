%%% Script to run the Spanish CCHFv model in the absence of infection

%close all
clearvars

%% ODE Solver
options = odeset('Refine', 1, 'NonNegative',1:15, 'RelTol', 1e-07, 'AbsTol', 1e-07);
% options = odeset('Refine', 1, 'NonNegative',1:15, 'RelTol', 1e-09, 'AbsTol', 1e-09);

tspan = 0:1:365;
%            Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec, Jan
month_vec = [1,   32,  60,  91,  121, 152, 182, 213, 244, 274, 305, 335, 366];

y0 = ones(1,15)*100;
for i = 1:50
    [t,y] = ode15s(@(t,y) Spanish_ticks(t,y), tspan, y0, options);
    y0 = y(end,:);
    y0(10) = 0.05*y0(9);
    y0(11) = 0;
    y0(12) = 0;
    y0([13, 14, 15])=0;
end

%% Data Processing

L_q = zeros(1,12);
N_q = zeros(1,12);
A_q = zeros(1,12);

for i = 1:12
    L_q(i) = y(month_vec(i),1) + y(month_vec(i+1),13)-y(month_vec(i),13);
    N_q(i) = y(month_vec(i),4) + y(month_vec(i+1),14)-y(month_vec(i),14);
    A_q(i) = y(month_vec(i),7) + y(month_vec(i+1),15)-y(month_vec(i),15);
end
L_q(1, 1:4) = 0;
L_q(1, 10:12) = 0;
N_q(1, 1:4) = 0;
N_q(1, 10:12) = 0;
A_q(1, 1:3) = 0;
A_q(1, 11:12) = 0;

%% Figure producing
%close all
%%% FINAL FIGURE - OPTION 1
figure
yyaxis left
plot(month_vec(1:end-1), A_q, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
ylabel({'Density of questing adult','\it{Hyalomma lusitanicum}, (per sq km)'})
ylim([0 6e+04])
yticks([0 2e+04 4e+04 6e+04])
yyaxis right
plot(month_vec(1:end-1), L_q, ':', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
hold on
plot(month_vec(1:end-1), N_q*10, '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
ylabel({'Density of questing immature', '\it{Hyalomma lusitanicum}, (per sq km)'})
ylim([0 9e+06])
yticks([0 3e+06 6e+06 9e+06])
xlim([0 365])
xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec', 'Jan'})
% legend('A_q','L_q', 'N_q (x10)')
grid on
ax = gca;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
set(gca,'box','off')
set(gca, 'position', [0.12 0.2 0.75 0.5])
ax.FontSize = 11;


%%% FINAL FIGURE - OPTION 2
% figure
% yyaxis left
% plot(month_vec(2:end), N_q)
% hold on
% plot(month_vec(2:end), A_q*10)
% ylabel({'Density of questing nymph and adult', '\it{Hyalomma lusitanicum}, (per sq km)'})
% yyaxis right
% plot(month_vec(2:end), L_q)
% ylabel({'Density of questing larvae','\it{Hyalomma lusitanicum}, (per sq km)'})
% xlim([0 365])
% xticks([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365])
% xticklabels({'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec', 'Jan'})
% legend('N_q', 'A_q (x10)','L_q')
% grid on
% ax = gca;
% ax.GridColor = [.5 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
% set(gca,'box','off')
% set(gca, 'position', [0.12 0.2 0.75 0.5])
% ax.FontSize = 11;