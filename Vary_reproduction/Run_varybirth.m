%%% Script to produce Figure_2a, where we measure the physical prevalence
% /seroprevalences levels and densities of infected questing ticks for 
% changes in reproduction rate.
% proportional difference 

clearvars

%%
% ode specifics
options = odeset('Refine', 1, 'NonNegative',1:33, 'RelTol', 1e-07, 'AbsTol', 1e-07);
tspan = 0:1:365;
month_vec = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366];

% scaling factor
k_vector = [0.1:0.01:0.5, 0.5:0.1:2];

%% What to keep track of

% Mean yearly questing tick prevalence
Prev_L = zeros(1, length(k_vector)); Prev_N = zeros(1, length(k_vector));
Prev_A = zeros(1, length(k_vector));

% Mean yearly host seroprevalence
Prev_S = zeros(1, length(k_vector)); Prev_U = zeros(1, length(k_vector));

% Peak infected questing tick densities
L_qi = zeros(1, length(k_vector)); N_qi = zeros(1, length(k_vector));
A_qi = zeros(1, length(k_vector));

%% Storing the data
for k = 1:length(k_vector)

    y0 = [10000*ones(1,24), 80, 20, 0, 10, 10, 0, zeros(1,3)];
    for i = 1:50
        [t,y] = ode15s(@(t,y) Spanish_ticks_INF(t, y, k_vector(k)), tspan, y0, options);
        y0 = y(end,:);
        y0(10) = 0.05*y0(9);
        y0(22) = 0.05*y0(21);
        y0([11,12,23,24]) = 0;
        y0([31, 32, 33])=0;
    end
    
    L = zeros(1,12);    N = zeros(1,12);     A = zeros(1,12);
    for i = 1:12
        L(i) = y(month_vec(i), 13) + (y(month_vec(i+1),31)-y(month_vec(i),31));
        N(i) = y(month_vec(i), 16) + (y(month_vec(i+1),32)-y(month_vec(i),32));
        A(i) = y(month_vec(i), 19) + (y(month_vec(i+1),33)-y(month_vec(i),33));
    end
    L(1, 1:4) = 0;
    L(1, 10:12) = 0;
    N(1, 1:4) = 0;
    N(1, 10:12) = 0;
    A(1, 1:3) = 0;
    A(1, 11:12) = 0;

    L_qi(k) = max(L);    N_qi(k) = max(N);    A_qi(k) = max(A);

    Prev_L(k) = mean(100*y(:,13)./(y(:,1) + y(:,13)));
    Prev_N(k) = mean(100*y(:,16)./(y(:,4) + y(:,16)));
    Prev_A(k) = mean(100*y(:,19)./(y(:,7) + y(:,19)));

    Prev_S(k) = mean(100*(y(:,26) + y(:,27))./(y(:,25) + y(:,26) + y(:,27)));
    Prev_U(k) = mean(100*(y(:,29) + y(:,30))./(y(:,28) + y(:,29) + y(:,30)));

    Prev_L(isnan(Prev_L)) = 0;
    Prev_N(isnan(Prev_N)) = 0;
    Prev_A(isnan(Prev_A)) = 0;
    Prev_S(isnan(Prev_S)) = 0;
    Prev_U(isnan(Prev_U)) = 0;
end

%% FINAL FIGURES
close all
fig = figure;
plot(k_vector, Prev_L, 'k:', 'LineWidth', 1);
hold on
plot(k_vector, Prev_N, 'k-.', 'LineWidth', 1);
plot(k_vector, Prev_A, 'k-', 'LineWidth', 1);
plot(k_vector, Prev_S, '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plot(k_vector, Prev_U, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
ylim([0 100])
xlim([0 2])
xticks([0.25 0.5 0.75 1 1.25 1.5 1.75 2])
ylabel('Average yearly sero/prevalences (%)')
xlabel('Tick reproduction scaling factor, k')
set(gca,'box','off')
%set(gca, 'position', [0.11 0.14 0.4 0.23])
grid on
ax = gca;
ax.GridColor = [0.5 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
%legend('Larvae', 'Nymphs','Adults','Small mammals','Ungulates', 'Location', 'East')
ax.FontSize = 11;

fig2 = figure;
plot(k_vector, A_qi*10, 'k-', 'LineWidth', 1.5);
hold on
plot(k_vector, N_qi*10, 'k-.',  'LineWidth', 1.5);
plot(k_vector, L_qi, 'k:', 'LineWidth', 1.5);
ylabel('Peak infected questing tick densities', 'Color', 'k')
xlim([0.5 1.5])
xticks([0.25 0.5 0.75 1 1.25 1.5 1.75])
ylim([0 6e+05])
yticks([0 1.5e+05 3e+05 4.5e+05 6e+05])
xlabel('Tick reproduction scaling factor, k')
set(gca,'box','off')
%set(gca, 'position', [0.11 0.14 0.4 0.23])
grid on
ax = gca;
ax.GridColor = [0.5 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
legend('Adults (x10)', 'Nymphs (x10)', 'Larvae', 'Location', 'NorthWest')
ax.FontSize = 11;
