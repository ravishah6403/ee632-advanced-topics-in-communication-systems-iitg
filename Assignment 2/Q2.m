clear; clc;

% Simulation Iteraions
N = 1e6;

% Noise Variance
var_n = 1;

% Power Range
P_dBm_range = 0:2:50;
P_num = numel(P_dBm_range);

% Channel Parameters
var_h1 = 1;
var_h2 = 1.5;

h1 = sqrt(var_h1 / 2) * (randn(1,N) + 1i*randn(1,N));
h2 = sqrt(var_h2 / 2) * (randn(1,N) + 1i*randn(1,N));

g1 = abs(h1).^2;
g2 = abs(h2).^2;

% Power Allocation
a_c = 0.5;
a_p1 = 0.3;
a_p2 = 0.2;

% SINR Thresholds
gamma_th_c = 0.9;
gamma_th_p = 0.5;

% Arrays to store results
R_c = zeros(1, P_num);
R_p1 = zeros(1, P_num);
R_p2 = zeros(1, P_num);

P_out1 = zeros(1, P_num);
P_out2 = zeros(1, P_num);

for i = 1:P_num
    P = (1/1000) * 10^(P_dBm_range(i) / 10);
    
    % SINR for D1
    gamma_c1 = (a_c * P .* g1) ./ ((a_p1 + a_p2) * P .* g1 + var_n);
    gamma_p1 = (a_p1 * P .* g1) ./ (a_p2 * P .* g1 + var_n);

    % SINR for D2
    gamma_c2 = (a_c * P .* g2) ./ ((a_p1 + a_p2) * P .* g2 + var_n);
    gamma_p2 = (a_p2 * P .* g2) ./ (a_p1 * P .* g2 + var_n);
    
    % Rates Calculation
    R_c(i) = mean(log2(1 + min(gamma_c1, gamma_c2)));
    R_p1(i) = mean(log2(1 + gamma_p1));
    R_p2(i) = mean(log2(1 + gamma_p2));

    % Outage Probabilty Calculation
    outage_D1 = (gamma_c1 < gamma_th_c) | (gamma_p1 < gamma_th_p);
    outage_D2 = (gamma_c2 < gamma_th_c) | (gamma_p2 < gamma_th_p);

    P_out1(i) = mean(outage_D1);
    P_out2(i) = mean(outage_D2);
end


f_2 = figure;

% Plot 1: Ergodic Rates
subplot(1, 2, 1);
plot(P_dBm_range, R_c, 'w-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(P_dBm_range, R_p1, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(P_dBm_range, R_p2, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 6);
title('Ergodic Rates vs Transmit Power');
xlabel('Transmit Power P (dBm)');
ylabel('Ergodic Rate (bits/s/Hz)');
grid on;
legend('Common Rate (R_c)', 'D1 Private Rate (R_{p,1})', 'D2 Private Rate (R_{p,2})', 'Location', 'best');

% Plot 2: Outage Probabilities
subplot(1, 2, 2);
semilogy(P_dBm_range, P_out1, 'b-s', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(P_dBm_range, P_out2, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 6);
title('Outage Probabilities vs Transmit Power');
xlabel('Transmit Power P (dBm)');
ylabel('Outage Probability');
grid on;
legend('D1 Outage Probability', 'D2 Outage Probability', 'Location', 'best');

saveas(f_2, 'Q2.fig');