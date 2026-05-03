clear; clc;

% System Parameters
N = 1e6;
var_n = 1;
P_S_dBm_range = 0:2:50;
P_S = 10.^(P_S_dBm_range./10);
P_num = numel(P_S);
P_Uf = 3;
P_Un = 6;
a_f = 0.7;
a_n = 0.3;


% Part (d) - Uplink
var_SUn = 1.3;
var_SUf = 1.7;
var_LI = 1;

% Channel Formation
h_SUn = sqrt(var_SUn / 2) * (randn(1, N) + 1i*randn(1, N));
h_SUf = sqrt(var_SUf / 2) * (randn(1, N) + 1i*randn(1, N));
h_LI = sqrt(var_LI / 2) * (randn(1, N) + 1i*randn(1, N));

g_SUn = abs(h_SUn).^2;
g_SUf = abs(h_SUf).^2;
g_LI = abs(h_LI).^2;

% SINR Thresholds
gamma_th_n_UL = 0.9;
gamma_th_f_UL = 0.7;

P_out_Un_sim = zeros(1, P_num);
P_out_Uf_sim = zeros(1, P_num);

for i = 1:P_num
    P = P_S(i);

    % SINR Calculations
    gamma_BS_sn = (P_Un * g_SUn) ./ (P_Uf * g_SUf + P * g_LI + var_n);
    gamma_BS_sf = (P_Uf * g_SUf) ./ (P * g_LI + var_n);

    P_out_Un_sim(i) = mean(gamma_BS_sn < gamma_th_n_UL);
    P_out_Uf_sim(i) = 1 - mean(gamma_BS_sn > gamma_th_n_UL & gamma_BS_sf > gamma_th_f_UL);
end

f_d = figure;
semilogy(P_S_dBm_range, P_out_Un_sim, 'g^'); hold on;
semilogy(P_S_dBm_range, P_out_Uf_sim, 'go');
grid on;
xlabel('P_S (dB)');
ylabel('Outage Probability');
legend('Monte Carlo (Near)', 'Monte Carlo (Far)');
title('NOMA Uplink Outage Probability');
saveas(f_d, 'Q1_d.fig');