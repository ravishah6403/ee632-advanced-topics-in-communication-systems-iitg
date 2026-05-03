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

% SINR Thresholds
gamma_th_n_DL = 0.9;
gamma_th_f_DL = 0.5;

% Part (b) - Downlink Near User
var_SDn = 1;
var_fn = 1.2;
var_nn = 1.1;

% Channel Formation
h_SDn = sqrt(var_SDn/2) * (randn(1,N) + 1i*randn(1,N));
h_fn = sqrt(var_fn/2) * (randn(1,N) + 1i*randn(1,N));
h_nn = sqrt(var_nn/2) * (randn(1,N) + 1i*randn(1,N));

g_SDn = abs(h_SDn).^2;
g_fn = abs(h_fn).^2;
g_nn = abs(h_nn).^2;

beta_SDn = var_SDn;
beta_fn = var_fn;
beta_nn = var_nn;


% Analytical Calculation
phi = max(gamma_th_n_DL ./ (a_n .* (P_S ./ var_n)), gamma_th_f_DL ./ ((a_f - gamma_th_f_DL*a_n) .* (P_S ./ var_n)));
P_out_Dn_al = 1 - (beta_SDn ./ (beta_SDn + beta_fn*(P_Uf/var_n).*phi)) .* (beta_SDn ./ (beta_SDn + beta_nn*(P_Un/var_n).*phi)) .* exp(-phi/beta_SDn);

% Part (c) - Downlink Far User
var_SDf = 1.5;
var_ff = 1.7;
var_nf = 1.5;

h_SDf = sqrt(var_SDf/2) * (randn(1,N) + 1i*randn(1,N));
h_ff = sqrt(var_ff/2) * (randn(1,N) + 1i*randn(1,N));
h_nf = sqrt(var_nf/2) * (randn(1,N) + 1i*randn(1,N));

g_SDf = abs(h_SDf).^2;
g_ff = abs(h_ff).^2;
g_nf = abs(h_nf).^2;

beta_SDf = var_SDf;
beta_ff = var_ff;
beta_nf = var_nf;



% Analytical Calculation
theta = gamma_th_f_DL ./ ((a_f - gamma_th_f_DL*a_n) .* (P_S ./ var_n));
P_out_Df_al = 1 - (beta_SDf ./ (beta_SDf + beta_ff*(P_Uf/var_n).*theta)) .* (beta_SDf ./ (beta_SDf + beta_nf*(P_Un/var_n).*theta)) .* exp(-theta/beta_SDf);

% Simulation Loop
P_out_Df_sim = zeros(1, P_num);
P_out_Dn_sim = zeros(1, P_num);
for i = 1:P_num
    P = P_S(i);

    % SINR for Near User in Downlink NOMA with SIC
    gamma_Dn_xf = (a_f * P .* g_SDn) ./ (a_n * P .* g_SDn + P_Uf .* g_fn + P_Un .* g_nn + var_n);
    gamma_Dn_xn = (a_n * P .* g_SDn) ./ (P_Uf .* g_fn + P_Un .* g_nn + var_n);

    % Near User Downlink Outage Probability
    P_out_Dn_sim(i) = 1 - mean(gamma_Dn_xn > gamma_th_n_DL & gamma_Dn_xf > gamma_th_f_DL);
    
    % SINR for Far User in Downlink NOMA
    gamma_Df_xf = (a_f * P .* g_SDf) ./ (a_n * P .* g_SDf + P_Uf .* g_ff + P_Un .* g_nf + var_n);

    % Far User Downlink Outage Probability
    P_out_Df_sim(i) = mean(gamma_Df_xf < gamma_th_f_DL);
end

% Plotting
f_bc = figure;
semilogy(P_S_dBm_range, P_out_Dn_al, 'r', LineWidth=1); hold on;
semilogy(P_S_dBm_range, P_out_Dn_sim, 'g^'); 
semilogy(P_S_dBm_range, P_out_Df_al, 'm', LineWidth=1);
semilogy(P_S_dBm_range, P_out_Df_sim, 'go');
grid on;
xlabel('P_S (dB)');
ylabel('Outage Probability');
legend('Analytical (Near)', 'Monte Carlo (Near)', 'Analytical (Far)', 'Monte Carlo (Far)');
title('NOMA Downlink Outage Probability');
saveas(f_bc, 'Q1_b_c.fig');