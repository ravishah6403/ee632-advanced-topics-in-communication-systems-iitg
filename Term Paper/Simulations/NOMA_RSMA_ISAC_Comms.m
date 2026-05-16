clear; clc; close all;

%% ========================================================================
%  1. SYSTEM PARAMETERS
%  ========================================================================
% Simulation & SNR
N        = 1e6;             % Monte-Carlo iterations
N0       = 1;               % Noise variance
P_S_dB   = 0:5:50;          % Transmit power in dB
P_S      = 10.^((P_S_dB)./10);% Transmit power in linear scale
gamma_S  = P_S ./ N0;       % Transmit SNR

% Channel Parameters (Nakagami-m)
m_SDn     = 2;              % Shape factor (Near User)
m_SDf     = 1;              % Shape factor (Far User)
Omega_SDn = 1;              % Spread factor (Near User)
Omega_SDf = 1;              % Spread factor (Far User)

% Target SINR Thresholds
gamma_thr  = 0.5;           
gamma_thn  = 0.6;           
gamma_thf  = 0.8;           
gamma_thc  = 0.6;           
gamma_thnp = 0.26;          
gamma_thfp = 0.42;          

% Hardware Impairment & Power Allocation
epsilon = 0.01;             % Hardware impairment level
sigma_SDf = [0; 0.1];       % Far user hardware noise [Ideal; Impaired]
sigma_SDn = [0; 0.1];       % Near user hardware noise [Ideal; Impaired]

% Power Coefficients (RSMA & NOMA)
a_c  = 0.3;                 % RSMA Common
a_np = 0.05;                % RSMA Near Private
a_fp = 0.15;                % RSMA Far Private
a_r  = 0.5;                 % Radar/Sensing
a_n  = 0.2;                 % NOMA Near
a_f  = 0.3;                 % NOMA Far

M = 2;                      % Number of Antennas             

%% ========================================================================
%  2. CHANNEL GENERATION
%  ========================================================================
g_SDf = gamrnd(m_SDf*M, Omega_SDf/m_SDf, 1, N);
g_SDn = gamrnd(m_SDn*M, Omega_SDn/m_SDn, 1, N);

%% ========================================================================
%  3. THEORETICAL OUTAGE PROBABILITY (OP)
%  ========================================================================
% --- RSMA Denominators (Interference Limits) ---
den_RS_r  = a_r - (a_np + a_fp + a_c)*gamma_thr;
den_RS_c  = a_c - (a_np + a_fp + epsilon*a_r)*gamma_thc;
den_RS_fp = a_fp - (a_np + epsilon*(a_c + a_r))*gamma_thfp;
den_RS_np = a_np - (a_fp + epsilon*(a_c + a_r))*gamma_thnp;

% --- RSMA Thresholds (phi -> Far, psi -> Near) ---
phi_r  = gamma_thr  * (sigma_SDf.^2 + 1./gamma_S) / den_RS_r;
phi_c  = gamma_thc  * (sigma_SDf.^2 + 1./gamma_S) / den_RS_c;
phi_fp = gamma_thfp * (sigma_SDf.^2 + 1./gamma_S) / den_RS_fp;

psi_r  = gamma_thr  * (sigma_SDn.^2 + 1./gamma_S) / den_RS_r;
psi_c  = gamma_thc  * (sigma_SDn.^2 + 1./gamma_S) / den_RS_c;
psi_np = gamma_thnp * (sigma_SDn.^2 + 1./gamma_S) / den_RS_np;

% Impossible SINR targets (Negative Denominators -> Outage = 1)
phi_r(phi_r < 0) = Inf; phi_c(phi_c < 0)   = Inf; phi_fp(phi_fp < 0) = Inf;
psi_r(psi_r < 0) = Inf; psi_c(psi_c < 0)   = Inf; psi_np(psi_np < 0) = Inf;

% Max thresholds for OP calculation
phi = max(max(phi_r, phi_c), phi_fp);
psi = max(max(psi_r, psi_c), psi_np);

P_RS_Df = gammainc((m_SDf*phi)/Omega_SDf, m_SDf*M, 'lower');
P_RS_Dn = gammainc((m_SDn*psi)/Omega_SDn, m_SDn*M, 'lower');

% --- NOMA Denominators ---
den_NO_1 = a_r - (a_n + a_f)*gamma_thr;
den_NO_2 = a_f - (a_n + epsilon*a_r)*gamma_thf; 
den_NO_3 = a_n - epsilon*(a_f + a_r)*gamma_thn;

% --- NOMA Thresholds ---
phi_d_1 = gamma_thr * (sigma_SDf.^2 + 1./gamma_S) / den_NO_1;
phi_d_2 = gamma_thf * (sigma_SDf.^2 + 1./gamma_S) / den_NO_2;

psi_d_1 = gamma_thr * (sigma_SDn.^2 + 1./gamma_S) / den_NO_1;
psi_d_2 = gamma_thf * (sigma_SDn.^2 + 1./gamma_S) / den_NO_2;
psi_d_3 = gamma_thn * (sigma_SDn.^2 + 1./gamma_S) / den_NO_3;

phi_d_1(phi_d_1 < 0) = Inf; phi_d_2(phi_d_2 < 0) = Inf;
psi_d_1(psi_d_1 < 0) = Inf; psi_d_2(psi_d_2 < 0) = Inf; psi_d_3(psi_d_3 < 0) = Inf;

% Max thresholds
phi_d = max(phi_d_1, phi_d_2);
psi_d = max(max(psi_d_1, psi_d_2), psi_d_3);

P_NO_Df = gammainc((m_SDf*phi_d)/Omega_SDf, m_SDf*M, 'lower'); 
P_NO_Dn = gammainc((m_SDn*psi_d)/Omega_SDn, m_SDn*M, 'lower');

%% ========================================================================
%  4. THEORETICAL ASYMPTOTIC ER (High SNR Ceiling)
%  ========================================================================
% --- NOMA Asymptotic Limits ---
SINR_NO_xn_inf = a_n ./ (epsilon*(a_r + a_f) + sigma_SDn.^2 ./ g_SDn);
R_NO_Dn_h = mean(log2(1 + SINR_NO_xn_inf), 2);

SINR_NO_xf_SDn_inf = a_f ./ (a_n + epsilon*a_r + sigma_SDn.^2 ./ g_SDn);
SINR_NO_xf_SDf_inf = a_f ./ (a_n + epsilon*a_r + sigma_SDf.^2 ./ g_SDf);
R_NO_Df_h = mean(log2(1 + min(SINR_NO_xf_SDn_inf, SINR_NO_xf_SDf_inf)), 2);

% --- RSMA Asymptotic Limits ---
SINR_RS_xc_SDn_inf = a_c ./ (a_np + a_fp + epsilon*a_r + sigma_SDn.^2 ./ g_SDn);
SINR_RS_xc_SDf_inf = a_c ./ (a_np + a_fp + epsilon*a_r + sigma_SDf.^2 ./ g_SDf);
R_RS_C_h = mean(log2(1 + min(SINR_RS_xc_SDn_inf, SINR_RS_xc_SDf_inf)), 2);

SINR_RS_xnp_inf = a_np ./ (epsilon*(a_c + a_r) + a_fp + sigma_SDn.^2 ./ g_SDn);
R_RS_Dn_h = (R_RS_C_h / 2) + mean(log2(1 + SINR_RS_xnp_inf), 2);

SINR_RS_xfp_inf = a_fp ./ (epsilon*(a_c + a_r) + a_np + sigma_SDf.^2 ./ g_SDf);
R_RS_Df_h = (R_RS_C_h / 2) + mean(log2(1 + SINR_RS_xfp_inf), 2);

%% ========================================================================
%  5. MONTE-CARLO SIMULATIONS
%  ========================================================================
num_snr = numel(P_S);
P_out_NO_Dn = zeros(2, num_snr);
P_out_NO_Df = zeros(2, num_snr);
P_out_RS_Dn = zeros(2, num_snr);
P_out_RS_Df = zeros(2, num_snr);

R_NO_Dn = zeros(2, num_snr);
R_NO_Df = zeros(2, num_snr);
R_RS_Dn = zeros(2, num_snr);
R_RS_Df = zeros(2, num_snr);

disp('Running Monte-Carlo Loop...');
for i = 1:num_snr
    g_S = gamma_S(i);
    
    % -----------------------------------------------------------
    % NOMA Simulation
    % -----------------------------------------------------------
    % Near User
    gamma_NO_xr_SDn = (g_SDn*a_r*g_S) ./ ((g_SDn*a_n + g_SDn*a_f + sigma_SDn.^2)*g_S + 1);
    gamma_NO_xf_SDn = (g_SDn*a_f*g_S) ./ ((g_SDn*a_n + g_SDn*epsilon*a_r + sigma_SDn.^2)*g_S + 1);
    gamma_NO_xn_SDn = (g_SDn*a_n*g_S) ./ ((g_SDn*a_f + g_SDn*a_r)*epsilon*g_S + sigma_SDn.^2*g_S + 1);
    
    NO_Dn_success = (gamma_NO_xr_SDn > gamma_thr) & (gamma_NO_xf_SDn > gamma_thf) & (gamma_NO_xn_SDn > gamma_thn);
    P_out_NO_Dn(:, i) = 1 - mean(NO_Dn_success, 2);
    R_NO_Dn(:, i)     = mean(log2(1 + gamma_NO_xn_SDn), 2);
    
    % Far User
    gamma_NO_xr_SDf = (g_SDf*a_r*g_S) ./ ((g_SDf*a_n + g_SDf*a_f + sigma_SDf.^2)*g_S + 1);
    gamma_NO_xf_SDf = (g_SDf*a_f*g_S) ./ ((g_SDf*a_n + g_SDf*epsilon*a_r + sigma_SDf.^2)*g_S + 1);
    
    NO_Df_success = (gamma_NO_xr_SDf > gamma_thr) & (gamma_NO_xf_SDf > gamma_thf);
    P_out_NO_Df(:, i) = 1 - mean(NO_Df_success, 2);
    R_NO_Df(:, i)     = mean(log2(1 + min(gamma_NO_xf_SDn, gamma_NO_xf_SDf)), 2);

    % -----------------------------------------------------------
    % RSMA Simulation
    % -----------------------------------------------------------
    % Near User
    gamma_RS_xr_SDn  = (a_r*g_SDn*g_S)  ./ ((a_np*g_SDn + a_fp*g_SDn + a_c*g_SDn + sigma_SDn.^2)*g_S + 1);
    gamma_RS_xc_SDn  = (a_c*g_SDn*g_S)  ./ ((a_np*g_SDn + a_fp*g_SDn + epsilon*a_r*g_SDn + sigma_SDn.^2)*g_S + 1);
    gamma_RS_xnp_SDn = (a_np*g_SDn*g_S) ./ ((epsilon*a_c*g_SDn + a_fp*g_SDn + epsilon*a_r*g_SDn + sigma_SDn.^2)*g_S + 1);
    
    RS_Dn_success = (gamma_RS_xr_SDn > gamma_thr) & (gamma_RS_xc_SDn > gamma_thc) & (gamma_RS_xnp_SDn > gamma_thnp);
    P_out_RS_Dn(:, i) = 1 - mean(RS_Dn_success, 2);
    
    % Far User
    gamma_RS_xr_SDf  = (a_r*g_SDf*g_S)  ./ ((a_np*g_SDf + a_fp*g_SDf + a_c*g_SDf + sigma_SDf.^2)*g_S + 1);
    gamma_RS_xc_SDf  = (a_c*g_SDf*g_S)  ./ ((a_np*g_SDf + a_fp*g_SDf + epsilon*a_r*g_SDf + sigma_SDf.^2)*g_S + 1);
    gamma_RS_xfp_SDf = (a_fp*g_SDf*g_S) ./ ((epsilon*a_c*g_SDf + a_np*g_SDf + epsilon*a_r*g_SDf + sigma_SDf.^2)*g_S + 1);
    
    RS_Df_success = (gamma_RS_xr_SDf > gamma_thr) & (gamma_RS_xc_SDf > gamma_thc) & (gamma_RS_xfp_SDf > gamma_thfp);
    P_out_RS_Df(:, i) = 1 - mean(RS_Df_success, 2);
    
    % RSMA Ergodic Rates
    RS_C_ER       = min(log2(1 + gamma_RS_xc_SDn), log2(1 + gamma_RS_xc_SDf));
    R_RS_Dn(:, i) = mean((RS_C_ER/2) + log2(1 + gamma_RS_xnp_SDn), 2);
    R_RS_Df(:, i) = mean((RS_C_ER/2) + log2(1 + gamma_RS_xfp_SDf), 2);
end
disp('Simulation Complete.');

%% ========================================================================
%  6. PLOTTING: OUTAGE PROBABILITY
%  ========================================================================
figure('Name', 'Outage Probability', 'Position', [100, 100, 800, 600]);

% --- Case 1: Ideal Hardware (sigma = 0) ---
h_NO_th_0 = semilogy(P_S_dB, P_NO_Df(1,:), 'k-', 'LineWidth', 1.8); hold on;
            semilogy(P_S_dB, P_NO_Dn(1,:), 'k-', 'LineWidth', 1.8); 
h_NO_dn_0 = semilogy(P_S_dB, P_out_NO_Dn(1,:), 'k^', 'MarkerSize', 8, 'LineWidth', 1.2);
h_NO_df_0 = semilogy(P_S_dB, P_out_NO_Df(1,:), 'ko', 'MarkerSize', 8, 'LineWidth', 1.2);

h_RS_th_0 = semilogy(P_S_dB, P_RS_Df(1,:), 'r-', 'LineWidth', 1.8);
            semilogy(P_S_dB, P_RS_Dn(1,:), 'r-', 'LineWidth', 1.8);
h_RS_dn_0 = semilogy(P_S_dB, P_out_RS_Dn(1,:), 'r^', 'MarkerSize', 8, 'LineWidth', 1.2);
h_RS_df_0 = semilogy(P_S_dB, P_out_RS_Df(1,:), 'ro', 'MarkerSize', 8, 'LineWidth', 1.2);

% --- Case 2: Impaired Hardware (sigma = 0.1) ---
h_NO_th_1 = semilogy(P_S_dB, P_NO_Df(2,:), 'k--', 'LineWidth', 1.8);
            semilogy(P_S_dB, P_NO_Dn(2,:), 'k--', 'LineWidth', 1.8);
h_NO_dn_1 = semilogy(P_S_dB, P_out_NO_Dn(2,:), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
h_NO_df_1 = semilogy(P_S_dB, P_out_NO_Df(2,:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

h_RS_th_1 = semilogy(P_S_dB, P_RS_Df(2,:), 'r--', 'LineWidth', 1.8);
            semilogy(P_S_dB, P_RS_Dn(2,:), 'r--', 'LineWidth', 1.8);
h_RS_dn_1 = semilogy(P_S_dB, P_out_RS_Dn(2,:), 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
h_RS_df_1 = semilogy(P_S_dB, P_out_RS_Df(2,:), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Formatting
ylim([1e-7, 1]); xlim([0, 50]);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('SNR (dB)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Outage Probability', 'FontSize', 14, 'FontWeight', 'bold');
title('Outage Probability vs SNR', 'FontSize', 16);

% OP Legend
lgd_op = legend([h_NO_th_0, h_RS_th_0, h_NO_th_1, h_RS_th_1, ...
                 h_NO_dn_0, h_RS_dn_0, h_NO_df_0, h_RS_df_0], ...
    {'NOMA Theory (\sigma = 0)',   'RSMA Theory (\sigma = 0)', ...
     'NOMA Theory (\sigma = 0.1)', 'RSMA Theory (\sigma = 0.1)', ...
     'NOMA Sim D_n (\sigma = 0)',  'RSMA Sim D_n (\sigma = 0)', ...
     'NOMA Sim D_f (\sigma = 0)',  'RSMA Sim D_f (\sigma = 0)'}, ...
    'Location', 'southwest', 'NumColumns', 2, 'FontSize', 11);

%% ========================================================================
%  7. PLOTTING: ERGODIC RATE
%  ========================================================================
figure('Name', 'Ergodic Rate', 'Position', [920, 100, 800, 600]);

% --- Case 1: Ideal Hardware (sigma = 0) ---
h_er_NO_dn_0 = plot(P_S_dB, R_NO_Dn(1,:), 'k-^', 'LineWidth', 1.8, 'MarkerSize', 8); hold on;
h_er_NO_df_0 = plot(P_S_dB, R_NO_Df(1,:), 'k-o', 'LineWidth', 1.8, 'MarkerSize', 8);
h_er_RS_dn_0 = plot(P_S_dB, R_RS_Dn(1,:), 'r-^', 'LineWidth', 1.8, 'MarkerSize', 8);
h_er_RS_df_0 = plot(P_S_dB, R_RS_Df(1,:), 'r-o', 'LineWidth', 1.8, 'MarkerSize', 8);

% --- Case 2: Impaired Hardware (sigma = 0.1) ---
h_er_NO_dn_1 = plot(P_S_dB, R_NO_Dn(2,:), 'k--^', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
h_er_NO_df_1 = plot(P_S_dB, R_NO_Df(2,:), 'k--o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
h_er_RS_dn_1 = plot(P_S_dB, R_RS_Dn(2,:), 'r--^', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
h_er_RS_df_1 = plot(P_S_dB, R_RS_Df(2,:), 'r--o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% --- Asymptotic Lines ---
yline(R_NO_Dn_h(1), 'k:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off');
yline(R_RS_Dn_h(1), 'r:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off');
yline(R_NO_Dn_h(2), 'k:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off'); 
yline(R_RS_Dn_h(2), 'r:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off'); 

yline(R_NO_Df_h(1), 'k:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off');
yline(R_RS_Df_h(1), 'r:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off');
yline(R_NO_Df_h(2), 'k:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off'); 
yline(R_RS_Df_h(2), 'r:', 'LineWidth', 2, 'Alpha', 0.6, 'HandleVisibility', 'off'); 

% Formatting
xlim([0, 50]);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('SNR (dB)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Ergodic Rate (bps/Hz)', 'FontSize', 14, 'FontWeight', 'bold');
title('Ergodic Rate vs SNR', 'FontSize', 16);

% ER Legend
lgd_er = legend([h_er_NO_dn_0, h_er_RS_dn_0, h_er_NO_df_0, h_er_RS_df_0, ...
                 h_er_NO_dn_1, h_er_RS_dn_1, h_er_NO_df_1, h_er_RS_df_1], ...
    {'NOMA D_n (\sigma = 0)',   'RSMA D_n (\sigma = 0)', ...
     'NOMA D_f (\sigma = 0)',   'RSMA D_f (\sigma = 0)', ...
     'NOMA D_n (\sigma = 0.1)', 'RSMA D_n (\sigma = 0.1)', ...
     'NOMA D_f (\sigma = 0.1)', 'RSMA D_f (\sigma = 0.1)'}, ...
    'Location', 'northwest', 'NumColumns', 2, 'FontSize', 11);