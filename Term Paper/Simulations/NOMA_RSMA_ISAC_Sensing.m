clear; clc; close all;

%% ========================================================================
%  1. SHARED SENSING PARAMETERS & TARGET SETUP
%  ========================================================================
N0 = 1;

% Antennas & Arrays
M = 16;  % Transmit antennas
N = 16;  % Receive antennas
lambda = 0.01; 
d = lambda / 2;

% Target Setup
theta_deg = [10, 45]; 
theta = theta_deg * pi / 180;
dist = [50, 10]; % True physical distances in meters
K = 2;

% Power allocations
a_r = 0.5; a_c = 0.3; a_np = 0.05; a_fp = 0.15; % RSMA
a_n = 0.2; a_f = 0.3;                           % NOMA

% -----------------------------------------------------------
% PHYSICALLY CORRECT RADAR PATH LOSS
% -----------------------------------------------------------
% Radar echo power decays proportionally to 1 / d^4
scale_factor = 1e4; 
beta_squared = scale_factor ./ (dist.^4 * M * N);
beta_k = sqrt(beta_squared);

% Build Steering Vectors & Target Response Matrix (G)
G_k = zeros(N, M, K);
R_t_all = cell(K, 1);

for t = 1:K
    a_rx = exp(1j * 2*pi*d/lambda * (0:N-1)' * sin(theta(t)));
    b_tx = exp(1j * 2*pi*d/lambda * (0:M-1)' * sin(theta(t)));
    G_k(:,:,t) = beta_k(t) * a_rx * b_tx.';
    
    % Store the isolated spatial covariance matrix for each target
    R_t_all{t} = G_k(:,:,t) * G_k(:,:,t)';
end

% Optimizer options for root finding
options = optimset('Display', 'off', 'TolX', 1e-4);

%% ========================================================================
%  2. CALCULATE Figure 8 (PoD vs SNR)
%  ========================================================================
disp('Calculating PoD vs SNR...');

% Expanded SNR range to capture both the strong and weak targets!
P_S_dB  = 10:2:80;           
P_S     = 10.^(P_S_dB./10);  
P_fa_target  = 1e-5;              

PoD_RS_fig8 = zeros(K, length(P_S));
PoD_NO_fig8 = zeros(K, length(P_S));

for t = 1:K
    R_t = R_t_all{t};
    for i = 1:length(P_S)
        curr_PS = P_S(i);
        
        % --- RSMA ---
        ncp_fa_RS = sqrt(abs(trace((a_np + a_fp)*curr_PS * R_t)) / N0);
        obj_fa_RS = @(xi) marcumq(ncp_fa_RS, sqrt(2*xi/N0), 3) - P_fa_target;
        xi_RS = fzero(obj_fa_RS, [0, 1e8], options); 
        ncp_d_RS = sqrt(abs(trace(curr_PS * R_t)) / N0);
        PoD_RS_fig8(t, i) = marcumq(ncp_d_RS, sqrt(2*xi_RS/N0), 5);
        
        % --- NOMA ---
        ncp_fa_NO = sqrt(abs(trace(a_n * curr_PS * R_t)) / N0);
        obj_fa_NO = @(xi) marcumq(ncp_fa_NO, sqrt(2*xi/N0), 3) - P_fa_target;
        xi_NO = fzero(obj_fa_NO, [0, 1e8], options);
        ncp_d_NO = sqrt(abs(trace(curr_PS * R_t)) / N0);
        PoD_NO_fig8(t, i) = marcumq(ncp_d_NO, sqrt(2*xi_NO/N0), 3);
    end
end

%% ========================================================================
%  3. CALCULATE Figure 9 (PoD vs PoFA - ROC Curves)
%  ========================================================================
disp('Calculating PoD vs PoFA at 40 dBm...');

% 1. CORRECT POWER CONVERSION: 40 dBm -> Watts
P_S_dBm = 40;                        % Power in dBm
P_S_W   = 10^((P_S_dBm - 30)/10);    % Convert dBm to Watts (10 Watts)
curr_PS = P_S_W;                     % Linear transmit power (assuming N0 = 1W)

PoFA_vec = linspace(0.001, 0.999, 21); 

PoD_RS_fig9 = zeros(K, length(PoFA_vec));
PoD_NO_fig9 = zeros(K, length(PoFA_vec));

for t = 1:K
    R_t = R_t_all{t};
    
    % --- RSMA Non-centrality parameters ---
    ncp_fa_RS = sqrt(abs(trace((a_np + a_fp)*curr_PS * R_t)) / N0);
    ncp_d_RS  = sqrt(abs(trace(curr_PS * R_t)) / N0);
    
    % --- NOMA Non-centrality parameters ---
    % PHYSICALLY CORRECT: Using a_n (near user power) for the H0 interference, 
    % rejecting the paper's Eq. 59 typo which incorrectly used a_f.
    ncp_fa_NO = sqrt(abs(trace(a_n * curr_PS * R_t)) / N0);
    ncp_d_NO  = sqrt(abs(trace(curr_PS * R_t)) / N0);
    
    for i = 1:length(PoFA_vec)
        curr_PoFA = PoFA_vec(i);
        
        % --- RSMA Detection ---
        obj_fa_RS = @(xi) marcumq(ncp_fa_RS, sqrt(2*xi/N0), 3) - curr_PoFA;
        xi_RS = fzero(obj_fa_RS, [0, 1e8], options); 
        PoD_RS_fig9(t, i) = marcumq(ncp_d_RS, sqrt(2*xi_RS/N0), 5);
        
        % --- NOMA Detection ---
        obj_fa_NO = @(xi) marcumq(ncp_fa_NO, sqrt(2*xi/N0), 3) - curr_PoFA;
        xi_NO = fzero(obj_fa_NO, [0, 1e8], options);
        PoD_NO_fig9(t, i) = marcumq(ncp_d_NO, sqrt(2*xi_NO/N0), 3);
    end
end
disp('Calculations Complete!');
%% ========================================================================
%  4. PLOTTING
%  ========================================================================

% --- FIGURE 8 ---
figure('Name', 'Physically Correct: PoD vs SNR', 'Position', [100, 100, 700, 500]);
plot(P_S_dB, PoD_NO_fig8(1,:), 'k-s', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
plot(P_S_dB, PoD_RS_fig8(1,:), 'r-s', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(P_S_dB, PoD_NO_fig8(2,:), 'k-o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(P_S_dB, PoD_RS_fig8(2,:), 'r-o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlim([10, 80]); ylim([0, 1]); yticks(0:0.25:1.00);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('SNR (dB)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability of Detection (PoD)', 'FontSize', 14, 'FontWeight', 'bold');
title('PoD vs BS SNR', 'FontSize', 14);
legend('NOMA T1', 'RSMA T1', 'NOMA T2 ', 'RSMA T2', ...
       'Location', 'southeast', 'FontSize', 11);

% --- FIGURE 9 ---
figure('Name', 'Physically Correct: ROC Curves', 'Position', [820, 100, 700, 500]);
plot(PoFA_vec, PoD_NO_fig9(1,:), 'k-s', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k'); hold on;
plot(PoFA_vec, PoD_RS_fig9(1,:), 'r-s', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(PoFA_vec, PoD_NO_fig9(2,:), 'k-o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(PoFA_vec, PoD_RS_fig9(2,:), 'r-o', 'LineWidth', 1.8, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlim([0, 1]); ylim([0, 1]); xticks(0:0.25:1.00); yticks(0:0.25:1.00);
grid on; set(gca, 'FontSize', 12, 'LineWidth', 1.2);
xlabel('Probability of False Alarm (PoFA)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability of Detection (PoD)', 'FontSize', 14, 'FontWeight', 'bold');
title('PoD of the BS vs PoFA', 'FontSize', 14);
legend('NOMA T1', 'RSMA T1', ...
       'NOMA T2', 'RSMA T2', ...
       'Location', 'southeast', 'FontSize', 11);