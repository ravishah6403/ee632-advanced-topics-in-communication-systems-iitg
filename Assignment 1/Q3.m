clear; clc;

N = 10000;
tot_var = 132;

X = sqrt(tot_var / 2) * (randn(N, 1) + 1i*randn(N, 1));

amplitude = abs(X);
phase = angle(X);

r_vals = linspace(0, max(amplitude), 500);
r_pdf = (r_vals / (tot_var / 2)) .* exp(-(r_vals.^2)/tot_var);

fig = figure;

subplot(2,1,1);

histogram(amplitude, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;
plot(r_vals, r_pdf, 'r', 'LineWidth', 2);

title('Amplitude Distribution (Rayleigh)');
xlabel('Magnitude'); ylabel('PDF');
legend('Simulated', 'Analytical');
grid on;

theta_vals = linspace(-pi, pi, 500);
theta_pdf = ones(size(theta_vals)) * (1 / (2*pi));

subplot(2,1,2);

histogram(phase, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;
plot(theta_vals, theta_pdf,'r', 'LineWidth', 2);

title('Phase Distribution (Uniform)');
xlabel('Phase'); ylabel('PDF');
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
legend('Simulated', 'Analytical');
grid on;

saveas(fig, 'Q3_figure', 'fig');