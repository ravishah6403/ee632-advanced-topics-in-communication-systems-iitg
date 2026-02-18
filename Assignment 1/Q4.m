clear; clc;

N = 10000;
mu_real = 4;
mu_imag = 8;
var = 16;

X = (mu_real + sqrt(var)*randn(N, 1)) + 1i*(mu_imag + sqrt(var)*randn(N, 1));

amplitude = abs(X);
phase = angle(X);

s = sqrt(mu_real^2 + mu_imag^2);


r_vals = linspace(0, max(amplitude), 500);
r_pdf = (r_vals./var) .* exp(-(s^2 + r_vals.^2)/(2*var)) .* besseli(0, (s * r_vals)/var);

fig = figure;

subplot(2,1,1);

histogram(amplitude, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;
plot(r_vals, r_pdf, 'r', 'LineWidth', 2);

title('Amplitude Distribution (Rician)');
xlabel('Magnitude'); ylabel('PDF');
legend('Simulated', 'Analytical');
grid on;

theta_vals = linspace(0, 2*pi, 500);
m = mu_real*cos(theta_vals) + mu_imag*sin(theta_vals);
theta_pdf = (1/(2*pi)) * exp(-s^2/(2*var)) + (m./sqrt(2*pi*var)) .* exp((m.^2 - s^2)/(2*var)) .* qfunc(-m./sqrt(var));

subplot(2,1,2);

histogram(phase, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;
plot(theta_vals, theta_pdf,'r', 'LineWidth', 2);

title('Phase Distribution');
xlabel('Phase'); ylabel('PDF');
xticks([-pi -pi/2 0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
legend('Simulated', 'Analytical');
grid on;

saveas(fig, 'Q4_figure', 'fig');