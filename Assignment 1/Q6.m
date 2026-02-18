clear; clc;

N = 10000;
n = 10;
var = 4;

X = sqrt(var/2)*(randn(n, N) + 1i*randn(n, N));

Y = sum(abs(X).^2, 1);

fig = figure;

histogram(Y, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;

c_vals = linspace(0, max(Y), 500);
c_pdf = (1 ./ ((var).^(n) .* gamma(n))) .* c_vals.^((n)-1) .* exp(-c_vals ./ (var));

plot(c_vals, c_pdf, 'r', 'LineWidth', 2);

title(sprintf('Chi-Squared Distribution (N = %d, n = %d)', N, n));
xlabel('y'); ylabel('PDF');
legend('Simulated', 'Analytical');
grid on;

saveas(fig, 'Q6_figure', 'fig');