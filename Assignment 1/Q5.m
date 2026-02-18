clear; clc;

N = 10000;
n = 80;
lamda = 0.4;

X = sqrt(1/(2*lamda))*randn(2*n, N);
Y = sum(X.^2, 1);

fig = figure;
histogram(Y, 100, 'Normalization', 'pdf', 'FaceColor', [0.3 0.7 0.9], 'EdgeAlpha', 0.1);
hold on;

g_vals = linspace(0, max(Y), 500);
g_pdf = (1 ./ ((1/lamda).^n .* gamma(n))) .* g_vals.^(n-1) .* exp(-g_vals .* lamda);

plot(g_vals, g_pdf, 'r', 'LineWidth', 2);

title('Gamma Distribution');
xlabel('y'); ylabel('PDF');
legend('Simulated', 'Analytical');
grid on;

saveas(fig, 'Q5_figure', 'fig');