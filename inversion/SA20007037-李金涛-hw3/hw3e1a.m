% copyright
% author: Jintao Li

clear; clc; close all;
addpath("fig2svg");

E = 0;
Sigma = 0.1;

x = [6.0000 10.1333 14.2667 18.4000 22.5333 26.6667]';
t = [3.4935 4.2853 5.1374 5.8181 6.8632 8.1841]';

G = [ones(6, 1), x];

m_L2 = inv(G' * G) * G' * t;

r = t - G * m_L2;

%% 1-a
figure(1)
plot(x, t, "ro");
hold on;
t_m = G * m_L2;
plot(x, t_m, "-b");
xlabel("Distances/km");
ylabel("The first arrival times/s")
legend(["data", "the fitted model"], 'Location', 'southeast')

figure(2)
plot(x, r, "-o");
xlabel("Distances/km");
ylabel("The differece of times/s")
title("The residual")


%% 1-b
C = Sigma^2 .* inv(G' * G);
rho = C(1,2) ./ sqrt(C(1,1) .* C(2,2));
correlation_matrix = [1 rho; rho 1];

%% 1-c
DeltaS = 5.99;
plot_ellipse(DeltaS, C, m_L2)
xlabel("t_0 / s")
ylabel("s_2 / km")
% fig2svg('hw3e1c.svg')
lambda = eig(inv(C));
intv = sqrt(DeltaS) .* lambda .^ 0.5;


%% 1-d
chiValue = sum(r.^2 ./ Sigma^2);
chi = chi2cdf(chiValue, 4);
p = chi2cdf(chiValue, 4, 'upper');

%% 1-e
chi_sim = zeros(1, 1000);
p_sim = zeros(1, 1000);
chiv = zeros(1, 1000);
for i = 1:1:1000
    noise = 0.1*randn(6, 1);
    t_noise = t + noise;
    m_noise = inv(G' * G) * G' * t_noise;
    r = t_noise - G * m_noise;
    chiv(i) = sum(r .^2 ./ Sigma^2);
    chi_sim(i) = chi2cdf(chiv(i), 4);
    p_sim(i) = chi2cdf(chiv(i), 4,'upper');
end

n_bins = 0:1:30;y1=chi2pdf(n_bins,4);
histogram(chiv);
hold on
plot(n_bins,y1*1000)
% fig2svg("hw3e1e.svg")


%% 1-f
histogram(p_sim-p);
xlabel("difference (p_M - p) ")
ylabel("times")
% fig2svg("hw3e1f.svg")






