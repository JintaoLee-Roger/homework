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

figure(1)
plot(x, t, "ro");
hold on;
t_m = G * m_L2;
plot(x, t_m, "-b");
hold on;

m1 = IRLS(t, G);


plot(x, G * m1, "--g");

xlabel("Distances/km");
ylabel("The first arrival times/s")
legend(["data", "L2 solution", "IRLS"], 'Location', 'southeast')
fig2svg("hw3e1g.svg")

%% 1-h
q = 1000;
m_all = zeros(q, 2);

for i = 1:1:q
    noise = Sigma .* randn(6, 1);
    t_noise = t + noise;
    m_all(i, :) = (IRLS(t_noise, G))';
end

A = m_all - mean(m_all);

CovML1 = A' * A ./ q;

conf = 1.96 .* diag(CovML1) .^ 0.5;




