clear;clc;close all;

x = (-0.95:0.1:0.95)';
y = x;

G = zeros(20,20);
for i=1:20
   G(:, i) = x.^(i-1);
end

m = G \ y;

m_true = zeros(20, 1);
m_true(2) = 1;

y_pre = G*m;

figure(1)
plot(x, y, "ro")
hold on;
plot(x, y_pre, "-b")
xlabel("x")
ylabel("y")
legend(["the data", "the fitted model"], "Location", "southeast")

figure(2)

r = m - m_true;

a = 0:1:19;
plot(a, r, "-o")
xlabel("a")
ylabel("the difference")
fig2svg("hw3e52.svg")
