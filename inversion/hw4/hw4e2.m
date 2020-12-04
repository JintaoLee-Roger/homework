% Copyright
% Author: Jintao Li
% Date: 11-19-2020
% 
% Discription
% Inversion homework
% Chapter 3: Rank Deficiency and Ill-conditioning
% Exercise 2 in page 87

clear; clc; close all;

m_true = [-1, 1, -1, 1, -1, 1, -1, 1, -1]';
G = [1,0,0,1,0,0,1,0,0; 0,1,0,0,1,0,0,1,0; 0,0,1,0,0,1,0,0,1; ...
     1,1,1,0,0,0,0,0,0; 0,0,0,1,1,1,0,0,0; 0,0,0,0,0,0,1,1,1; ...
     sqrt(2),0,0,0,sqrt(2),0,0,0,sqrt(2); ...
     0,0,0,0,0,0,0,0,sqrt(2)];
 
% Forward, i.e. generate synthetic data
d = G * m_true;

% generalized inverse solution
diag_s = svd(G);
m = pinv(G, 0.001) * d;

% comparion between m and m_true
m = reshape(m, 3, 3)';
m_true = reshape(m_true, 3, 3)';

% model resolution matrix
R_m = pinv(G) * G;

% plot
figure(1)
imagesc(m)
xticks([1, 2, 3])
yticks([1, 2, 3])
colorbar;
temp1=caxis;
title("Generalized inverse solution")
xlabel("j")
ylabel("i")

figure(2)
imagesc(m_true)
xticks([1, 2, 3])
yticks([1, 2, 3])
caxis(temp1)
colorbar;
title("The true model")
xlabel("j")
ylabel("i")

figure(3)
imagesc(reshape(diag(R_m), 3,3))
xticks([1,2,3])
yticks([1,2,3])
colorbar;
title("Model resolution matrix (the diagonal elements)")


 