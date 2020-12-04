% Copyright
% Author: Jintao Lee
% Date: 11-19-2020
%
% Discription
% Inversion homework
% Chapter 3: Rank Deficiency and Ill-conditioning
% Exercise 4 in page 88

clear;clc;close all;
load rowscan.mat;
load colscan.mat;
load diag1scan.mat; % SW -> NE
load diag2scan.mat; % NW -> SE

%% Use the row and column scans only
G1 = zeros(32, 256);
d1 = zeros(32, 1);
for i = 1:1:32
    temp = zeros(16, 16);
    if i < 17
        temp(i, :) = 1;
    else
        temp(:, i-16) = 1;
    end   
    G1(i, :) = reshape(temp', 1, 256);
end
d1(1:16, :) = rowscan(:,5);
d1(17:32, :) = colscan(:,5);

rank1 = rank(G1);
model1_null_space = null(G1);
data1_null_space = null(G1');
s1 = svd(G1);

m1 = pinv(G1, 0.1) * d1;

% model resolution matrix and data space resolution matrix
R_m1 = pinv(G1, 0.1) * G1;
R_d1 = G1 * pinv(G1, 0.1);


%% Use the complete data set
G2 = zeros(94, 256);
d2 = zeros(94, 1);
d2(1:32, :) = d1;
G2(1:32, :) = G1;
for i = 1:1:31
    if i < 17
        temp = diag(sqrt(2).*ones(1,i), -16+i);
    else
        temp = diag(sqrt(2).*ones(1,32-i), -16+i);
    end
    G2(32+i, :) = reshape(flipud(temp)', 1, 256);
    G2(32+i+31, :) = reshape(temp', 1, 256);
end
d2(33:63, :) = diag1scan(:, 5);
d2(64:94, :) = diag2scan(:, 5);

rank2 = rank(G2);
model2_null_space = null(G2);
data2_null_space = null(G2');
s2 = svd(G2);

m2 = pinv(G2, 0.1) * d2;

% model resolution matrix and data space resolution matrix
R_m2 = pinv(G2, 0.1) * G2;
R_d2 = G2 * pinv(G2, 0.1);


%% the 4.f problem
d_11 = G1 * m1;
d_12 = G1 * (m1 + model1_null_space(:, 1));
res1 = d_12 - d_11;

d_21 = G2 * m2;
d_22 = G2 * (m2 + model2_null_space(:, 1));
res2 = d_22 - d_21;


%%
m1 = reshape(m1, 16,16)';
m2 = reshape(m2, 16,16)';

%% plot

figure(1)
imagesc(m1)
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])
colorbar;
temp1=caxis;
title("Model-only row and column data")
xlabel("j")
ylabel("i")

figure(2)
imagesc(m2)
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])
colorbar;
caxis(temp1)
title("Model-all data")
xlabel("j")
ylabel("i")


%%
figure(3)
imagesc(R_m2)
xticks(1:16:256)
yticks(1:16:256)
temp2=caxis;
colorbar;
title("Model resolution 2")
xlabel("j")
ylabel("i")

figure(4)
imagesc(R_m1)
xticks(1:16:256)
yticks(1:16:256)
caxis(temp2)
colorbar;
title("Model resolution 1")
xlabel("j")
ylabel("i")

%%
figure(5)
imagesc(reshape(diag(R_m2), 16,16)')
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])
caxis([0, 1])
colorbar;
title("Model resolution--all data")
xlabel("j")
ylabel("i")

figure(6)
imagesc(reshape(diag(R_m1), 16,16)')
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])
caxis([0, 1])
colorbar;
title("Model resolution--only row and column")
xlabel("j")
ylabel("i")

%% 
figure(7)
imagesc(model1_null_space)
xticks([1, 25, 50, 75, 100, 125, 150, 175, 200, 225])
yticks([1, 32, 64, 96, 128, 160, 192, 224, 256])
% caxis([0, 1])
colorbar;
title("Model null space--only row and column scans")
xlabel("j")
ylabel("i")

figure(8)
imagesc(model2_null_space)
xticks([1, 25, 50, 75, 100, 125, 150])
yticks([1, 32, 64, 96, 128, 160, 192, 224, 256])
% caxis([0, 1])
colorbar;
title("Model null space--all data set")
xlabel("j")
ylabel("i")

%%
figure(9)
sum1 = sum(model2_null_space);
max(sum1)
min(sum1)
imagesc(reshape(model2_null_space(:, 1), 16,16)')
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])

colorbar;
title("One model null space--all data")
xlabel("j")
ylabel("i")

figure(10)
sum2 = sum(model1_null_space);
max(sum2)
min(sum2)
imagesc(reshape(model1_null_space(:, 1), 16,16)')
xticks([1, 4, 8, 12, 16])
yticks([1, 4, 8, 12, 16])

colorbar;
title("one mdoel null space--only row and column")
xlabel("j")
ylabel("i")

%% 
figure(11)
imagesc(data1_null_space)
yticks([1, 8, 16, 24, 32])
% caxis([0, 1])
colorbar;
title("data null space--only row and column scans")
xlabel("j")
ylabel("i")

figure(12)
imagesc(data2_null_space)
xticks([1,3,5,7])
yticks([1, 32, 64, 94])
% caxis([0, 1])
colorbar;
title("data null space--all data set")
xlabel("j")
ylabel("i")

