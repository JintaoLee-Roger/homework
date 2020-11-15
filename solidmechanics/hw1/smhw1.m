% copyright
% Author: Jintao Li
% Date: 2020-11-15

clear;clc;close all;

U = 0.05 / ((356*24+6)*60*60); % m/s
mu = 1e14 ; % Mpa/s
% mu = 1e20; % pa/s

%% calculate consistant A,B,C,D
% A1,B1,C1,D1 of upper plate
c = [1 0 0 0; ...
    0 1 1 0; ...
    0 0 0.25 -(pi/6 + sqrt(3)/4); ...
    0 0 (pi/6 - sqrt(3)/4), -0.25] \ ...
    [0;0;sqrt(3) * U /2; 0.5*U];
A1 = c(1); B1 = c(2); C1 = c(3); D1 = c(4);

% A2,B2,C2,D2 of down plate
c = [1,0,pi,0; ...
    0,1,1,pi; ...
    0,0,0.25,-(sqrt(3)/4 - 5*pi/6); ...
    0,0,-(sqrt(3)/4 + 5*pi/6),-0.25] \ ...
    [0;-U;sqrt(3) * U /2; 0.5*U];
A2 = c(1); B2 = c(2); C2 = c(3); D2 = c(4);

%% calculate u, v, P
% TODO: range
x = -75:1:250;
y = 0:1:300;
[X, Y] = meshgrid(x, y);
X = X*1000;
Y = Y*1000;
[n1,n2] = size(X);
% TODO: be consistent with range
actan_yx = atan(Y./X)+[ones(n1,75)*pi, zeros(n1, n2-75)];

u_upper = -B1 - D1 .* actan_yx + ...
    (C1 .* X + D1 .* Y) .* (-X) ./ (X.^2 + Y.^2);
v_upper = A1 + C1 .* actan_yx + ...
    (C1 .* X + D1 .* Y) .* (-Y) ./ (X.^2 + Y.^2);
P_upper = -2.*mu.*((C1 .* X + D1 .* Y)) ./ (X.^2 + Y.^2);

u_down = -B2 - D2 .* actan_yx + ...
    (C2 .* X + D2 .* Y) .* (-X) ./ (X.^2 + Y.^2);
v_down = A2 + C2 .* actan_yx + ...
    (C2 .* X + D2 .* Y) .* (-Y) ./ (X.^2 + Y.^2);
P_down = -2.*mu.*((C2 .* X + D2 .* Y)) ./ (X.^2 + Y.^2);

u = u_down; 
v = v_down; 
P = P_down; 

for i = 1:1:n1
    for j = 76:1:n2
        if (sqrt(3) * Y(i, j)) < X(i, j)
            u(i, j) = u_upper(i, j);
            v(i, j) = v_upper(i, j);
            P(i, j) = P_upper(i, j);
        end
    end
end

% there is NaN number at point (0,0), i.e. u(76, 1) and so on.
% and replace it with 0

P(isnan(P))=0;
u(isnan(u))=0;
v(isnan(v))=0;

%% calculate tau
% calculate du/dx, du/dy, dv/dx, dv/y
[u_x, u_y] = gradient(u, 1000, 1000);
[v_x, v_y] = gradient(v, 1000, 1000);
[P_x, P_y] = gradient(P, 1000, 1000);

tau_xx = 2 * mu .* u_x;
tau_yy = 2 * mu .* v_y;
tau_xy = mu .* (u_y + v_x);

sigma_xx = P - tau_xx;
sigma_yy = P - tau_yy;
sigma_xy = tau_xy;

% calculate the angles and sigma1 sigma2 tau_max
alpha1 = 0.5 .* atan(2*sigma_xy ./ (sigma_xx - sigma_yy));
alpha2 = alpha1 + pi/2;
sigma1 = 0.5 .* (sigma_xx + sigma_yy) + ...
    ((0.5 .* (sigma_xx - sigma_yy)) .^2 + sigma_xy .^2).^0.5;
sigma2 = 0.5 .* (sigma_xx + sigma_yy) - ...
    ((0.5 .* (sigma_xx - sigma_yy)) .^2 + sigma_xy .^2).^0.5;

tau_max = 0.5 * abs(sigma1 - sigma2);


%% plot 
% 
figure(1)

pcolor(X/1000, Y/1000, P);
shading interp; % plot color by interpolation
set(gca,'YDir','reverse'); % reverse Y coordinate
cb=colorbar;
set(gca,'CLim',[-50,25]);
set(get(cb,'Title'),'string','MPa');
xlabel("x/km")
ylabel("h/km")
title("Velocity and pressure fields wth U=5 cm/year and theta=\pi/6")
hold on;

dx=-75:15:250;
dy=0:15:300;
[Dx,Dy]=meshgrid(dx,dy);
u_p=u(1:15:end,1:15:end);
v_p=v(1:15:end,1:15:end);
quiver(Dx,Dy,u_p,v_p); 
hold on;
quiver(40,12,5*8,0*8,'autoscale','off','color','w');
text(50,5,'5 cm/year','color','w')
hold off; 

% 
figure(2)
% tau
pcolor(X/1000, Y/1000, tau_max);
shading interp; % plot color by interpolation
set(gca,'YDir','reverse'); % reverse Y coordinate
cb=colorbar;
set(gca,'CLim',[0,30]);
set(get(cb,'Title'),'string','MPa');
xlabel("x/km")
ylabel("h/km")
title("Principal stress and maximum shear stress")
hold on;

% sigma1 and sigma2
scale = 100;
dx=-75:20:250;
dy=0:20:300;
[Dx,Dy]=meshgrid(dx,dy);

% sigma1_x sigma1_y sigma2_x sigma2_y
% Use the four arrows to represent the cross
sigma1_x = sigma1 .* cos(alpha1);
sigma1_y = sigma1 .* sin(alpha1);

sigma2_x = sigma2 .* cos(alpha2);
sigma2_y = sigma2 .* sin(alpha2);


sigma1_x_p=sigma1_x(1:20:end,1:20:end);
sigma1_y_p=sigma1_y(1:20:end,1:20:end);
q = quiver(Dx,Dy,sigma1_x_p,sigma1_y_p, 10, 'r');
q.ShowArrowHead = 'off';
hold on;

sigma1_x_p=-sigma1_x(1:20:end,1:20:end);
sigma1_y_p=-sigma1_y(1:20:end,1:20:end);
q = quiver(Dx,Dy,sigma1_x_p,sigma1_y_p, 10, 'r');
q.ShowArrowHead = 'off';
hold on;

sigma2_x_p=sigma2_x(1:20:end,1:20:end);
sigma2_y_p=sigma2_y(1:20:end,1:20:end);
q = quiver(Dx,Dy,sigma2_x_p,sigma2_y_p, 8, 'r'); 
q.ShowArrowHead = 'off';
hold on;

sigma2_x_p=-sigma2_x(1:20:end,1:20:end);
sigma2_y_p=-sigma2_y(1:20:end,1:20:end);
q = quiver(Dx,Dy,sigma2_x_p,sigma2_y_p, 8, 'r'); 
q.ShowArrowHead = 'off';
hold on;

q = quiver(200,250,sigma1(145,60)*5,0, 8, 'w');
text(175,260,'3 Mpa','color','w')
q.ShowArrowHead = 'off';
hold off; 


















