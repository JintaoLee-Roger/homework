clear;clc;close all;
addpath("lib/")

load("ifk.mat");
x = (0.5:1:19.5)*0.05;
G = zeros(20,20);
for i=1:20
    for j=1:20
        G(i,j)=x(j)*exp(-x(j)*x(i))*0.05;
    end
end 

%% c. Zeroth-order Tikhonov regularization
L = eye(20);

%% L-curve criterion
% Calculate the svd
[U,S,V]=svd(G);

% calculate and plot the L-curve, and find its corner.
s=diag(S);
[rho,eta,reg_param]=l_curve_tikh_svd(U,s,d,1000);

% estimate the L-curve corner in log-log space using a spline 
% fit to find where maximum negative curvature occurs
[alpha1,ireg_corner] = l_curve_corner(rho,eta,reg_param);
rho_corner=rho(ireg_corner);
eta_corner=eta(ireg_corner);

disp(['alpha from L curve corner is ', num2str(alpha1)]);

% model
m_tikh=(G'*G+alpha1^2*L)\G'*d;

% plot L curve and add the corner marker
figure(1)
clf
loglog(rho,eta);
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Norm ||m||_2')
axis tight
hold on
H=loglog(rho_corner,eta_corner,'ko');
set(H,'markersize',8)
title("L-curve criterion (Zeroth-Order)")
hold off

% plot the L curve predicted model
figure(2)
clf
plotconst(m_tikh, 0,1);
xlabel("x");
ylabel("m(x)");
title("Model by L-curve criterion (zeroth-order)")
ylim([0.1, 1.1])

%% the discrepancy principle
discrep=(1e-4)*sqrt(20);
alpha_disc=interp1(rho,reg_param,discrep);
disp(['alpha from the discrepancy principle is ', num2str(alpha_disc)])

% get the model and residual
m_disc=(G'*G+alpha_disc^2*L)\G'*d;
r_spike_disc=norm(G*m_disc-d);

% plot the discrepancy principle predicted model
figure(3)
clf
plotconst(m_disc,0,1);
xlabel('x');
ylabel('m(x)');
title("Model by the discrepancy principle (zeroth-order)")
ylim([0.1 1.1])

%% GCV 
% generalized svd
[U0,V0,X0,Lam0,MU0]=gsvd(G,L);
lam=sqrt(diag(Lam0'*Lam0));
mu=sqrt(diag(MU0'*MU0));
p=rank(L);
sm0=[lam(1:p),mu(1:p)];

% get the gcv values varying alpha
[alpha0,g0,reg_param0]=gcval(U0,sm0,d,1000);
[ming0,~] = min(g0);

figure(4)
clf
loglog(reg_param0,g0,'-k');
xlabel('\alpha');
ylabel('g(\alpha)');
ax=axis;
hold on
H=loglog(alpha0,ming0,'ok',[alpha0,alpha0],[ming0/1000,ming0],'-k');
set(H,'markersize',8);
title("GCV (zeroth-order)")
hold off
axis(ax);

m0=(G'*G+alpha0^2*(L'*L))\G'*d;
figure(5)
clf
plotconst(m0,0,1);
xlabel('x')
ylabel('m(x)')
title("Model by GCV (zeroth-order)")

%% d. First-Order Tikhonov regularization 
% Get roughening matrix and generalized svd (using cgsvd.m)
L1 = get_l_rough(20,1);
[U1,V1,X1,Lam1,MU1]=gsvd(G,L1);
lam=sqrt(diag(Lam1'*Lam1));
mu=sqrt(diag(MU1'*MU1));
p=rank(L1);
sm1=[lam(1:p),mu(1:p)];

%% GCV
% get the gcv values varying alpha
[alpha1,g1,reg_param1]=gcval(U1,sm1,d,1000);
[ming1,~] = min(g1);

figure(6)
clf
loglog(reg_param1,g1,'-k');
xlabel('\alpha');
ylabel('g(\alpha)');
ax=axis;
hold on
H=loglog(alpha1,ming1,'ok',[alpha1,alpha1],[ming1/1000,ming1],'-k');
set(H,'markersize',8);
title("GCV (First-order)")
hold off
axis(ax);

m1=(G'*G+alpha1^2*(L1'*L1))\G'*d;
figure(7)
clf
plotconst(m1,0,1);
xlabel('x')
ylabel('m(x)')
title("Model by GCV (First-order)")

%% L-curve criterion
[rho1, eta1, reg_param1, m1] = l_curve_tikh_gsvd(U1, d, X1, ...
    Lam1, MU1, G, L1, 1000, 1e-9, 1e-1);

% estimate the L-curve corner in log-log space using a spline 
% fit to find where maximum negative curvature occurs
[alpha1,ireg_corner1] = l_curve_corner(rho1,eta1,reg_param1);
rho_corner1=rho1(ireg_corner1);
eta_corner1=eta1(ireg_corner1);

disp(['alpha from L curve corner is ', num2str(alpha1)]);

% model
m_tikh1=(G'*G+alpha1^2*(L1'*L1))\G'*d;

% plot L curve and add the corner marker
figure(8)
clf
loglog(rho1,eta1);
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Norm ||L_1m||_2')
axis tight
hold on
H=loglog(rho_corner1,eta_corner1,'ko');
set(H,'markersize',8)
title("The L-curve criterion (First-order)")
hold off

% plot the L curve predicted model
figure(9)
clf
plotconst(m_tikh1, 0, 1);
xlabel("x");
ylabel("m(x)");
title("Model by the L-curve criterion (First-order)")
ylim([0.1, 1.1])

%% the discrepancy principle
discrep=(1e-4)*sqrt(20);
alpha_disc1=interp1(rho1,reg_param1,discrep);
disp(['alpha from the discrepancy principle is ', num2str(alpha_disc1)])

% get the model and residual
m_disc1=(G'*G+alpha_disc1^2*(L1'*L1))\G'*d;

% plot the discrepancy principle predicted model
figure(10)
clf
plotconst(m_disc1,0,1);
xlabel('x');
ylabel('m(x)');
title("Model by the discrepancy principle (First-order)")
ylim([0.1 1.1])









