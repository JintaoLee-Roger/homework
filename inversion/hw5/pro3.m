clear;clc;close all;

load("crosswell.mat")
p = rank(G);

%% a.
[U0,S0,V0] = svd(G);
s = diag(S0);
[rho,eta,reg_param]=l_curve_tsvd(U0,s,dn);
[alpha0,ireg_corner] = l_curve_corner(rho,eta,reg_param);
rho_corner=rho(ireg_corner);
eta_corner=eta(ireg_corner);
disp(['alpha from L curve corner is ', num2str(alpha0)]);

figure(1)
clf
loglog(rho,eta);
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Norm ||m||_2')
% axis tight
% hold on
% H=loglog(rho_corner,eta_corner,'ko');
% set(H,'markersize',8)
% title("L-curve criterion curve")
% hold off

% model
m_tikh=(G'*G+alpha0^2*eye(size(G'*G)))\G'*dn;
m_tikh = reshape(m_tikh, 16,16)';

figure(2)
imagesc(m_tikh)
title("Model by TSVD")

%% b
L = eye(size(G));
% Calculate the svd
[U,S,V]=svd(G);

% calculate and plot the L-curve, and find its corner.
s=diag(S);
[rho,eta,reg_param]=l_curve_tikh_svd(U,s,dn,1000);

% estimate the L-curve corner in log-log space using a spline 
% fit to find where maximum negative curvature occurs
[alpha1,ireg_corner] = l_curve_corner(rho,eta,reg_param);
rho_corner=rho(ireg_corner);
eta_corner=eta(ireg_corner);

disp(['alpha from L curve corner is ', num2str(alpha0)]);

% model
m_tikh=(G'*G+alpha1^2*(L'*L))\G'*dn;
m_tikh = reshape(m_tikh, 16, 16)';

% plot L curve and add the corner marker
figure(3)
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
figure(4)
clf
imagesc(m_tikh)
xlabel("x");
ylabel("m(x)");
title("Model by L-curve criterion (zeroth-order)")

%% c
L2=zeros(14*14,256);
k=1;
for i=2:15
    for j=2:15
        M=zeros(16,16);
        M(i,j)=-4;
        M(i,j+1)=1;
        M(i,j-1)=1;
        M(i+1,j)=1;
        M(i-1,j)=1;
        L2(k,:)=reshape(M,256,1)';
        k=k+1;
    end
end

[U2,V2,X2,Lam2,MU2]=gsvd(G,L2);
lam=sqrt(diag(Lam2'*Lam2));
mu=sqrt(diag(MU2'*MU2));
p=rank(L2);
sm2=[lam(1:p),mu(1:p)];

[rho2, eta2, reg_param2, m2] = l_curve_tikh_gsvd(U2, dn, X2, ...
    Lam2, MU2, G, L2, 1000, 1e-9, 1e-1);

% estimate the L-curve corner in log-log space using a spline 
% fit to find where maximum negative curvature occurs
[alpha2,ireg_corner2] = l_curve_corner(rho2,eta2,reg_param2);
rho_corner2=rho2(ireg_corner2);
eta_corner2=eta2(ireg_corner2);

disp(['alpha from L curve corner is ', num2str(alpha2)]);

% model
m_tikh2=(G'*G+alpha2^2*(L2'*L2))\G'*dn;

% plot L curve and add the corner marker
figure(5)
clf
loglog(rho2,eta2);
xlabel('Residual Norm ||Gm - d||_2')
ylabel('Solution Norm ||L_2m||_2')
axis tight
hold on
H=loglog(rho_corner2,eta_corner2,'ko');
set(H,'markersize',8)
title("The L-curve criterion (Second-order)")
hold off

% plot the L curve predicted model
figure(6)
clf
imagesc(reshape(m_tikh2, 16,16)')
title("Model by the L-curve criterion (Second-order)")


