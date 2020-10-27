clear;clc;
%addpath('/home/jtli/code/matlab/plot2svg/src')
addpath('/home/jtli/code/matlab/export_fig')
G = 0.2 .* tril(ones(100));
z = transpose(0.2:0.2:20);
d = log(1+z/25)/40;

m_pre = G \ d;
m_true = 1 ./ (1000 + 40*z);
err = abs(m_pre-m_true);

d_noise = d+5*10^(-5).*randn(100, 1);
m_pre_noise = G \ d_noise;
err_add_noise = abs(m_pre_noise-m_true);

figure(1);
plot(z, m_pre)
hold on
plot(z, m_true)
hold on 
plot(z, m_pre_noise)
xlabel('z/m')
ylabel('t/s')
title('N = 100')
legend('m_{true}', 'm_{pre}', 'm_{pre\_add\_noise}')
 %plot2svg('fig1.svg')
export_fig fig1.pdf -r300 -a2 -opengl

figure(2)
plot(z, err)
hold on
plot(z, err_add_noise)
xlabel('z/m')
ylabel('error/s')
title('N = 100')
legend('err_{m_{true} - m_{pre}}', 'err_{m_{true} - m_{pre\_add\_noise}}')
export_fig fig2.pdf -r300 -a2 -opengl
%plot2svg('fig2.svg')

z2 = transpose(linspace(5,20,4));
G2 = 5 .* tril(ones(4));
d2 = log(1+z2/25)/40;
 
m2_pre = G2 \ d2;
m2_true = 1 ./ (1000 + 40*z2);
d2_noise = d2+5*10^(-5).*randn(4, 1);
m2_pre_noise = G2 \ d2_noise;
err2 = abs(m2_pre-m2_true);
err2_add_noise = abs(m2_pre_noise-m2_true);


figure(3);
plot(z2, m2_pre)
hold on
plot(z2, m2_true)
hold on 
plot(z2, m2_pre_noise)
xlabel('z/m')
ylabel('t/s')
title('N = 4')
legend('m_{true}', 'm_{pre}', 'm_{pre\_add\_noise}')
export_fig fig3.pdf -r300 -a2 -opengl
%plot2svg('fig3.svg')

figure(4)
plot(z2, err2)
hold on
plot(z2, err2_add_noise)
xlabel('z/m')
ylabel('error/s')
title('N = 4')
legend('err_{m_{true} - m_{pre}}', 'err_{m_{true} - m_{pre\_add\_noise}}')
export_fig fig4.pdf -r300 -a2 -opengl
%plot2svg('fig4.svg')

figure(5)
plot(z, err)
hold on 
plot(z2, err2)
xlabel('z/m')
ylabel('error/s')
title('error without noise')
legend('N=100', 'N=4')
export_fig fig5.pdf -r300 -a2 -opengl

figure(5)
plot(z, err)
hold on 
plot(z2, err2)
xlabel('z/m')
ylabel('error/s')
title('error with noise')
legend('N=100', 'N=4')
export_fig fig6.pdf -r300 -a2 -opengl
