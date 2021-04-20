clear;close all;clc;

%%% parameters
E = 70; % GPa, Young's modulus
poisson = 0.25; % Poisson's ratios
rou_basalt = 2700; % kg/m^3, basalt density
rou_crust = 2900; % kg/m^3, crust density
th_crust = [30, 50, 70]; % km, crust thickness
v = 300; % km/Ma, Basalts spread velocity
h = 2; % km, basalts accumulated height
T = 2; % Ma, Total time of volcanic eruption

%%% 
% D = E*T^3/(12*(1-poisson^2));
rou_prime = rou_mantle + E*th_crust/
x = [150, 300, 450];
t = 0.1:0.1:T;
steps = 0.005:0.005:2;

