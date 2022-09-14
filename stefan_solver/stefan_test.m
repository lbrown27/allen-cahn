%% Stefan Plotter
clear all;
close all;
clc;

alpha_guess = .13;
k_l = .005;
k_s = 1;
L = .53;
T_l = .53;
T_s = -.01;
start_loc = .2;


alpha_guess = -.00000014;
therm_cond_ice = 2.25;
therm_cond_water = .5918;
density_ice = 898;
rho_water = 998;
heat_capacity_ice = 2018;
cp_water = 4200;

%% Experimental function solve



k_l = therm_cond_water;
k_s = therm_cond_ice;
L = 334000;
T_l = .2;
T_s = -5;
L_domain = 100;
x = 0:.01:L_domain;

start_loc = L_domain /2;

alpha = find_alpha_fast(k_l, k_s,L,T_l,T_s, rho_water, cp_water);

%alpha2 = find_alpha(k_l, k_s,L,T_l,T_s,alpha_guess);
t = 0;
T = stefan_temp_field(T_l, k_l,alpha, start_loc,t,T_s, k_s,x,L_domain);

plot(x,T);
hold on;
t_next = 20;
T_next = stefan_temp_field(T_l, k_l,alpha, start_loc,t_next,T_s, k_s,x,L_domain);
int_loc = interface_location(start_loc, alpha,t_next,L_domain)
int_dx = int_loc - start_loc;
figure(1);
plot(x,T_next);
drawnow;