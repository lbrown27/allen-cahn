clear all;
close all;
clc;

num_pts = 100000;

x = linspace(0,10,num_pts)';
dx = x(2)-x(1);
y = x.^3;

stag_x = stagger(x,num_pts);

deriv_num = div_stag(stag_x.^3,num_pts,dx);
deriv_num_2 = div_stag(stagger(x.^3,num_pts),num_pts,dx);

deriv_ana = 3*x.^2;

diff1= deriv_num - deriv_ana
diff2 = deriv_num_2-deriv_ana
max(diff1(2:end-1))