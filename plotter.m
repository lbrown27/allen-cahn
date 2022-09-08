clear all;
close all;
clc;

num_gridpoints = [200 400 800 1600];
error = [1.7955e-5 4.4703e-6 1.11e-6 2.7852e-7];
figure();
loglog(num_gridpoints, error,'.');

n_log = log(num_gridpoints);
error_log = log(error);
hold on;
L = error(1)/num_gridpoints(1)^2;
loglog(num_gridpoints,1./(num_gridpoints.^2/L),'.');