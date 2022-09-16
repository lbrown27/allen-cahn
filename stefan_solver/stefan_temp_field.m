function T_ana = stefan_temp_field(T_l, k_l,alpha, start_loc,t,T_s, k_s,x,L_domain,T_M)
%% Calculates the analytical temperature field for 1D Stefan problem, such
% that the left side is ice and the right side is water.
% PARAMS:
% T_l: Temperature of the liquid
% k_l: Thermal diffusivity of liquid
% alpha: Freezing speed parameter of problem, see the find_alpha_fast
% function
% start_loc: Starting location of the interface
% t: time from initial t0 at which you want the temp field
% T_s: Initial temperature of the solid
% k_s: Thermal diffusivity of the solid
% x: domain discretization vector
% L_domain: Length of the domain
x_flipped = L_domain - x;
start_loc = L_domain - interface_location(start_loc, alpha,0,L_domain);
%start_loc = L_domain - start_loc;

T_l = T_l - T_M;
T_s = T_s - T_M;
t_1 = -T_l * erfc(alpha / sqrt(k_l)) / (2-erfc(alpha/sqrt(k_l)));
t_2 = T_l * erfc((x_flipped - start_loc)/(2*sqrt(k_l * t)))/(2 - erfc(alpha/sqrt(k_l)));
t_3 = T_s;
t_4 = -T_s * erfc((x_flipped - start_loc)/(2*sqrt(k_s * t)))/(erfc(alpha/sqrt(k_s)));

T_ana = (t_1 + t_2).*(x_flipped < start_loc) + (t_3 + t_4).*(x_flipped >= start_loc);
T_ana = T_ana + T_M;


end