function T = T_init(c,x,pc)
%c = (tanh((x - pc.l / 2)/(sqrt(2) * pc.ksi_c)) - 1)/2;
T = zeros(pc.N + 2,1);
%if x(i) < pc.ksi_c
% c(i) = 1 / pc.ksi_c * x(i) - 1;
thick_factor = .1;
width = 2.5 *pc.ksi_c;
T = (pc.wall_T +pc.T_M)/2+(pc.T_M - pc.wall_T)/2*tanh((x-(pc.x_init - width/2))/(thick_factor*pc.ksi_c))+(pc.init_T +pc.T_M)/2  - pc.T_M + (pc.init_T - pc.T_M)/2 * tanh((x-(pc.x_init + width/2))/(thick_factor*pc.ksi_c));
%end
%T = (pc.wall_T + pc.init_T)/2 + (pc.init_T - pc.wall_T)/2 * tanh((x-pc.x_init)/(2*pc.ksi_c));

end