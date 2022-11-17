function T = T_initialization(c,rho, cp,x_coll,pc)
sharp_terms = (x_coll > pc.x_init) .* (pc.rho_water * pc.L+ (pc.cp_water * (pc.init_T - pc.T_M)*pc.rho_water)) + ...
    (x_coll <= pc.x_init) .* (pc.cp_ice*(pc.wall_T - pc.T_M));
T = (sharp_terms - c .*pc.L.*rho+cp.*rho.*pc.T_M)./(cp.*rho);
T(1602) = .9637+(1.663-.9637)*.75;
T(1601) = .9637 + (1.663-.9637)*.25;

T_arezoo = stefan_temp_field(pc.init_T, pc.thermal_diff_water,pc.alpha, pc.x_init,.01,pc.wall_T, pc.thermal_diff_ice,x_coll,pc.l,pc.T_M);
%T = T_arezoo;
end