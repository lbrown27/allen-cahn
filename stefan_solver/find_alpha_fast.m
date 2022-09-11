function alpha = find_alpha_fast(k_l, k_s,L,T_l,T_s, rho_water, cp_water)
k_s = k_s / (rho_water * cp_water);
k_l = k_l / (rho_water * cp_water);
L = L / cp_water;
T_l = T_l - 273.15;
T_s = T_s - 273.15;
fcn = @(alpha) 1/(sqrt(pi) * L)*(sqrt(k_s)*T_s * exp(-alpha^2/k_s)/erfc(alpha/sqrt(k_s)) + ...
    sqrt(k_l)*T_l * exp(-alpha^2/k_l)/(2-erfc(alpha/sqrt(k_l))))-alpha;

alpha = fzero(fcn, -.0001);

end