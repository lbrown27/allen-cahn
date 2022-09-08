function res = RHS_PE(rho_new,pc, u_n,c_new,T_new,eta_new, u_star,rho_n, rho_old)
% PARAMS: rho_f - rho future
%         rho_n - rho now
%         rho_p - rho past
    fp = rhs_ac(c_new,T_new,u_n,eta_new,rho_new, pc)*(pc.rho_water - pc.rho_ice);
    fp = (3 * rho_new - 4 * rho_n + rho_old)/(2 * pc.dt);
    fp = (rho_new - rho_n)/pc.dt;
    sp = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1; % indices for flux terms
    sp(i) = (u_star(i_plus) * (rho_new(i+1) + rho_new(i)) - u_star(i_minus) * (rho_new(i-1) + rho_new(i)))/(2 * pc.dx);
end
    res = fp + sp;

res = res / (pc.dt);
end