function res = RHS_PE(rho_new,pc, u_n,c_new,T_new,eta_new, u_star,rho_n, rho_old,x_coll)
% PARAMS: rho_f - rho future
%         rho_n - rho now
%         rho_p - rho past
persistent n
n=1;
if isempty(n)
    n = 1;
    fp = (rho_new - rho_n)/pc.dt;
else
    fp = rhs_ac_wrapper(c_new,T_new,u_n,eta_new,rho_new, x_coll, pc,1)*(pc.rho_water - pc.rho_ice);
end
%fp = (3 * rho_new - 4 * rho_n + rho_old)/(2 * pc.dt);
%fp = (rho_new - rho_n)/pc.dt;
sp = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1; % indices for flux terms
    sp(i) = (u_star(i_plus) * (rho_new(i+1) + rho_new(i)) - u_star(i_minus) * (rho_new(i-1) + rho_new(i)))/(2 * pc.dx);
end
res = fp + sp;

res = res / (pc.dt);
end