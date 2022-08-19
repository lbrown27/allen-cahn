function ans = RHS_PE(rho_f, rho_n, rho_p, pc, u_star)
% PARAMS: rho_f - rho future
%         rho_n - rho now
%         rho_p - rho past
for i = 3:pc.N - 1
    i_plus = i;
    i_minus = i - 1; % indices for flux terms
    
    ans = 3 * (3 * rho_f - 4 * rho_n + rho_p)/(2 * pc.dt) ...
        + 3 * (u_star(i_plus) * (rho_f(i+1) + rho_n(i)) - u_star(i_minus) * (rho_f(i-1) + rho_f(i)));
    ans = ans / (2 * pc.dx);
end
end