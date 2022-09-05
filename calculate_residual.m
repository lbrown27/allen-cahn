function residual = calculate_residual(c,T,u_new,eta,rho, pc)
d_rhou_dx = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
   d_rhou_dx(i) = ((rho(i+1) + rho(i)) * u_new(i) - (rho(i) + rho(i-1)) * u_new(i - 1))/(2*pc.dx);
end
ac_d_rho = rhs_ac(c,T,u_new,eta,rho, pc) * (pc.rho_water - pc.rho_ice);
residual = ac_d_rho + d_rhou_dx;
residual(1) = 0;
residual(end) = 0;
end