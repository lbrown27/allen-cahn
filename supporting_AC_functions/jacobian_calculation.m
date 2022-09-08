function dF = jacobian_calculation(u_n,rho_new,eta_new,pc,c_new)
%% This function calculates the jacobian so that can use in the crank-nicolson
%  method for implicit time advancement

du = 1e-16;
current_RU = RU_n_func(rho_new, eta_new, pc, c_new,u_n);
dF = zeros(pc.N+2,pc.N+2);
for i = 1:pc.N+2
    u_probe = u_n;
    u_probe(i) = u_probe(i) + du;
dF(:,i) = (RU_n_func(rho_new, eta_new, pc, c_new, u_probe) - current_RU)/du; % one column of the jacobian matrix.
end
end