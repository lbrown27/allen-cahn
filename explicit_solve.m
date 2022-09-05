function u_star = explicit_solve(rho_n, rho_new, eta_n, pc, c_n, u_n,rho_old,eta_old, c_old,u_old)
    %% performs the explicit Adams-Bashforth method as described in the zheng icing paper.
       RU_n =  RU_n_func(rho_n, eta_n, pc,c_n,u_n);
    RU_n_old = RU_n_func(rho_old, eta_old, pc, c_old, u_old);
    
    u_star = (4 .* pc.dt .* RU_n  -  2 .* pc.dt .* RU_n_old + 4 .* rho_n .* u_n - rho_old .* u_old)./(3.* rho_new);


end
