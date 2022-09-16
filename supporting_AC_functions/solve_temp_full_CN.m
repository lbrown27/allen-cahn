function T_new = solve_temp_full_CN(rho_cp_n, rho_cp_new,u_new,u_n,k_new, k_n,T_n,c_new, eta_new,rho_new,c_n,eta_n,rho_n,pc)

    E_new = E_matrix_shifted(rho_cp_new, u_new, k_new, pc);
    E_n = E_matrix_shifted(rho_cp_n, u_n, k_n, pc);
    
    % Crank nicolson
    % next line technically splits method but i think effect is small.
    h_prime = h_prime_func(c_new,T_n,pc.T_M);

    A_temp = sparse(diag(rho_cp_new) - .5 * pc.dt * E_new+.5*diag(pc.Mc*pc.L^2*pc.rho_water^2*h_prime*pc.dt/pc.T_M));
   % A_temp = sparse(diag(rho_cp_new) - .5 * pc.dt * E_new);
    B = (.5 * pc.dt * E_n) * T_n;
    C = -.5*rhs_ac_arezoo(c_n,T_n,u_n,eta_n,rho_n, pc,0) * pc.dt * pc.rho_water * pc.L;% .* (T_n < 273.15); % IS THIS LINE GOOD?
    D = rho_cp_n .* T_n;
F = -.5*rhs_ac_arezoo_for_implicit_T(c_new,u_new,T_n,eta_new,rho_new, pc,0)*pc.rho_water * pc.L*pc.dt;
    Visc = zeros(pc.N + 2,1);
    for i = 2: pc.N + 1
        Visc(i) = 4/3 * eta_new(i) .* ((u_new(i) - u_new(i - 1)) / pc.dx)^2;
    end    
  
    % Set boundary conditions
    RHS = B + C + D + F+ Visc;
    A_temp(pc.N + 2,:) = 0;
    A_temp(pc.N + 2, pc.N + 2) = 1;
    A_temp(pc.N + 2, pc.N + 1) = -1;
    RHS(pc.N + 2) = 0;
    
    
    A_temp(1,:) = 0;
    
if strcmp(pc.left_BC,'Dirichlet')
    A_temp(1,1) = 1;
    A_temp(1,2) = 1;
    RHS(1) = 2*pc.wall_T;
elseif strcmp(pc.left_BC,'Neumann')
    A_temp(1,1) = 1;
    A_temp(1,2) = -1;
    RHS(1) = 0;
else
    error("Error: Boundary condition not supported.\n")
end
        
    % solve:
    T_new = A_temp\RHS;



end