function T_new = solve_temp_CN(rho_cp_n, rho_cp_new,u_new,u_n,k_new, k_n,T_n,c_new, eta_new,rho_new,wall_temp,pc)

    E_new = E_matrix_shifted(rho_cp_new, u_new, k_new, pc);
    E_n = E_matrix_shifted(rho_cp_n, u_n, k_n, pc);
    
    % Crank nicolson
    A_temp = sparse(diag(rho_cp_new) - .5 * pc.dt * E_new);
    B = (.5 * pc.dt * E_n) * T_n;
    C = -rhs_ac_arezoo(c_new,T_n,u_new,eta_new,rho_new, pc) * pc.dt * pc.rho_water * pc.L;% .* (T_n < 273.15); % IS THIS LINE GOOD?
    D = rho_cp_n .* T_n;
    F = zeros(pc.N + 2,1);
    for i = 2: pc.N + 1
        i_plus = i;
        i_minus = i-1;
        F(i) =  (-(u_new(i_plus)*(c_new(i+1) + c_new(i)) - u_new(i_minus) * (c_new(i) + c_new(i - 1)))/(2*pc.dx)) * pc.rho_water * pc.L* pc.dt;
    end
    Visc = zeros(pc.N + 2,1);
    for i = 2: pc.N + 1
        Visc(i) = 4/3 * eta_new(i) .* ((u_new(i) - u_new(i - 1)) / pc.dx)^2;
    end    
  
    % Set boundary conditions
    RHS = B + C + D + Visc+ F;
    A_temp(pc.N + 2,:) = 0;
    A_temp(pc.N + 2, pc.N + 2) = 1;
    A_temp(pc.N + 2, pc.N + 1) = -1;
    RHS(pc.N + 2) = 0;
    
    
    A_temp(1,:) = 0;
    A_temp(1,1) = 1;
    A_temp(1,2) = 1;
    RHS(1) = 2*wall_temp;
        
    % solve:
    T_new = A_temp\RHS;



end