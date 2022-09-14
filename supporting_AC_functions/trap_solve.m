function res = trap_solve(rho_new, eta_new, pc,u_n, rho_n, eta_n, c_n, c_new)
%% See notes in graph paper notebook.

%jacobian = jacobian_calculation(u_n,rho_new,eta_new,pc,c_new);
simple_jacobian = simple_jacobian_calculation(u_n,rho_new,eta_new,pc,c_new);
jacobian = simple_jacobian;
LHS = diag(rho_new) - pc.dt * .5 * jacobian;
A = rho_n .* u_n;
B = pc.dt * .5 * RU_n_func(rho_new, eta_new, pc, c_new,u_n);
C = pc.dt * .5 * RU_n_func(rho_n, eta_n, pc, c_n, u_n);
D = - jacobian * u_n * pc.dt * .5;

RHS = A + B + C + D;
diff = simple_jacobian - jacobian;
max_diff = max(max(simple_jacobian - jacobian));

%fprintf("max jacobian difference = %d \n",max_diff);

% setting bc's:
% LHS(1,:) = 0;
% LHS(1,1) = 1;
% RHS(1) = 0;
% 
% LHS(pc.N+ 2, :) = 0;
% LHS(pc.N+ 2, pc.N +2) = 1;
% LHS(pc.N + 2, pc.N + 1) = -1;
% RHS(pc.N + 2) = 0;

res = LHS \ RHS;


end