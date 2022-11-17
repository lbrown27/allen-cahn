function res = backward_euler_momentum_solve(rho_new, eta_new, pc,u_n, rho_n, eta_n, c_n, c_new)
%% See notes in graph paper notebook.

%jacobian = jacobian_calculation(u_n,rho_new,eta_new,pc,c_new);
jacobian = simple_jacobian_calculation(u_n,rho_new,eta_new,pc,c_new);
LHS = diag(rho_new) - .5* pc.dt * jacobian;

A = rho_n .* u_n;
 B = pc.dt * .5 * RU_n_func_immersed_boundary(rho_new, eta_new, pc, c_new,u_n);
 C = pc.dt * .5 * RU_n_func_immersed_boundary(rho_n, eta_n, pc, c_n, u_n);

%B = pc.dt * .5 * RU_n_func(rho_new, eta_new, pc, c_new,u_n);
%C = pc.dt * .5 * RU_n_func(rho_n, eta_n, pc, c_n, u_n);
D = -pc.dt*.5 * jacobian *u_n;

RHS = A + B + C+D;

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