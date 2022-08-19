function P = matrix_solve(P_old,RHS,pc)
% performs gauss seidel to solve the pressure poisson equation
A = gallery('tridiag', pc.N, 1,-2,1);
A(1,1) = -1;
A(end,end) = 1;
A(end,end-1) = 0;
RHS = transpose(RHS);
RHS(end) = 0; % set final pressure = 0 outflow condition
P = A \ (pc.dx ^ 2 * RHS);
%P(j) = - .5 * (pc.dx ^2 * RHS(j) - P_old(j + 1) - P(j - 1));

end