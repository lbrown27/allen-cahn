function P = matrix_solve(RHS,pc,A)
RHS(1) = 0; % set P(1) = P(2)
RHS(end) = 0; % set final pressure = 0 outflow condition
P = A \ (pc.dx ^ 2 * RHS);
end