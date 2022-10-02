function P = matrix_solve(RHS,pc,A)
%% This function solves the poisson equation of form del^2 P = RHS.
% PARAMS: 
% RHS: vector of right-hand-sides of the equation. See solver notes.
% pc: storage struct
% A: The del^2 operator in matrix form.

RHS(1) = 0; % set P(1) = P(2)
RHS(end) = 0; % set final pressure = 0 outflow condition
P = A \ (pc.dx ^ 2 * RHS);
end