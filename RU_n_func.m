function ans = RU_n_func(rho, eta, pc, c,u)
% Calculates the RHS of the intermediate velocity, u*, calculation.

% initialize vectors for the three variables
A = ones(1,pc.N);
B = ones(1,pc.N);
C = ones(1,pc.N);

% go through for loop for i = 3 to N - 1 because the loop will access i + 1
% and i - 2 for each element.
for i = 3:pc.N-1
    % define appropriate indices for fractionally-indexed flux vectors
    i_plus = i;
    i_minus = i - 1;
    % calculate the three terms.
    A(i) = -1 * ((rho(i+1) + rho(i))*u(i_plus)^2 - (rho(i) + rho(i - 1)) * u(i_minus)^2)/(2 * pc.dx);
    B(i) = 2/(3*pc.dx^2) ...
        * (eta(i+1) * (u(i_plus + 1) - u(i_plus)) - eta(i - 1) * (u(i_minus) - u(i_minus - 1)));
    C(i) = pc.ksi_c * (c(i+1) - c(i-1))/(2 * pc.dx);
   
end

ans = A + B + C;
ans(pc.N) = ans(pc.N - 1); % outflow bc, the velocity at the boundary will be the same as N - 1. 
ans(1) = 0; % wall bc, the wall will be stationary.
for i = 2 % go through the loop for i = 2 since this was skipped before.
       i_plus = i;
    i_minus = i - 1;
    A(i) = -1 * ((rho(i+1) + rho(i))*u(i_plus)^2 - (rho(i) + rho(i - 1)) * u(i_minus)^2)/(2 * pc.dx);
    B(i) = 2/(3*pc.dx^2) ...
        * (eta(i+1) * (u(i_plus + 1) - u(i_plus)) - eta(i - 1) * (u(i_minus) - 0));
    C(i) = pc.ksi_c * (c(i+1) - c(i-1))/(2 * pc.dx);
    ans(2) = A(2) + B(2) + C(2);
end
end

