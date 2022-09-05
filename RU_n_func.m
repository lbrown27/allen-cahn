function RUn_ans = RU_n_func(rho, eta, pc, c,u)
% Calculates the RHS of the intermediate velocity, u*, calculation.

% initialize vectors for the three variables
A = ones(pc.N+2,1);
B = ones(pc.N+2,1);
C = ones(pc.N+2,1);

for i = 2:pc.N+1
    % define appropriate indices for fractionally-indexed flux vectors
    i_p = i;
    i_m = i - 1;
    % calculate the three terms.
    A(i) = -1/(4*pc.dx) * (rho(i+1)*(u(i_p+1) + u(i_p))^2 - rho(i)*(u(i_p) + u(i_m))^2);
    B(i) = 4/(3 * pc.dx^2)*(eta(i+1)*(u(i_p+1) - u(i_p)) - eta(i) * (u(i_p) - u(i_m)));
    C(i) = 0;%-pc.sigma_c * pc.ksi_c * (c(i + 1) - c(i))/pc.dx; %TURNED OFF SURFACE TENSION ADDED NEGATIVE
end
RUn_ans = A + B + C;
i_p = pc.N + 2;
i_m = pc.N + 1;
RUn_ans(pc.N+2) = -1/(4*pc.dx) * (rho(pc.N+2)*(2*u(i_p))^2-rho(pc.N + 2)*(u(i_p) + u(i_m))^2) + ...
    4/(3 * pc.dx^2)*(- eta(pc.N + 2) * (u(i_p) - u(i_m))); 

% RUn_ans(pc.N + 2) = RUn_ans(pc.N + 1);
i_p = 1;

RUn_ans(1) = -1/(4*pc.dx) * (rho(2)*(u(i_p+1) + u(i_p))^2 - rho(1)*(u(i_p))^2) + ...
    4/(3 * pc.dx^2)*(eta(2)*(u(i_p+1) - u(i_p)) - eta(1) * (u(i_p))) + ...
    0;%-pc.sigma_c * pc.ksi_c * (c(2) - c(1))/pc.dx;% TURNED OFF SURFACE TENSION HERE TOO ADDED NEGATIVE

%RUn_ans(1) = 0;
end

