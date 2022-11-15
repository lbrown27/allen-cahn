function RUn_ans = RU_n_func_immersed_boundary(rho, eta, pc, c,u)
% Calculates the RHS of the intermediate velocity, u*, calculation.

% initialize vectors for the three variables
A = ones(pc.N+2,1);
B = ones(pc.N+2,1);
C = ones(pc.N+2,1);
D = ones(pc.N+2,1);
cd = (max(rho)/pc.dt + max(eta)/pc.dx^2);
e_d = 10^-3;
for i = 2:pc.N+1
    % define appropriate indices for fractionally-indexed flux vectors
    i_p = i;
    i_m = i - 1;
    % calculate the three terms.
    A(i) = -1/(4*pc.dx) * (rho(i+1)*(u(i_p+1) + u(i_p))^2 - rho(i)*(u(i_p) + u(i_m))^2);
    B(i) = 4/(3 * pc.dx^2)*(eta(i+1)*(u(i_p+1) - u(i_p)) - eta(i) * (u(i_p) - u(i_m)));
    C(i) = 0;%-pc.sigma_c * pc.ksi_c * (c(i + 1) - c(i))/pc.dx; %TURNED OFF SURFACE TENSION ADDED NEGATIVE
    D(i) = cd * (1-(c(i+1)+c(i))/2).^2/((c(i+1)+c(i)/2).^3 + e_d).*(-u(i_p));
end
RUn_ans = A + B + C+D;
i_p = pc.N + 2;
i_m = pc.N + 1;
RUn_ans(pc.N+2) = -1/(4*pc.dx) * (rho(pc.N+2)*(2*u(i_p))^2-rho(pc.N + 2)*(u(i_p) + u(i_m))^2) + ...
    4/(3 * pc.dx^2)*(- eta(pc.N + 2) * (u(i_p) - u(i_m)))+ cd * (1-(c(i_p)+c(i_m))/2).^2/((c(i_p)+c(i_m)/2).^3 + e_d)*(-u(i_p));
; 

% RUn_ans(pc.N + 2) = RUn_ans(pc.N + 1);
i_p = 1;

RUn_ans(1) = -1/(4*pc.dx) * (rho(2)*(u(i_p+1) + u(i_p))^2 - rho(1)*(u(i_p))^2) + ...
    4/(3 * pc.dx^2)*(eta(2)*(u(i_p+1) - u(i_p)) - eta(1) * (u(i_p))) + ...
    cd * (1-(c(i_p)+0)/2).^2/((c(i_p)+0/2).^3 + e_d)*(-u(1));
;%-pc.sigma_c * pc.ksi_c * (c(2) - c(1))/pc.dx;% TURNED OFF SURFACE TENSION HERE TOO ADDED NEGATIVE

%RUn_ans(1) = 0;
end

