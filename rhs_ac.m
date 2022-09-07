function f = rhs_ac(c,T,u,eta,rho, pc)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);
fc_deriv = dfc_dc(c, pc, T,rho);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;
    A = -(u(i_plus)*(c(i+1) + c(i)) - u(i_minus) * (c(i) + c(i - 1)))/(2*pc.dx); % ADVECTION TERM
    B = pc.extra_factor * pc.Mc * 6 * pc.sigma_c * pc.ksi_c *(c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2; % DIFFUSION-LIKE TERM
    
    C = -1*pc.Mc*fc_deriv(i); % FREE ENERGY DERIVATIVE
    D = (u(i_plus) - u(i_minus))/pc.dx * c(i);
    %B = 0;
    %B = pc.Mc * 6 * pc.ksi_c^2 *(c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2; % DIFFUSION-LIKE TERM, CHANGED!

    f(i) = A+B+ C + D;
    
end
end