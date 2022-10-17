function f = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);

%% Calculate the free energy derivative
g_prime = (4.*c.^3 - 6.*c.^2 + 2.*c);
h_prime = h_prime_func(c);

%A = 3 .* pc.sigma_c ./(pc.ksi_c) .* g_prime; % HILL TERM W/ g'(c)
A = -pc.Mc *pc.lambda /pc.ksi_c^2 * g_prime;
B =  -pc.Mc * pc.rho_water .* pc.L .* ((pc.T_M - T)./pc.T_M).*h_prime; % LATENT HEAT TERM W/ h'(c)
%B = -pc.mu / pc.ksi_c .*(pc.T_M - T).*(c.*(1-c)); % LATENT HEAT TERM W/ h'(c) according to BOETTIGER GEOMETRIC.
fc_deriv = A + B;
%fc_deriv = dfc_dc(c, pc, T,rho);
% for i = 2:pc.N + 1
%     i_plus = i;
%     i_minus = i - 1;
%     if (adv == 1)
%     Advection = -(u(i_plus)*(c(i+1) + c(i)) - u(i_minus) * (c(i) + c(i - 1)))/(2*pc.dx); % ADVECTION TERM
%     else
%         Advection = 0;
%     end

if (adv ==  1)
    Advection =  - div(c,u,pc);
else
    Advection = zeros(pc.N + 2,1);
end

for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;
    laplacian_c = (c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2;
    B = pc.Mc * pc.lambda* laplacian_c; % DIFFUSION-LIKE TERM, added 2 to convert from sqrt(2) stable soln to 2.

    %C = -1*pc.Mc*fc_deriv(i); % FREE ENERGY DERIVATIVE
    div_u = div_stag(u,pc);
    D = div_u(i) * c(i);
    %B = 0;
    %B = pc.Mc * 6 * pc.si_c^2 *(c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2; % DIFFUSION-LIKE TERM, CHANGED!

    f(i) = Advection(i)+B+ fc_deriv(i)+ D;
end
end
