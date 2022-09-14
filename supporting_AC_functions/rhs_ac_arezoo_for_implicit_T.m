function f = rhs_ac_arezoo_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv);
%% does not include some of the phase change term so that it can be pulled over to lhs.
f = zeros(pc.N + 2,1);

%% Calculate the free energy derivative
g_prime = (4.*c.^3 - 6.*c.^2 + 2.*c);
h_prime = h_prime_func(c,T,pc.T_M);

%A = 3 .* pc.sigma_c ./(pc.ksi_c) .* g_prime; % HILL TERM W/ g'(c)
A = -pc.Mc *pc.lambda /pc.ksi_c^2 * g_prime;
B =  -pc.Mc * pc.rho_water .* pc.L *h_prime; % LATENT HEAT TERM W/ h'(c)

fc_deriv = A + B;
%fc_deriv = dfc_dc(c, pc, T,rho);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;

    laplacian_c = (c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2;
    B = pc.Mc * pc.lambda* laplacian_c; % DIFFUSION-LIKE TERM
    
    %C = -1*pc.Mc*fc_deriv(i); % FREE ENERGY DERIVATIVE
    D = (u(i_plus) - u(i_minus))/pc.dx * c(i);
    %B = 0;
    %B = pc.Mc * 6 * pc.ksi_c^2 *(c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2; % DIFFUSION-LIKE TERM, CHANGED!

    f(i) =B+ fc_deriv(i)+ D;
    
end
end