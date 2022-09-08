function res = dfc_dc(c, pc, T,rho)
A = 3 .* pc.sigma_c ./(pc.ksi_c) .* (4.*c.^3 - 6.*c.^2 + 2.*c);
B =  rho .* pc.L .* ((pc.T_M - T)./pc.T_M).*(30 .*(c-1).^2.*c.^2); % LATENT HEAT TERM, added negative

%A = zeros(pc.N + 2,1);
%B = zeros(pc.N + 2, 1);
res = A + B;
end