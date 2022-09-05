function result = rho_flux(rho,pc)
result = zeros(pc.N + 2,1);
for i = 2: pc.N + 1
    result(i) = (rho(i) + rho(i-1))/2;
end
result(1) = (pc.rho_ice + rho(1))/2;
result(pc.N + 2) = rho(pc.N + 1);
end