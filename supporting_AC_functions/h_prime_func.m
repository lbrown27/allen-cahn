function h_prime = h_prime_func(c)
% h_prime = (c ==1) .* (T<=T_M) + (c == 0).*(T>=T_M);
% h_prime = h_prime + (h_prime == 0).*(30 .*(c-1).^2.*c.^2);
h_prime = (30 .*(c-1).^2.*c.^2);

end