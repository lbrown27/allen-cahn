function f = rhs_ac(c,T,u,pc)
% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. 
f = zeros(pc.N,1);
for i = 2:pc.N-1
    i_plus = i;
    i_minus = i - 1;
    A = -(u(i_plus)*(c(i+1) + c(i)) - u(i_minus) * (c(i) + c(i - 1)))/(2*pc.dx);
    B = pc.Mc * 6 * pc.sigma_c * pc.ksi_c *(c(i+1) - 2 * c(i) + c(i-1))/pc.dx^2;
    C = -pc.Mc*(3 * pc.sigma_c /pc.ksi_c * (4*c(i)^3 + 6*c(i)^2 + 2*c(i)) ...
        + pc.L * ((pc.T_M - T(i))/pc.T_M)*(30 *(c(i)+1)^2*c(i)^2));
    f(i) = A+B+C;
end
f(1) = -1;
f(end) = 0;
end