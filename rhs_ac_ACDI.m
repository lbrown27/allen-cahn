function f = rhs_ac_ACDI(c,T,u,eta,rho, pc,adv)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);
if (adv == 1)
    
    A = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;
    A(i) = - (u(i_plus)*(c(i+1) + c(i)) - u(i_minus) * (c(i) + c(i - 1)))/(2*pc.dx); % ADVECTION TERM
end
else
    A = zeros(pc.N + 2,1);
end
B = div_stag(pc.gas_pedal * pc.ksi_c * grad(c,pc),pc);
C = div_stag(pc.gas_pedal * -1 * stagger_me(c,pc) .* (1 - stagger_me(c,pc)),pc);
%func = c.^2 .* (1 - c).^2;
func = 30.*(c- 1).^2.*c.^2;
D = -pc.gas_pedal / (pc.gamma*3*sqrt(2)) .* func.*(pc.T_M - T);
f = A + B + C + D;
end