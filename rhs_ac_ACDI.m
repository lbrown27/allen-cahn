function f = rhs_ac_ACDI(c,T,u,eta,rho, pc,adv)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);
if (adv == 1)
    A = -div(c,u,pc);
else
    A = zeros(pc.N + 2,1);
end
B = div_stag(pc.gas_pedal * pc.ksi_c * grad(c,pc),pc);
C = div_stag(pc.gas_pedal * -1 * stagger_me(c,pc) .* (1 - stagger_me(c,pc)),pc);
D = -pc.gas_pedal * pc.gamma * c .* (1 - c).*(pc.T_M - T);
f = A + B + C + D;
end
