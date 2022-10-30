function f = rhs_ac_ACDI_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv)
%% does not include some of the phase change term so that it can be pulled over to lhs.
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
D = -pc.gas_pedal * pc.gamma * c.^2 .* (1 - c).^2.*(pc.T_M);
f = A + B + C + D;
end