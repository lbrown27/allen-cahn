function [f,Advection,Diffusion,Sharpening,Phase_Change] = rhs_ac_ACDI(c,T,u,eta,rho,x_coll, pc,adv)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);
distance_fn = distance_field(c,x_coll,pc);

if (adv == 1)
    
    Advection = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;
    Advection(i) = - (u(i_plus)*(c(i+1) + c(i)) - u(i_minus) * (c(i) + c(i - 1)))/(2*pc.dx); % ADVECTION TERM
end
else
    Advection = zeros(pc.N + 2,1);
end
Diffusion = div_stag(pc.gas_pedal * pc.ksi_c * grad(c,pc),pc); %diffusion term
%C = div_stag(pc.gas_pedal * -1 * stagger_me(c,pc) .* (1 - stagger_me(c,pc)),pc); % CDI term
Sharpening = div_stag(pc.gas_pedal * -1 * stagger_me(c,pc).^2 .* (1 - stagger_me(c,pc)).^2,pc); % EXPERIMENTAL CDI TERM!!!!!!!!!
%C = div_stag(-pc.gas_pedal*.25*(1-tanh(distance_fn/(2*pc.ksi_c)).^2),pc);% ACDI Term! not CDI
%func = c.^2 .* (1 - c).^2;
func = 30.*(c- 1).^2.*c.^2;
Phase_Change = -pc.gas_pedal / (pc.gamma*3*sqrt(2)) .* func.*(pc.T_M - T);
f = Advection + Diffusion + Sharpening + Phase_Change;
end