function [f,Advection,Diffusion,Sharpening,Phase_Change] = rhs_ac_CDI(c,T,u, pc,adv)
%% returns a vector which is the rhs in the allen-cahn equation
% discretizaition. Valid for all points except the boundary points (1 and N+2)
f = zeros(pc.N + 2,1);

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
stag_c = stagger(c,pc.N + 2);
Diffusion = div_stag(pc.gas_pedal * pc.ksi_c * grad(c,pc),pc.N + 2, pc.dx); %diffusion term

Sharpening= pc.gas_pedal * (-1) * div_stag(stag_c .* (1-stag_c),pc.N + 2,pc.dx); % ORIGINAL!

exponent = 2;
func = (1-c).^exponent.*c.^exponent*.0683*2;
Phase_Change = -pc.gas_pedal / (pc.gamma) .* func.*(pc.T_M - T); % best result was obtained by dividing by 3 * sqrt(5).
f = Advection + Diffusion + Sharpening + Phase_Change;
end