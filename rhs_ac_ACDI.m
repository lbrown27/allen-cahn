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
stag_c = stagger(c,pc.N + 2);
Diffusion = div_stag(pc.gas_pedal * pc.ksi_c * grad(c,pc),pc.N + 2, pc.dx); %diffusion term
%Sharpening = div_stag(pc.gas_pedal * -1 * 4 * stag_c.^2 .* (1 - stag_c).^2,pc); % CDI term. ADDED A FACTOR OF FOUR.



deriv_c = div_stag(stag_c,pc.N + 2,pc.dx);
big_factor = 1/pc.ksi_c;

theta = max(abs(deriv_c + big_factor/pc.l))*pc.ksi_c;
stag_distance = stagger(distance_fn,pc.N + 2);
Sharpening = pc.gas_pedal * (-1) * (1/theta) * (4 * stag_c.^3 -6*stag_c.^2+2*stag_c).*div_stag(stag_c + big_factor *stag_distance/pc.l,pc.N + 2,pc.dx); % TERM ADDED IN THE DIV_STAG!@!!!!!!!! 
theta_orig = max(abs(deriv_c))*pc.ksi_c;
Sharpening_ORIGINAL = pc.gas_pedal * (-1) * (1/theta_orig) * (4 * stag_c.^3-6*stag_c.^2+2*stag_c).*div_stag(stag_c,pc.N + 2,pc.dx); % ORIGINAL!
Sharpening_ORIGINAL_UNSCALED= pc.gas_pedal * (-1) * (4 * stag_c.^3-6*stag_c.^2+2*stag_c).*div_stag(stag_c,pc.N + 2,pc.dx); % ORIGINAL!
Sharpening_ORIGINAL_UNSCALED= pc.gas_pedal * (-1) * div_stag(stag_c .* (1-stag_c),pc.N + 2,pc.dx); % ORIGINAL!
%% DEBUGGING ZONE:
% original copied from allen-cahn method.
% 
% g_prime = (4.*c.^3 - 6.*c.^2 + 2.*c); % please delete when dne experimenting with line below.
% Sharpening_ac = -pc.Mc *pc.lambda /pc.ksi_c^2 * g_prime; % This term is borrowed from ac. Please delete!!!!!
% figure(56);
% plot(x_coll, Sharpening);
% hold on;
% plot(x_coll,Sharpening_ac);
% plot(x_coll, Sharpening_ORIGINAL);
% plot(x_coll,Sharpening_ORIGINAL_UNSCALED)
% legend("CDI method","Allen-Cahn Method", "Original CDI","Original CDI Unscaled");

% deriv_c = div_stag(stag_c,pc.N + 2,pc.dx);
% theta = max(abs(deriv_c))*pc.ksi_c;
% 
% figure(13);
% plot(x_coll,deriv_c/theta);
% 
% hold on;
% 
% plot_thickness_const = ones(pc.N + 2,1)/pc.ksi_c;
% plot(x_coll,plot_thickness_const);



%% DEBUGGING ZONE OVER.



%Sharpening_equivalent = pc.gas_pedal *(-4)*(4*stag_c.^3 - 6*stag_c.^2 + 2.*stag_c)*
%experimental_function(c,pc);


%% Resume non-experimental:

%Sharpening = div_stag(pc.gas_pedal * -1 * stagger_me(c,pc).^2 .* (1 - stagger_me(c,pc)).^2,pc); % EXPERIMENTAL CDI TERM. Has squared function. 
%C = div_stag(-pc.gas_pedal*.25*(1-tanh(distance_fn/(2*pc.ksi_c)).^2),pc);% ACDI Term! not CDI
%func = c.^2 .* (1 - c).^2;

Sharpening = Sharpening_ORIGINAL_UNSCALED;
exponent = 5.3219/2;
exponent = 2;
func = (1-c).^exponent.*c.^exponent;
Phase_Change = -pc.gas_pedal / (pc.gamma*3*sqrt(2)) .* func.*(pc.T_M - T);
f = Advection + Diffusion + Sharpening + Phase_Change;
end