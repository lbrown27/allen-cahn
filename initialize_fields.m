function [c_n,T_n,u_n,rho_n,eta_n,k_n,rho_old] = initialize_fields(pc, x_coll)

% initialize the velocity and the phase fields:
u_n = zeros(pc.N+2,1);
u_n(1) = 0;
u_n(pc.N + 2 ) = u_n(pc.N + 1);

c_n = ones(pc.N + 2,1);

% Add these two lines if you want the profile to start with a diffuse htan
% profile.

c_n_experimental = c_init(c_n,x_coll,pc);
c_n = c_n_experimental;

% Set boundary conditions
c_n(pc.N + 2) = c_n(pc.N + 1);
c_n(1) = - c_n(2);

% Initialize the temperature field and set the wall temperature.
% Also set the wall temp at i = 1/2
T_n =ones(pc.N + 2,1).*( pc.wall_T.*(x_coll < pc.l/2) + (x_coll>=pc.l/2).* pc.init_T);
T_n = T_init(T_n,x_coll,pc);
%T_n = ones(pc.N + 2,1)*pc.init_T;
T_n(1) = 2 * pc.wall_T - T_n(2); % wall temperature

% calculate the initial density in the domain based off phase field
% variable
rho_n = c_n*pc.rho_water + (1-c_n) * pc.rho_ice;
rho_old = rho_n; % due to initial condition, this is the same


eta_n = c_n*pc.eta_water +(1-c_n) * pc.eta_ice;
eta_old = eta_n;
eta_new = eta_n;
k_n = c_n*pc.k_water + (1- c_n) * pc.k_ice;

end