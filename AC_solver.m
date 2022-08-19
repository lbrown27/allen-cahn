clear all;
close all;
clc;

% number of gridpoints. n = 1 represents the ice wall, and n = end
% represents outflow boundary conditions.
N = 40;
% length of domain
l = .0004;
% calculate dx based off l and N
dx = l / (N-1);

x = 0:dx:l;

%% Define physical constants
Mc = 2 * 10^-5;
sigma_c = .0317;
ksi_c = 2.64 * 10^-4;
T_M = 273.15;

eta_water = 1e-3;
eta_ice = 100;

k_water = 0.5918;
k_ice = 2.25;

rho_water = 998;
rho_ice = 898;

cp_water = 4200;
cp_ice = 2018;

L = 334000; %latent heat of freezing

% Initialize the temperature field and set the wall temperature.
% First cell is ice.
T_n = (273.15 + 10)*ones(N,1);
T_n(1) = 273.15;

% set the first time step
dt = 10^-10;

% initialize the velocity and the phase fields:
u_n = zeros(1,N);
u_old = zeros(1,N); % store the old velocity

c_n = zeros(1,N);
c_old = c_n;
% set the BC that the wall is ice.
c_n(1) = -1;

% calculate the initial density in the domain based off phase field
% variable
rho_n = (c_n+1)*rho_water - c_n * rho_ice;
rho_old = rho_n; % due to initial condition, this is the same


eta_n = (c_n+1)*eta_water - c_n * eta_ice;
eta_old = eta_n;

k_n = (c_n+1)*k_water - c_n * k_ice;

% initially, the pressure field is 0 because there are no pressure
% gradients in the initial field and we know the outflow BC is P = 0.
P = zeros(1,N);

% link physical constants (pc) into one structure to ease function calls
pc.N = N;
pc.l = l;
pc.dx = dx;
pc.Mc = Mc;
pc.sigma_c = sigma_c;
pc.ksi_c = ksi_c;
pc.T_M = T_M;
pc.L = L;
pc.dt = dt;

num_iterations = 10000;
print_interval = 500;
physical_time = 0;
refinement_factor = 1;

for count = 1:num_iterations * refinement_factor
    %% Step 1: Calculate c at time n + 1 using the Allen-Cahn equation.
    % copy old phase field to storage
    c_new = c_step(c_n, T_n,u_n,c_old, pc);
    
    %% Step 2: Obtain density, thermal conductivity, and viscosity at time n+1:
    
    rho_new = (c_new+1)*rho_water - c_new * rho_ice; %this is a new value since calculated with c new.
    eta_new = (c_new+1)*eta_water - c_new * eta_ice;
    k_new = (c_new+1)*k_water - c_new * k_ice;
    
    %% Step 3: Find the approximate velocity u* at time step n+1 using the momentum equation SANS pressure term.
    
    RU_n =  RU_n_func(rho_n, eta_n, pc,c_n,u_n);
    RU_n_old = RU_n_func(rho_old, eta_old, pc, c_old, u_old);
    A = jacobian_calculation(u_old, u_n);
    u_star = (4 .* pc.dt .* RU_n  -  2 .* pc.dt .* RU_n_old + 4 .* rho_n .* u_n - rho_old .* u_old)./(3.* rho_new);
    
    
    %% Step 4: Calculate the pressure field by solving a poisson equation
    P_new = matrix_solve(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
    %P_new = gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
    %P_new_acc = accelerated_gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
    
    %% Step 5: Correct velocity to satisfy continuity equation, given the pressure field
    pressure_g = pressure_grad(P_new,pc);
    u_new = u_star - pc.dt * pressure_g ./ rho_new;
    u_new(pc.N) = u_new(pc.N - 1);
    
    %% Step 6: Solve the energy equation
    % TODO: Add viscous dissipation!
    rho_cp_n = rho_water * cp_water * (1 + c_n) + rho_ice * cp_ice * (- c_n);
    
    rho_cp_new = rho_water * cp_water * (1 + c_new) + rho_ice * cp_ice * (- c_new);
    
    E_new = E_matrix(rho_cp_new, u_new, k_new, pc);
    E_n = E_matrix(rho_cp_n, u_n, k_n, pc);
    
    A = diag(rho_cp_new) - .5 * pc.dt * E_new;
    B = (.5 * pc.dt * E_n) * T_n;
    B = transpose(B);
    C = -(3 * c_new - 4 * c_n + c_old)/2 * rho_ice * pc.L;
    D = rho_cp_n .* transpose(T_n);
    
    % A * T_{n+1} = B + C + D
    
    RHS = B + C + D;
    
    RHS = transpose(RHS);
    
    % TODO: make this faster
    T_new = A\RHS;
    T_new(1) = 273.15;
    T_new(end) = T_n(end);
    
    %% Step 7: correct the ice velocity to zero:
    %u_new_before_correction = u_new;
    u_new = u_new .* (1 + c_new).* rho_water ./ ((1 + c_new) .* rho_water - c_new .* rho_ice);
    %u_diff = u_new - u_new_before_correction;
    
    %% Step 8: Determine next time step size
    t_v = .25 * pc.dx ^ 2 *rho_ice/ max(eta_new);
    t_c = pc.dx / max(abs(u_new));
    t_s = .5 * sqrt(rho_water / sigma_c * pc.dx^3);
    t_visc_jump = pc.dx ^ 2 * (rho_ice + rho_water) * .5 /(eta_ice - eta_water);
    A = [t_v,t_c,t_s,t_visc_jump];
    % calculate the next time step based off solver info
    new_timestep = .5 * min(A)/refinement_factor;
    
    %% Step 9: Set up solver for next loop.
    
    c_old = c_n;
    c_n = c_new;
    eta_old = eta_n;
    eta_n = eta_new;
    k_n = k_new;
    P = P_new;
    rho_old = rho_n;
    rho_n = rho_new;
    T_n = T_new;
    u_old = u_n;
    u_n = u_new;
    
    physical_time = physical_time + pc.dt;
    pc.dt = new_timestep;
    %fprintf("%d \n", count);
    if (mod(count,print_interval) == 0)
        
        fprintf("Time step: %d \n Iteration: %d \n",new_timestep,count);
        fprintf("Physical time: %d \n",physical_time);
        subplot(2,3,1);
        plot(x,c_new); % plot phase field after one step
        title("phase field");
        
        subplot(2,3,2);
        plot(x,T_new);
        title('Temp');
        
        subplot(2,3,3);
        plot(x,u_new);
        title("velocity");
        
        subplot(2,3,4);
        plot(x,P_new);
        title('Pressure');
        
        subplot(2,3,5);
        plot(x,RU_n);
        title('RU_n');
        drawnow();
        fprintf("graphs updated. \n");
        
    end
end
subplot(2,2,1);
plot(x,c_new); % plot phase field after one step
title("phase field");

subplot(2,2,2);
plot(x,T_new);
title('Temp');

subplot(2,2,3);
plot(x,u_new);
title("velocity");

subplot(2,2,4);
plot(x,P_new);
title('Pressure');

figure(2);
plot(x,RU_n);
title('RU_n');
% subplot(3,3,5);
% plot(x,eta_new);
% title('Eta');

% subplot(3,3,6);
% plot(x,k_new);
% title('Thermal Conductivity');

% subplot(2,3,5);
% plot(x,rho_new);
% title('Density');