clear all;
%close all;
clc;

% number of gridpoints. n = 1 represents the ice wall, and n = end
% represents outflow boundary conditions.
% length of domain
% calculate dx based off l and N


%% Define physical constants
pc = init_diml();
x_coll = -pc.dx / 2:pc.dx:pc.l + pc.dx / 2;
x_stag = 0:pc.dx:pc.l + pc.dx;

x_coll = transpose(x_coll);
x_stag = transpose(x_stag);


% initialize the velocity and the phase fields:
u_n = zeros(pc.N+2,1);
u_n(1) = 0;
u_n(pc.N + 2 ) = u_n(pc.N + 1);

u_old = zeros(pc.N+2,1); % store the old velocity


c_n = ones(pc.N + 2,1) * 0;
c_n(1) = -1;
%c_n(2) = -.5;
c_n(pc.N + 2) = c_n(pc.N + 1);
%c_n_experimental = c_init(c_n,x_coll,pc);
%c_n = c_n_experimental; % ADDED THIS LINE, ITS EXPERIMENTAL!
c_n(1) = -2 - c_n(2);
c_old = c_n;

% set the BC that the wall is ice.


%figure(2);
%plot(x_coll, c_n);
% Initialize the temperature field and set the wall temperature.
% First cell is ice.
T_init = 273.15 + .02;
T_n = T_init*ones(pc.N+ 2,1);
wall_temp = -5 +273.15;
T_n(1) = 2 * wall_temp - T_n(2); % wall temperature

% calculate the initial density in the domain based off phase field
% variable
rho_n = (c_n+1)*pc.rho_water - c_n * pc.rho_ice;
rho_old = rho_n; % due to initial condition, this is the same


eta_n = (c_n+1)*pc.eta_water - c_n * pc.eta_ice;
eta_old = eta_n;
eta_new = eta_n;
k_n = (c_n+1)*pc.k_water - c_n * pc.k_ice;

% initially, the pressure field is 0 because there are no pressure
% gradients in the initial field and we know the outflow BC is P = 0.
P = zeros(pc.N + 2,1);
P(1) = P(2);
P(pc.N + 2) = 0 - P(pc.N + 1);

num_iterations = 2500000;
print_interval = 20000;
physical_time = 0;
refinement_factor = 1;
% error_max keeps track of largest discrepancy between the explicit and
% implicit solving method for u*.
%error_max = 0;

% generate the matrix once that matrix_solve needs for pressure
A_pres = gallery('tridiag', pc.N+2, 1,-2,1);
A_pres(1,1) = -1;
A_pres(end,end) = 1;
A_pres(end,end-1) = 0;

count = 0;
trigger = 0;
%for count = 1:num_iterations * refinement_factor
while physical_time < 20
    %% Step 1: Calculate c at time n + 1 using the Allen-Cahn equation.
    % copy old phase field to storage
    
    if count > 1
    c_new = c_step(c_n, T_n,u_n,c_old,eta_new, eta_old,rho_n, pc);
    else 
      c_new = c_step_fe(c_n, T_n,u_n,c_old,eta_new, eta_old,rho_n, pc);
        
    end
    c_diff = c_new - c_n;
    [t_a, t_d, t_f] = calculate_timescales(pc, u_n,rho_n);
    %% Step 2: Obtain density, thermal conductivity, and viscosity at time n+1:
    
    rho_new = (c_new+1)*pc.rho_water - c_new * pc.rho_ice; %this is a new value since calculated with c new.
    %rho_new = ones(pc.N + 2, 1) * pc.rho_water; % ADDED THIS LINE IN DEBUGGING, INVALID FOR TYPICAL PROBLEM
    rho_new_flux = rho_flux(rho_new,pc);
    eta_new = (c_new+1)*pc.eta_water - c_new * pc.eta_ice;
    k_new = (c_new+1)*pc.k_water - c_new * pc.k_ice;
    
    %% Step 3: Find the approximate velocity u* at time step n+1 using the momentum equation SANS pressure term.
    
    %u_star = explicit_solve(rho_n, rho_new, eta_n, pc, c_n, u_n,rho_old,eta_old, c_old,u_old);
    imp_u_star = backward_euler_momentum_solve(rho_new, eta_new, pc,u_n,rho_n, eta_n, c_n, c_new); % NEEDS to be tested to see if BC's work!
    u_star = imp_u_star; % this line switches the code from explicit to implicit!
    %diff_exp_imp = u_star - imp_u_star;
    %% Step 4: Calculate the pressure field by solving a poisson equation
    RHS_Pres = RHS_PE(rho_new,pc,u_n,c_new, T_n,eta_new, u_star);
    P_new = matrix_solve(P,RHS_Pres,pc,A_pres);
    %P_new = gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
    %P_new_acc = accelerated_gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
    
    %% Step 5: Correct velocity to satisfy continuity equation, given the pressure field
    pressure_g = pressure_grad(P_new,pc);
    u_new = u_star - pc.dt * pressure_g ./ rho_new_flux; % need to multiply by 2/3 here if using explicit.
    u_new(1) = 0;
    u_new(pc.N + 2) = u_new(pc.N + 1);
    %u_new(pc.N) = u_new(pc.N - 1);
    %% Step 6: Solve the energy equation
    % TODO: Add viscous dissipation!
    rho_cp_n = pc.rho_water * pc.cp_water * (1 + c_n) + pc.rho_ice * pc.cp_ice * (- c_n);
    
    rho_cp_new = pc.rho_water * pc.cp_water * (1 + c_new) + pc.rho_ice * pc.cp_ice * (- c_new);
    
    %% Step 6.5: correct the temperature and phase field because if T < melting, temp shouldn't decrease until all ice is frozen.
    
    E_new = E_matrix(rho_cp_new, u_new, k_new, pc);
    E_n = E_matrix(rho_cp_n, u_n, k_n, pc);
    
    % Crank nicolson
    A_temp = diag(rho_cp_new) - .5 * pc.dt * E_new;
    B = (.5 * pc.dt * E_n) * T_n;
    C = (-rhs_ac(c_new,T_n,u_new,eta_new,rho_new, pc) * pc.dt * pc.rho_ice * pc.L);% .* (T_n < 273.15); % IS THIS LINE GOOD?
    D = rho_cp_n .* T_n;
    
    Visc = zeros(pc.N + 2,1);
    for i = 2: pc.N + 1
        Visc(i) = 4/3 * eta_new(i) .* ((u_new(i) - u_new(i - 1)) / pc.dx)^2;
    end
    Visc(pc.N + 2) = 0;
    Visc(1) = 0;
    
    %fprintf("Visc max: %d \n", max(Visc));
    % Backward Euler!
    %      A_temp = diag(rho_cp_new) - pc.dt * E_new;
    %     B = zeros(pc.N +2, 1);%(.5 * pc.dt * E_n) * T_n;
    %     C = -(3 * c_new - 4 * c_n + c_old)/2 * pc.rho_ice * pc.L;
    %     D = rho_cp_n .* T_n;
    
    % A * T_{n+1} = B + C + D
    
    RHS = B + C + D + Visc;
    A_temp(pc.N + 2,:) = 0;
    A_temp(pc.N + 2, pc.N + 2) = 1;
    A_temp(pc.N + 2, pc.N + 1) = -1;
    RHS(pc.N + 2) = 0;
    
    
    A_temp(1,:) = 0;
    A_temp(1,1) = 1;
    A_temp(1,2) = 1;
    RHS(1) = 2*wall_temp;
    % TODO: make this faster
    T_new = A_temp\RHS;
    
    %% Step 7: correct the ice velocity to zero:
    %u_new_before_correction = u_new;
    u_new = u_new .* (1 + c_new).* pc.rho_water ./ ((1 + c_new) .* pc.rho_water - c_new .* pc.rho_ice);
    %u_diff = u_new - u_new_before_correction;
    
    %% Step 8: Determine next time step size
    t_v = .5 * pc.dx ^ 2 *pc.rho_ice/ max(eta_new)*10;
    t_c = pc.dx / max(abs(u_new));
    t_s = .5 * sqrt(pc.rho_water / pc.sigma_c * pc.dx^3);
    %t_visc_jump = pc.dx ^ 2 * (rho_ice + rho_water) * .5 /(eta_ice - eta_water);
    TIME = [t_v,t_c,t_s];
    % calculate the next time step based off solver info
    new_timestep =  .5* min(TIME)/refinement_factor;
    new_timestep = 5e-7;
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
    %new_timestep = 10^-7;
    pc.dt = new_timestep;
    %fprintf("%d \n", count);
    count = count + 1;
    %if ((count < print_interval) && (mod(count,print_interval) == 0)) || ((count > print_interval)&& (mod(count,10) == 0))
    if (mod(count,print_interval) == 0)
        fprintf("Time step: %d \n Iteration: %d \n",new_timestep,count);
        fprintf("Physical time: %d \n",physical_time);
        figure(1);
        n_plots = 7;
        subplot(n_plots,1,1);
        plot(x_coll,c_new); % plot phase field after one step
        title("phase field");
        xlim([x_coll(1), x_stag(end)])
        
        subplot(n_plots,1,2);
        plot(x_coll,T_new);
        title('Temp');
        xlim([x_coll(1), x_stag(end)])
        
        subplot(n_plots,1,3);
        plot(x_stag,u_new);
        title("velocity");
        xlim([x_coll(1), x_stag(end)])
        
        subplot(n_plots,1,4);
        plot(x_coll,P_new);
        title('Pressure');
        xlim([x_coll(1), x_stag(end)])
      subplot(n_plots,1,5);
        plot(x_stag,u_star);
        title('U Star');
        xlim([x_coll(1), x_stag(end)])
        
        subplot(n_plots,1,6);
        plot(x_coll, c_diff);
        title('c diff');
        xlim([x_coll(1), x_stag(end)])
        
        
        
        subplot(n_plots,1,7);
        plot(x_stag,RU_n_func(rho_new, eta_new, pc, c_new,u_new));
        title('RUn function');
        xlim([x_coll(1), x_stag(end)])
        drawnow();
        %fprintf("graphs updated. \n");
        %fprintf("c_new(2): %.16f \n",c_new(2));
        recommended_timestep = ac_rec_time(u_new,pc, c_new,T_n);
        %fprintf("TIME SCALES: \n");
%         fprintf("convective: %.16f \n", t_a);
%         fprintf("diffusive: %.16f \n",t_d);
%         fprintf("function: %.16f \n", t_f);
%save(['backup_' num2str(count)]);
    end
    if (trigger == 0)
        if (max(abs(u_new)) > 0)
            trigger = count;
        end
    end
    
    magic = calculate_residual(c_new,T_new,u_new,eta_new,rho_new, pc);
end
subplot(2,2,1);
plot(x_coll,c_new); % plot phase field after one step
title("phase field");

subplot(2,2,2);
plot(x_coll,T_new);
title('Temp');

subplot(2,2,3);
plot(x_stag,u_new);
title("velocity");

subplot(2,2,4);
plot(x_coll,P_new);
title('Pressure');