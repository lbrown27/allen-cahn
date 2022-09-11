clear all;
close all;
clc;
addpath('supporting_AC_functions')
addpath('stefan_solver')
% number of gridpoints. n = 1 represents the ice wall, and n = end
% represents outflow boundary conditions.
% length of domain
% calculate dx based off l and N


%% Define physical constants
pc = init_diml();
x_coll = transpose(-pc.dx / 2:pc.dx:pc.l + pc.dx / 2);
x_stag = transpose(0:pc.dx:pc.l + pc.dx);

%% User input information:
num_iterations = 2500000;
print_interval = 10000;
new_timestep = 10^-7;
pc.dt = new_timestep;
wall_temp = 273.15-5;
temp_init = 273.15+.2;
%% Initialize field

[c_n,T_n,u_n,rho_n,eta_n,k_n,rho_old] = initialize_fields(pc,wall_temp, temp_init,x_coll);

alpha = find_alpha_fast(pc.k_water, pc.k_ice,pc.L,temp_init,wall_temp, pc.rho_water, pc.cp_water)
mass_orig = sum(rho_n(2:pc.N + 1))*pc.dx;


physical_time = 0;

% error_max keeps track of largest discrepancy between the explicit and
% implicit solving method for u*.
%error_max = 0;

% generate the matrix that matrix_solve needs for pressure
A_pres = gallery('tridiag', pc.N+2, 1,-2,1);
A_pres(1,1) = -1;
A_pres(end,end) = 1;
A_pres(end,end-1) = 0;

count = 0;
imcount = 0;
violates_bounds_count = 0;
t_vec = 0;
loc_vec = find_interface_loc(c_n, x_coll,pc);
loc_ana = loc_vec;
f= figure(1);
vel_on = 0;
c_diff_count = 0;
%for count = 1:num_iterations

%% RESTART WHERE LEFT OFF.
while physical_time < 20
    %while find_interface_loc(c_n, x_coll,pc) < pc.l
    %% Step 1: Calculate c at time n + 1 using the Allen-Cahn equation.
    % copy old phase field to storage
    
    if count > 1
        c_new = c_step(c_n, T_n,u_n,c_old,eta_new, eta_old,rho_n, pc);
    else
        c_new = c_step_fe(c_n, T_n,u_n,eta_n,rho_n, pc);
    end
    c_diff = c_new - c_n;
    for i = 2:pc.N+1
        if c_diff(i) ~= 0
            c_diff_count = c_diff_count + 1;
        end
    end
    
    %[t_a, t_d, t_f] = calculate_timescales(pc, u_n,rho_n);
    %% Step 2: Obtain density, thermal conductivity, and viscosity at time n+1:
    
    rho_new = c_new*pc.rho_water +(1 - c_new) * pc.rho_ice; %this is a new value since calculated with c new.
    rho_new_flux = rho_flux(rho_new,pc);
    eta_new = c_new*pc.eta_water + (1 - c_new) * pc.eta_ice;
    k_new = c_new*pc.k_water +(1 - c_new) * pc.k_ice;
    rho_cp_n = pc.rho_water * pc.cp_water * c_n + pc.rho_ice * pc.cp_ice * (1 - c_n);
    rho_cp_new = pc.rho_water * pc.cp_water * c_new + pc.rho_ice * pc.cp_ice * (1 - c_new);
    %% Step 3: Find the approximate velocity u* at time step n+1 using the momentum equation SANS pressure term.
    if vel_on == 1
        %u_star = explicit_solve(rho_n, rho_new, eta_n, pc, c_n, u_n,rho_old,eta_old, c_old,u_old);
        u_star = backward_euler_momentum_solve(rho_new, eta_new, pc,u_n,rho_n, eta_n, c_n, c_new);
        
        %% Step 4: Calculate the pressure field by solving a poisson equation
        RHS_Pres = RHS_PE(rho_new,pc,u_n,c_new, T_n,eta_new, u_star,rho_n,rho_old);
        P_new = matrix_solve(RHS_Pres,pc,A_pres);
        
        %% Step 5: Correct velocity to satisfy continuity equation, given the pressure field
        pressure_g = pressure_grad(P_new,pc);
        u_new = u_star - pc.dt * pressure_g ./ rho_new_flux; % need to multiply by 2/3 here if using explicit for u_star
        u_new(1) = 0;
        u_new(pc.N + 2) = u_new(pc.N + 1);
        
    else
        u_new = zeros(pc.N + 2,1);
        P_new = zeros(pc.N + 2,1);
    end
    %% Step 6: Solve the energy equation
    T_new = solve_temp_CN(rho_cp_n, rho_cp_new,u_new,u_n,k_new, k_n,T_n, c_new,eta_new,rho_new,wall_temp,pc);
    %% Step 7: correct the ice velocity to zero:
    %u_new_before_correction = u_new;
    %u_new = u_new .* (c_new).* pc.rho_water ./ (c_new .* pc.rho_water + (1 - c_new) .* pc.rho_ice);
    %u_diff = u_new - u_new_before_correction;
    
    %% Step 8: Determine next time step size
    t_v = .5 * pc.dx ^ 2 *pc.rho_ice/ max(eta_new);
    t_c = pc.dx / max(abs(u_new));
    t_s = .5 * sqrt(pc.rho_water / pc.sigma_c * pc.dx^3);
    TIME = [t_v,t_c,t_s];
    % calculate the next time step based off solver info
    %new_timestep =  .5* min(TIME);
    %% Step 9: Set up solver for next loop.
    %percent_lost = mass_counter(rho_new, rho_n, pc, c_new, c_n,mass_orig, x_coll);
    
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
    %pc.dt = new_timestep;
    
    t_vec = [t_vec,physical_time];
    loc_vec = [loc_vec,find_interface_loc(c_n, x_coll,pc)];
    loc_ana = [loc_ana,interface_location(loc_vec(1), alpha,physical_time,pc.l)];
    count = count + 1;
    
    if (mod(count,print_interval) == 0)
        fprintf("Time step: %d \n Iteration: %d \n",new_timestep,count);
        fprintf("Physical time: %d \n",physical_time);
        f= figure(1);
        n_plots = 3;
        plot_count = 1;
        
        subplot(n_plots,1,plot_count);
        plot(x_coll,c_new); % plot phase field after one step
        title("phase field");
        xlim([x_coll(1), x_stag(end)])
        ylim([0 1])
        plot_count = plot_count + 1;
        
        subplot(n_plots,1,plot_count);
        plot(x_coll,T_new);
        title('Temp');
        xlim([x_coll(1), x_stag(end)])
        ylim([wall_temp - 2 temp_init])
        plot_count = plot_count + 1;
        
        
        
        %         subplot(n_plots,1,3);
        %         plot(x_stag,u_new);
        %         title("velocity");
        %         xlim([x_coll(1), x_stag(end)])
        %         %ylim([0 1e-4])
        %
        %         subplot(n_plots,1,4);
        %         plot(x_coll,P_new);
        %         title('Pressure');
        %         xlim([x_coll(1), x_stag(end)])
        %         %ylim([-1 10])
        %
        %
        subplot(n_plots,1,plot_count);
        plot(t_vec,loc_vec);
        hold on;
        plot(t_vec, loc_ana);
        %         hold off;
        %         alpha_num = mean((loc_vec - loc_vec(1))/(loc_ana - loc_vec(1)))*alpha
        legend('num','ana');
        title('Interface Location');
        %
        %
        %
        %         subplot(n_plots,1,6);
        %         hist = 4999;
        %         plot(t_vec(count - hist:end),loc_vec(count - hist:end));
        %         hold on;
        %         [A,B] = rootfit(t_vec(count - hist:end),loc_vec(count - hist:end));
        %         plot(t_vec(count - hist:end), A*sqrt(t_vec(count - hist:end)) + B);
        %         hold off;
        %         legend('num','ana');
        %         title('Fit to last 1000 pts');
        
        %xlim([x_coll(1), x_stag(end)])
        %ylim([0 pc.l])
        %       subplot(n_plots,1,5);
        %         plot(x_stag,u_star);
        %         title('U Star');
        %         xlim([x_coll(1), x_stag(end)])
        %
        %         subplot(n_plots,1,6);
        %         plot(x_coll, c_diff);
        %         title('c diff');
        %         xlim([x_coll(1), x_stag(end)])
        %max(abs(u_diff))
        %magic = calculate_residual(c_new,T_new,u_new,eta_new,rho_new, pc);
        
        %max(magic)
        
        
        
        %         subplot(n_plots,1,7);
        %         plot(x_stag,RU_n_func(rho_new, eta_new, pc, c_new,u_new));
        %         title('RUn function');
        %         xlim([x_coll(1), x_stag(end)])
        drawnow();
        fprintf("graphs updated. \n");
        %imcount = imcount + 1;
        %exportgraphics(f,['Img\img_' num2str(imcount) '.png'],'Resolution',300);
        %save('locs.mat','loc_vec');
        
        
        %fprintf("c_new(2): %.16f \n",c_new(2));
        %recommended_timestep = ac_rec_time(u_new,pc, c_new,T_n);
        %fprintf("TIME SCALES: \n");
        %         fprintf("convective: %.16f \n", t_a);
        %         fprintf("diffusive: %.16f \n",t_d);
        %         fprintf("function: %.16f \n", t_f);
        %save(['backup_' num2str(count)]);
        %fprintf("percent lost: %.16f \n",percent_lost);
    end
    for i = 2:pc.N
        if c_new(i) < -1 || c_new(i) > 0
            violates_bounds_count = violates_bounds_count + 1;
        end
    end
    
end
% figure(1);
% n_plots = 4;
% subplot(n_plots,1,1);
% plot(x_coll,c_new);
% title("phase field");
% xlim([x_coll(1), x_stag(end)])
% ylim([0 1])
%
% subplot(n_plots,1,2);
% plot(x_coll,T_new);
% title('Temp');
% xlim([x_coll(1), x_stag(end)])
% ylim([wall_temp - 2 temp_init])
%
% subplot(n_plots,1,3);
% plot(x_stag,u_new);
% title("velocity");
% xlim([x_coll(1), x_stag(end)])
% ylim([0 1e-3])
%
% subplot(n_plots,1,4);
% plot(x_coll,P_new);
% title('Pressure');
% xlim([x_coll(1), x_stag(end)])
% ylim([-1 10])