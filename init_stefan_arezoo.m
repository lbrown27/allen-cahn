function pc = init_stefan_arezoo()
%% Initializes constants of the problem according to the first test case.
% note that the values for parameters chosen in this test case are not
% realistic for a real icing problem.
%conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;% from paper
dt = 10^-10;
%fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);


% phase model. OPTIONS: "allen-cahn", "acdi", and "cdi".
pc.phase_model = "allen-cahn";

% vel_on will tell the code whether to couple the NS equations with the
% flow.
pc.vel_on = 1;

pc.N = 100;
pc.l = 1;
pc.dx = pc.l / (pc.N);
pc.sigma_c = .0317;
%pc.sigma_c = .0001325;

pc.ksi_c = 1.5*7.0711e-04*1000/pc.N;% from arezoo, initial multiple term added by me
%pc.ksi_c = 5*dx;%7.0711e-04;% from arezoo, initial multiple term added by me
%pc.ksi_c = 7.0711e-02;% ARTIFICIALLY CHANGED TO TEST T_INIT FUNCTION.

pc.x_init = .8; % initial location for the water-ice interface.
pc.thickness_num_pts = pc.ksi_c / pc.dx;
pc.T_M = 1;%273.15;
pc.L = .53; % latent heat of freezing
pc.dt = dt;
pc.rho_water = 1;

pc.gamma = 2.5 * 10^-4; %% MIGHT NEED TO TUNE THIS! ORIGINAL!
pc.gamma = pc.sigma_c * pc.T_M / (pc.rho_water*pc.L); % from definition given by Boettinger.
pc.mu = 4*10^3; %% might need to tune this too.


%% EXPERIMENTAL:

pc.gamma = 2.5 * 10^-4; %% MIGHT NEED TO TUNE THIS!
pc.mu = 4*10^3; %% might need to tune this too.
pc.mu = 4*10^3; %% might need to tune this too.

%% END EXPERIMENTAL.

pc.nu = 1; % used in the lazy definition for Mc and lambda, see the Arezoo paper.
pc.eta_water = 1;
pc.eta_ice = pc.eta_water;
pc.left_BC = 'Neumann';
pc.left_material = 'ice';
pc.k_water = 0.05;
pc.k_ice = 1; %change

pc.rho_ice = pc.rho_water*.9999999; % change!

pc.cp_water = 1;
pc.cp_ice = pc.cp_water;% change!

pc.thermal_diff_ice = pc.k_ice / (pc.rho_ice * pc.cp_ice);
pc.thermal_diff_water = pc.k_water / (pc.rho_water * pc.cp_water);

%pc.Mc = 1/(pc.nu*pc.ksi_c^2 *pc.L);% Lazy definition from paper.

%pc.lambda = 1/(pc.nu * pc.Mc); % Lazy definition from paper.

pc.lambda = 3 * sqrt(2) * pc.rho_water * pc.L * pc.gamma * pc.ksi_c/pc.T_M;
pc.Mc = pc.gamma * pc.mu/pc.lambda;


pc.wall_T = .1 + pc.T_M;
pc.init_T = 0.53 + pc.T_M;

% pc.wall_T = pc.T_M;
% pc.init_T = pc.T_M;

pc.alpha = find_alpha_fast(pc.k_water, pc.k_ice,pc.L,pc.init_T,pc.wall_T, pc.rho_water, pc.cp_water,pc.T_M);
pc.gas_pedal = pc.gamma * pc.mu / pc.ksi_c;

print_coefficients(pc);
end