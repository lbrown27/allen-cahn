function pc = init_stefan_dimensional()
%% Initializes constants of the problem according to the first test case.
% note that the values for parameters chosen in this test case are not
% realistic for a real icing problem.
%conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;% from paper
dt = 10^-9;
%fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);

% vel_on will tell the code whether to couple the NS equations with the
% flow.
pc.vel_on = 1;

pc.N = 1000;
pc.l = 1;
pc.dx = pc.l / (pc.N);
pc.sigma_c = .0317;
sigma_c_works = 0.317;
pc.ksi_c = 1.5*7.0711e-04*1000/pc.N;% from arezoo, initial multiple term added by me
ksi_c_works = 1.5*7.0711e-04*1000/pc.N;% from arezoo, initial multiple term added by me

%pc.ksi_c = 5*dx;%7.0711e-04;% from arezoo, initial multiple term added by me
%pc.ksi_c = 7.0711e-02;% ARTIFICIALLY CHANGED TO TEST T_INIT FUNCTION.

pc.x_init = .8; % initial location for the water-ice interface.
pc.thickness_num_pts = pc.ksi_c / pc.dx;
pc.T_M = 1;
T_M_works = 1;
pc.L = 334000; % latent heat of freezing
L_works = 1;
pc.dt = dt;
%pc.gamma = 2.5 * 10^-4; %% MIGHT NEED TO TUNE THIS!
%pc.mu = 4*10^3; %% might need to tune this too.
rho_works = 1;

%% EXPERIMENTAL:

%pc.gamma = 2.5 * 10^-4; %% MIGHT NEED TO TUNE THIS!

%% END EXPERIMENTAL.

pc.eta_water = 1e-3;
pc.eta_ice = 100;
pc.left_BC = 'Neumann';
pc.left_material = 'ice';
pc.k_water = .5918;
pc.k_ice = 2.25; %change

pc.rho_water = 998;
expansion_factor = 1;
pc.rho_ice = 998;%pc.rho_water*expansion_factor; % change!

% pc.rho_water = 1;
% expansion_factor = 1;
% pc.rho_ice = 1;%pc.rho_water*expansion_factor; % change!

pc.cp_water = 4200;
pc.cp_ice = 2018;% change!

pc.thermal_diff_ice = pc.k_ice / (pc.rho_ice * pc.cp_ice);
pc.thermal_diff_water = pc.k_water / (pc.rho_water * pc.cp_water);


pc.mu = 4*10^3/ksi_c_works * pc.ksi_c *rho_works / pc.rho_water; %% might need to tune this too.
pc.mu = 4*10^4/ksi_c_works * pc.ksi_c *rho_works / pc.rho_water; %% might need to tune this too.

%pc.gamma = 2.5 * 10^-4; %% MIGHT NEED TO TUNE THIS!

%pc.Mc = 1/(pc.nu*pc.ksi_c^2 *pc.L);% Lazy definition from paper.

%pc.lambda = 1/(pc.nu * pc.Mc); % Lazy definition from paper.

pc.gamma = 2.5 * 10^-4*pc.sigma_c/sigma_c_works*(L_works/T_M_works)*(pc.T_M/pc.L)/(pc.rho_water); %% MIGHT NEED TO TUNE THIS!

pc.lambda = 3 * sqrt(2) * pc.rho_water * pc.L * pc.gamma * pc.ksi_c/pc.T_M;


pc.Mc = pc.gamma * pc.mu/pc.lambda;


pc.wall_T = -5+ pc.T_M;
pc.init_T = 0.2 + pc.T_M;

%pc.alpha = find_alpha_fast(pc.k_water, pc.k_ice,pc.L,pc.init_T,pc.wall_T, pc.rho_water, pc.cp_water,pc.T_M);
pc.alpha = -1e-10;
end