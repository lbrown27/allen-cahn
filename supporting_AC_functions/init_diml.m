function pc = init_diml()
%% Initializes constants of the problem according to the first test case.
N = 100;
l = .001;
dx = l / (N);
Mc = 1*10^-5;
pc.vel_on = 1;

%conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;% from paper
ksi_c = 3 * dx;
T_M = 273.15;
L = 334000; %latent heat of freezing
dt = 10^-9;
pc.thickness_num_pts = ksi_c / dx;
fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);

pc.left_BC = 'Neumann';
pc.left_material = 'ice';

pc.N = N;
pc.l = l;
pc.x_init = .8*pc.l; % initial location for the water-ice interface.

pc.dx = dx;
pc.Mc = Mc;
pc.sigma_c = .0317;
pc.ksi_c = ksi_c;
pc.T_M = T_M;
pc.L = L;
pc.dt = dt;
pc.gamma = 1.3 * 10^-3; %% MIGHT NEED TO TUNE THIS!
pc.eta_water = 1e-3;
pc.eta_ice = 100;

pc.k_water = 0.005918;
pc.k_ice = .0225; %change

pc.rho_water = 998;
pc.rho_ice = 998; % change!

pc.cp_water = 4200;
pc.cp_ice = 2018;% change!

pc.wall_T = -5 + pc.T_M;
pc.init_T = 0.1 + pc.T_M;

pc.lambda = 3 * sqrt(2) * pc.rho_water*pc.L*pc.gamma *pc.ksi_c/pc.T_M;

pc.alpha = find_alpha_fast(pc.k_water, pc.k_ice,pc.L,pc.init_T,pc.wall_T, pc.rho_water, pc.cp_water,pc.T_M);
end