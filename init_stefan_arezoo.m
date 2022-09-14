function pc = init_stefan_arezoo()
%% Initializes constants of the problem according to the first test case.
N = 1000;
l = 1;
dx = l / (N);
sigma_c = .0317;
%conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;% from paper
dt = 10^-9;
%fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);

pc.N = N;
pc.l = l;
pc.dx = dx;
pc.sigma_c = sigma_c;
pc.ksi_c = 7.0711e-04;% from arezoo, initial multiple term added by me
%pc.ksi_c = 5*dx;%7.0711e-04;% from arezoo, initial multiple term added by me
%pc.ksi_c = 7.0711e-02;% ARTIFICIALLY CHANGED TO TEST T_INIT FUNCTION.

pc.x_init = pc.l/2
pc.thickness_num_pts = pc.ksi_c / dx;
pc.T_M = 1;
pc.L = .53; % latent heat of freezing
pc.dt = dt;
pc.gamma = .5 * 10^-5; %% MIGHT NEED TO TUNE THIS!
pc.mu = pc.gamma; %% might need to tune this too.
pc.nu = 1;
pc.eta_water = 1;
pc.eta_ice = pc.eta_water;

pc.k_water = 0.05;
pc.k_ice = 1; %change

pc.rho_water = 1;
pc.rho_ice = pc.rho_water; % change!

pc.cp_water = 1;
pc.cp_ice = pc.cp_water;% change!

pc.Mc = 1/(pc.nu*pc.ksi_c^2 *pc.L);

pc.lambda = 1/(pc.nu * pc.Mc);
pc.wall_T = 1.1;
pc.init_T = 1.53;
end