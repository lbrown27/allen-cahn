function pc = init_diml_arezoo()
%% Initializes constants of the problem according to the first test case.
N =100;
l = .001;
dx = l / (N);
sigma_c = .0317;
%conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;% from paper
ksi_c = 1 * dx;
T_M = 273.15;
L = 334000; %latent heat of freezing
dt = 10^-9;
pc.thickness_num_pts = ksi_c / dx;
fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);

pc.N = N;
pc.l = l;
pc.dx = dx;
pc.sigma_c = sigma_c;
pc.ksi_c = ksi_c;
pc.T_M = T_M;
pc.L = L;
pc.dt = dt;
pc.gamma = .5 * 10^-4; %% MIGHT NEED TO TUNE THIS!
pc.mu = pc.gamma; %% might need to tune this too.

pc.eta_water = 1e-3;
pc.eta_ice = 100;

pc.k_water = 0.005918;
pc.k_ice = .0225; %change

pc.rho_water = 998;
pc.rho_ice = 998; % change!

pc.cp_water = 4200;
pc.cp_ice = 2018;% change!

pc.lambda = 3 * sqrt(2) * pc.rho_water*pc.L*pc.gamma *pc.ksi_c/pc.T_M;

pc.Mc = pc.mu * pc.gamma / pc.lambda;

end