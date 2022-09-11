function pc = init_stefan_arezoo()
%% Initializes constants of the problem according to the first test case.
N = 100;
l = .001;
dx = l / (N);
Mc = 1*10^-3;
sigma_c = .0317;
ksi_c = 1 * dx;
T_M = 273.15;
L = 334000; %latent heat of freezing
dt = 10^-9;
pc.thickness_num_pts = ksi_c / dx;
pc.extra_factor = 1;

pc.N = N;
pc.l = l;
pc.dx = dx;
pc.Mc = Mc;
pc.sigma_c = sigma_c;
pc.ksi_c = ksi_c;
pc.T_M = T_M;
pc.L = L;
pc.dt = dt;

pc.eta_water = 1e-3;
pc.eta_ice = 100;

pc.k_water = 0.005918;
pc.k_ice = .0225; %change

pc.rho_water = 998;
pc.rho_ice = 998; % change!

pc.cp_water = 4200;
pc.cp_ice = 2018;% change!
end