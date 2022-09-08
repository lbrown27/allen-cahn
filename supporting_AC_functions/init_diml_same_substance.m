function pc = init_diml_same_substance()
%% Initializes constants of the problem according to the first test case.
N = 10;
l = .00001;
dx = l / (N);
Mc = 2 *10^-5;
sigma_c = .0317;
%onversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
%ksi_c = l * conversion_smallness_factor;
ksi_c = 26.4 * dx;
T_M = 273.15;
L = 334000; %latent heat of freezing
dt = 10^-7;

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
pc.eta_ice = 1e-3;

pc.k_water = 0.5918;
pc.k_ice = .5918;

pc.rho_water = 998;
pc.rho_ice = 998;

pc.cp_water = 4200;
pc.cp_ice = 4200;
end