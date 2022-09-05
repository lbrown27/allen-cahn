function pc = init_diml()
%% Initializes constants of the problem according to the first test case.
N = 15;
l = .00015;
dx = l / (N);
Mc = 1*10^-5;
sigma_c = .0317;
conversion_smallness_factor = 2.64*10^-4 / (4 * 10^-3);
ksi_c = l * conversion_smallness_factor;
%ksi_c = 26.4 * dx;
T_M = 273.15;
L = 334000; %latent heat of freezing
dt = 10^-9;
pc.thickness_num_pts = ksi_c / dx;
fprintf("Thickness is %f points wide (aim for at least 20!)\n",pc.thickness_num_pts);
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

pc.k_water = 0.5918;
pc.k_ice = 2.25; %change

pc.rho_water = 998;
pc.rho_ice = 898; % change!

pc.cp_water = 4200;
pc.cp_ice = 2018;% change!
end