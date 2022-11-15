% is exponent chosen for phase change term valid?
test_station = .05;
C = 1; % for n = 2
n_min_2 = 1 + log(pc.ksi_c * C/(3*sqrt(2)*pc.gamma))/log(1/test_station)


C = 140; % for n = 3
n_min_3 = 1 + log(pc.ksi_c * C/(3*sqrt(2)*pc.gamma))/log(1/test_station)

C = 630;
n_min_4 = 1 + log(pc.ksi_c * C/(3*sqrt(2)*pc.gamma))/log(1/test_station)
