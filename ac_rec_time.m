function recommended_timestep = ac_rec_time(u,pc, c,T)
k = 1;
lambda = zeros(pc.N + 2,1);
sigma_plus = zeros(pc.N + 2,1);
sigma_minus = zeros(pc.N + 2,1);
sigma_final_worst_case = 0;
dt = pc.dt;
im = sqrt(-1);
num = 100;
max_rk4_store_val = 0;
%while sigma_final_worst_case > 1
    
for j = 1:num
for i = 2:pc.N
    k = 2 *pi * j / (num * pc.dx);
i_plus = i;
i_minus = i - 1;
A = 6 * pc.sigma_c * pc.ksi_c;
part_1 = -(u(i_plus)*(exp(im * k * pc.dx) + 1) - u(i_minus) * (1 + exp(-im*k*pc.dx)))/(2*pc.dx);
    part_2 = pc.Mc* A*(2 * cos(k * pc.dx) - 2)/(pc.dx)^2; 
    part_3 = -(3 * pc.sigma_c / pc.ksi_c * (4 * c(i)^2 + 6 * c(i) + 2) + pc.L * pc.rho_water * (pc.T_M - T(i)) / pc.T_M * (30 * (c(i) + 1)^2 * c(i)));
lambda(i) = part_1 + part_2 + part_3;
    sigma_plus(i) = 1/6 * (4 + 4 *dt * lambda(i) + sqrt((4+ 4 * dt * lambda(i))^2 - 12 * (1 + 2 * dt * lambda(i))));
sigma_minus(i) = 1/6 * (4 + 4 *dt * lambda(i) + sqrt((4+4 * dt * lambda(i))^2 - 12 * (1 + 2 * dt * lambda(i))));

end
rk4 = 1 + lambda * dt + lambda .^2 * dt^2 / 2 + lambda .^3 *dt^3 / 6 + lambda .^4 * dt^4/ 24;

rk4 = abs(rk4);

max_rk_val = max(abs(rk4));
 
if max_rk_val > max_rk4_store_val
    max_rk4_store_val = max_rk_val;
end
    


max_val(j) = max(max(abs(sigma_plus), abs(sigma_minus)));

end
% NOW you can judge how good it is! Loops are done.
if (max(max_val) > sigma_final_worst_case)
sigma_final_worst_case = max(max_val);
end
recommended_timestep = 5;
end