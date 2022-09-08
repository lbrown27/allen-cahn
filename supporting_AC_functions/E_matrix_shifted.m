function E = E_matrix_shifted(rho_cp, u, k, pc)

 %diff_E = zeros(pc.N+2,pc.N+2);
%adv_E = zeros(pc.N+2,pc.N+2);
%E = sparse(zeros(pc.N + 2, pc.N +2));
%E_old = sparse(zeros(pc.N + 2, pc.N +2));

 %diff_E = sparse(zeros(pc.N+2,pc.N+2));
 %adv_E = sparse(zeros(pc.N+2,pc.N+2));
rho_cp_minus =zeros(size(rho_cp));
rho_cp_minus(2:end)=rho_cp(1:end-1);

rho_cp_plus =zeros(size(rho_cp));
rho_cp_plus(1:end - 1)=rho_cp(2:end);

k_minus =zeros(size(k));
k_minus(2:end)=k(1:end-1);

k_plus =zeros(size(k));
k_plus(1:end - 1)=k(2:end);

u_minus = zeros(size(u));
u_minus(2:end)=u(1:end-1);

minus = (rho_cp_minus .* u_minus)/(2 * pc.dx) + (k + k_minus) / (2 * pc.dx ^2);
normal = (-rho_cp .* u + rho_cp .* u_minus)/(2 * pc.dx) - (k_plus + 2 .* k + k_minus) / (2 * pc.dx ^2);
plus = (-rho_cp_plus .* u)/(2 * pc.dx) + (k_plus + k) / (2 * pc.dx ^2);
E = gallery('tridiag',minus(2:end),normal,plus(1:end-1));
% 
% for i = 2:pc.N + 1
%     i_plus = i;
%     i_minus = i - 1;
%    E_old(i,i + 1) = (-rho_cp(i+1) * u(i_plus))/(2 * pc.dx) + (k(i + 1) + k(i)) / (2 * pc.dx ^2);
%    E_old(i,i) = (-rho_cp(i) * u(i_plus) + rho_cp(i) * u(i_minus))/(2 * pc.dx) - (k(i + 1) + 2 * k(i) + k(i - 1)) / (2 * pc.dx ^2);
%    E_old(i,i-1) =  (rho_cp(i-1) * u(i_minus))/(2 * pc.dx) + (k(i) + k(i - 1)) / (2 * pc.dx ^2);
% end
%E = full(E);
%diff = E - E_old;
%max_diff =  full(max(max(abs(diff))));
%fprintf("max diff is %d",max_diff);
   

end