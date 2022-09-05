function E = E_matrix(rho_cp, u, k, pc)

 diff_E = zeros(pc.N+2,pc.N+2);
adv_E = zeros(pc.N+2,pc.N+2);
E = sparse(zeros(pc.N + 2, pc.N +2));
 %diff_E = sparse(zeros(pc.N+2,pc.N+2));
 %adv_E = sparse(zeros(pc.N+2,pc.N+2));
for i = 2:pc.N +1
    i_plus = i;
i_minus = i - 1; % indices for flux terms
% diff_E(i,i+1) = (k(i + 1) + k(i)) / (2 * pc.dx ^2);
% diff_E(i,i) = - (k(i + 1) + 2 * k(i) + k(i - 1)) / (2 * pc.dx ^2);
% diff_E(i,i-1) = (k(i) + k(i - 1)) / (2 * pc.dx ^2);
% 
% adv_E(i,i+1) = (-rho_cp(i+1) * u(i_plus))/(2 * pc.dx);
% adv_E(i,i) = (-rho_cp(i) * u(i_plus) + rho_cp(i) * u(i_minus))/(2 * pc.dx);% there was a minus here!
% adv_E(i,i-1) = (rho_cp(i-1) * u(i_minus))/(2 * pc.dx);

   E(i,i + 1) = (-rho_cp(i+1) * u(i_plus))/(2 * pc.dx) + (k(i + 1) + k(i)) / (2 * pc.dx ^2);
   E(i,i) = (-rho_cp(i) * u(i_plus) + rho_cp(i) * u(i_minus))/(2 * pc.dx) - (k(i + 1) + 2 * k(i) + k(i - 1)) / (2 * pc.dx ^2);
   E(i,i-1) =  (rho_cp(i-1) * u(i_minus))/(2 * pc.dx) + (k(i) + k(i - 1)) / (2 * pc.dx ^2);

   
  
end

%E = diff_E + adv_E;


%E(1,1) = rho_cp(1) * u(1)/ (2 * pc.dx) - (k(2) + k(1)) / (2 * pc.dx^2);
%E(1,2) = rho_cp(2) * u(1)/ (2 * pc.dx) + (k(2) + k(1)) / (2 * pc.dx^2);






% for i = 3:pc.N - 1
%     A = ((rho_cp_T(i+1) + rho_cp_T(i)) * u(i_plus) - (rho_cp_T(i) + rho_cp_T(i - 1)) * u(i_minus))/(2 * pc.dx);
%     B = ((k(i+1) + k(i)) * (T(i+1) + T(i)) - (k(i) + k(i-1)) * (T(i) - T(i - 1))) / (2 * pc.dx ^ 2);
%     E(i) = A + B;
% end

% 
end