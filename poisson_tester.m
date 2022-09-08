%% Poisson solver tester

%RHS_Pres = -pi^2 / pc.l ^2 * cos(pi * x_coll / pc.l);
P_new = matrix_solve(P,RHS_Pres,pc,A_pres);
%exact_soln = cos(pi * x_coll/pc.l)+1;
%figure(3);
%plot(x_coll, P_new);
%hold on;
%plot(x_coll, exact_soln);
%diff_P_analytical = exact_soln - P_new;
%max_error_poisson = max(P_new - exact_soln);
%rms_error_poisson = sqrt(sum(diff_P_analytical .^2)/pc.N)
%P_new = gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);
%P_new_acc = accelerated_gauss_seidel(P,RHS_PE(rho_new, rho_n, rho_old, pc, u_star),pc);