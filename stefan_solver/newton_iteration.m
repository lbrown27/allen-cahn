function alpha_next = newton_iteration(k_s, T_s, alpha, L, k_l, T_l)
%% This calculates the next alpha guess for the stefan_solver script using
% a newton iteration.
E_f = E_function(k_s, T_s, alpha, L, k_l, T_l);
E_d = E_deriv(k_s, T_s, alpha, L, k_l, T_l);
alpha_next = alpha -  E_f/ E_d;
end