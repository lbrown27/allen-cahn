function alpha = find_alpha(k_l, k_s,L,T_l,T_s,alpha_next)
%count = 0;
E_function(k_s, T_s, alpha_next, L, k_l, T_l)
while abs(E_function(k_s, T_s, alpha_next, L, k_l, T_l)) > eps(0)
    alpha_past = alpha_next;
    alpha_next =  newton_iteration(k_s, T_s, alpha_past, L, k_l, T_l);
    %count = count + 1;
end
alpha = alpha_next;
end
