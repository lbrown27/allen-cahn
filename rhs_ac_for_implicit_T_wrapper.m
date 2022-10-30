function f = rhs_ac_for_implicit_T_wrapper(c,u,T,eta_new,rho_new, pc,adv)
f = rhs_ac_arezoo_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv);

%f = rhs_ac_ACDI_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv);
end