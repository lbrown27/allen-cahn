function f = rhs_ac_for_implicit_T_wrapper(c,u,T,eta_new,rho_new, pc,adv)
if (strcmp(pc.phase_model,"allen-cahn"))
    f = rhs_ac_arezoo_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv);
elseif (strcmp(pc.phase_model,"acdi"))
    f = rhs_ac_ACDI_for_implicit_T(c,u,T,eta_new,rho_new, pc,adv);
else
    error("Phase model specified is not valid. \n");
end
end