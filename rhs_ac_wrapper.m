function f = rhs_ac_wrapper(c,T,u,eta,rho, pc,adv)
if strcmp(pc.phase_model,"allen-cahn")
    f = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);
elseif strcmp(pc.phase_model, "acdi")
    f = rhs_ac_ACDI(c,T,u,eta,rho, pc,adv);
else
    error("Error: phase model specified is not valid. \n");
end
end