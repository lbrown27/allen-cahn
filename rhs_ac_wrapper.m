function f = rhs_ac_wrapper(c,T,u,eta,rho, x_coll, pc,adv)

%% I made the code more verbose so that the acdi and A-C methods can be directly compared.
fprintf("Difference between f_ac - f_acdi: \n");

    f_ac = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);

    f_acdi = rhs_ac_ACDI(c,T,u,eta,rho, x_coll,pc,adv);
diff = f_ac - f_acdi

    
if strcmp(pc.phase_model,"allen-cahn")
    f = f_ac;
elseif strcmp(pc.phase_model, "acdi")
    f = f_acdi;
else
    error("Error: phase model specified is not valid. \n");
end
end