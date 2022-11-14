function f = rhs_ac_wrapper(c,T,u,eta,rho, x_coll, pc,adv)

%% I made the code more verbose so that the acdi and A-C methods can be directly compared.
%fprintf("Difference between f_ac - f_acdi: \n");

    [f_ac,Advection_ac,Diffusion_ac,Sharpening_ac,Phase_Change_ac] = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);

   [f_acdi,Advection_acdi,Diffusion_acdi,Sharpening_acdi,Phase_Change_acdi]  = rhs_ac_ACDI(c,T,u,eta,rho, x_coll,pc,adv);
diff = f_ac - f_acdi;
Advection_diff = Advection_ac - Advection_acdi;
Diffusion_diff = Diffusion_ac - Diffusion_acdi;
Sharpening_diff = Sharpening_ac - Sharpening_acdi;
Phase_Change_diff = Phase_Change_ac - Phase_Change_acdi;
    
if strcmp(pc.phase_model,"allen-cahn")
    f = f_ac;
elseif strcmp(pc.phase_model, "acdi")
    f = f_acdi;
else
    error("Error: phase model specified is not valid. \n");
end
end