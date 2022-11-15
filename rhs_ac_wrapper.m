function f = rhs_ac_wrapper(c,T,u,eta,rho, x_coll, pc,adv)
if strcmp(pc.phase_model,"allen-cahn")
    [f,Advection_ac,Diffusion_ac,Sharpening_ac,Phase_Change_ac] = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);
elseif strcmp(pc.phase_model, "cdi")
    [f,Advection_acdi,Diffusion_acdi,Sharpening_acdi,Phase_Change_acdi]  = rhs_ac_CDI(c,T,u,pc,adv);
elseif strcmp(pc.phase_model, "acdi")
    [f,Advection_acdi,Diffusion_acdi,Sharpening_acdi,Phase_Change_acdi]  = rhs_ac_CDI(c,T,u,pc,adv);
else
    error("Error: phase model specified is not valid. \n");
end
end