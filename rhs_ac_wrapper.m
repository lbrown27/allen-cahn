function f = rhs_ac_wrapper(c,T,u,eta,rho, x_coll, pc,adv)
if strcmp(pc.phase_model,"allen-cahn")
    [f,Advection,Diffusion,Sharpening,Phase_Change] = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);
elseif strcmp(pc.phase_model, "cdi")
    [f,Advection,Diffusion,Sharpening,Phase_Change]  = rhs_ac_CDI(c,T,u,pc,adv);
elseif strcmp(pc.phase_model, "acdi")
    [f,Advection,Diffusion,Sharpening,Phase_Change]  = rhs_ac_CDI(c,T,u,pc,adv);
else
    error("Error: phase model specified is not valid. \n");
end


%% DEBUGGING MODE:

[f_ac,Advection_ac,Diffusion_ac,Sharpening_ac,Phase_Change_ac] = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);
[f_cdi,Advection_cdi,Diffusion_cdi,Sharpening_cdi,Phase_Change_cdi]  = rhs_ac_CDI(c,T,u,pc,adv);
[f_acdi,Advection_acdi,Diffusion_acdi,Sharpening_acdi,Phase_Change_acdi]  = rhs_ACDI(c,T,u,pc,adv);

debug = 0;
if debug == 1
figure;
subplot(2,2,1);
plot(x_coll, f_ac);
title(append(pc.phase_model , " f"));
hold on;
xline(pc.x_init);

subplot(2,2,2);
plot(x_coll,Diffusion)
hold on;
xline(pc.x_init);
title(append(pc.phase_model , " Diffusion"));

subplot(2,2,3);
plot(x_coll,Sharpening);
hold on;
xline(pc.x_init);
title(append(pc.phase_model , " Sharpening"));

subplot(2,2,4);
plot(x_coll, Phase_Change);
title(append(pc.phase_model , " Phase Change"));
hold on;
xline(pc.x_init);
end
end