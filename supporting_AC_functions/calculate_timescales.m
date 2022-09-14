function [t_advective,t_diffusive, t_func]  = calculate_timescales(pc, u,rho)
t_advective = pc.dx / max(abs(u));
t_diffusive = pc.dx ^ 2 / (pc.Mc * 6 * pc.sigma_c * pc.ksi_c);
t_func = pc.ksi_c / ((3 * pc.ksi_c + 30 * pc.L* max(rho) * pc.ksi_c)* pc.Mc);
end