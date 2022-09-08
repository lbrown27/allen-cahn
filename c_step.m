function c_next = c_step(c, T,u,c_past,eta, eta_old,rho, pc)
% Steps the phase field variable forward in time.
rhs_cn = rhs_ac(c,T,u,eta,rho,pc);
rhs_cn_past = rhs_ac(c,T,u,eta_old,rho,pc);
c_next = (4 * pc.dt * rhs_cn - 2 * pc.dt * rhs_cn_past + 4 * c - c_past)/3;
c_next(1) = - c_next(2);
c_next(pc.N + 2) = c_next(pc.N + 1);
end