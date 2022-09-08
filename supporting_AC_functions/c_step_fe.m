function c_next = c_step_fe(c, T,u,eta, rho, pc)
% Steps the phase field variable forward in time.
rhs_cn = rhs_ac(c,T,u,eta,rho,pc);
c_next = rhs_cn * pc.dt + c;
c_next(1) = - c_next(2);
c_next(pc.N + 2) = c_next(pc.N + 1);
end