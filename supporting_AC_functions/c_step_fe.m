function c_next = c_step_fe(c, T,u,eta, rho, pc)
% Steps the phase field variable forward in time.
rhs_cn = rhs_ac_arezoo(c,T,u,eta,rho,pc,1);

rhs_cn = rhs_ac_arezoo(c,T,u,eta,rho,pc,1);

c_next = rhs_cn * pc.dt + c;
c_next(1) = - c_next(2);
c_next(pc.N + 2) = c_next(pc.N + 1);
end