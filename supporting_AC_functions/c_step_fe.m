function c_next = c_step_fe(c, T,u,eta, rho, pc)
% Steps the phase field variable forward in time.
rhs_cn = rhs_ac_arezoo(c,T,u,eta,rho,pc,1);

rhs_cn = rhs_ac_arezoo(c,T,u,eta,rho,pc,1);

c_next = rhs_cn * pc.dt + c;

if strcmp(pc.left_BC,'Dirichlet')
    c_next(1) = - c_next(2);
elseif strcmp(pc.left_BC,'Neumann')
    c_next(1) = c_next(2);
else
    error("Error: Boundary condition not supported.\n")
end

c_next(pc.N + 2) = c_next(pc.N + 1);
end