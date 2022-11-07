function c_next = c_step(c, T,u,c_past,eta, eta_old,rho,x_coll, pc)
% Steps the phase field variable forward in time.
%rhs_cn = rhs_ac(c,T,u,eta,rho,pc,1);
%rhs_cn_past = rhs_ac(c,T,u,eta_old,rho,pc,1);

rhs_cn = rhs_ac_wrapper(c,T,u,eta,rho,x_coll,pc,1);
rhs_cn_past = rhs_ac_wrapper(c,T,u,eta_old,rho,x_coll,pc,1);
c_next = (4 * pc.dt * rhs_cn - 2 * pc.dt * rhs_cn_past + 4 * c - c_past)/3;
if strcmp(pc.left_BC,'Dirichlet')
    c_next(1) = - c_next(2);
elseif strcmp(pc.left_BC,'Neumann')
    c_next(1) = c_next(2);
else
    error("Error: Boundary condition not supported.\n")
end
c_next(pc.N + 2) = c_next(pc.N + 1);
end