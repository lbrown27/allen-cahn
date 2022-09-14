function dF = simple_jacobian_calculation(u_n,rho_new,eta_new,pc,c_new)
%% This function calculates the jacobian simply.
%  method for implicit time advancement
dF = sparse(zeros(pc.N + 2));

for i = 2: pc.N + 1
    dF(i,i) = (rho_new(i+1)*(u_n(i+1)+u_n(i)) - rho_new(i+1)*(u_n(i) + u_n(i-1)))/(2* pc.dx)- 4/(3 * pc.dx^2) * (eta_new(i+1) + eta_new(i));
    dF(i,i-1) = -rho_new(i) * (u_n(i) + u_n(i -1)) / (2*pc.dx) + 4/3 * (eta_new(i)/pc.dx ^ 2);
    dF(i,i+1) = rho_new(i+1) * (u_n(i+1) + u_n(i)) / (2*pc.dx) + 4/3 * (eta_new(i+1)/pc.dx ^ 2);
    
end

    dF(1,1) = (rho_new(2)*(u_n(2)+u_n(1)) - rho_new(2)*(u_n(1)))/(2* pc.dx)- 4/(3 * pc.dx^2) * (eta_new(2) + eta_new(1));
dF(1,2) = rho_new(2) * (u_n(2) + u_n(1)) / (2*pc.dx) + 4/3 * (eta_new(2)/pc.dx ^ 2);

dF(pc.N + 2, pc.N +2) = (rho_new(pc.N + 2)*(u_n(pc.N + 2)*2) - rho_new(pc.N +2)*(u_n(pc.N +2) + u_n(pc.N +1)))/(2* pc.dx)- 4/(3 * pc.dx^2) * (eta_new(pc.N +2));
dF(pc.N + 2, pc.N +1) = -rho_new(pc.N +2) * (u_n(pc.N +2) + u_n(pc.N +1)) / (2*pc.dx) + 4/3 * (eta_new(pc.N +2)/pc.dx ^ 2);
end