function f = rhs_ac_wrapper(c,T,u,eta,rho, pc,adv)
f = rhs_ac_arezoo(c,T,u,eta,rho, pc,adv);
%f = rhs_ac_ACDI(c,T,u,eta,rho, pc,adv);
end