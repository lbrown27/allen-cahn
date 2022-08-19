function c_next = c_step(c, T,u,c_past, pc)
% Steps the phase field variable forward in time.
for i = 1:pc.N
    rhs_cn = rhs_ac(c,T,u,pc);
    rhs_cn_past = rhs_ac(c,T,u,pc);
    c_next(i) = (4 * pc.dt * rhs_cn(i) - 2 * pc.dt * rhs_cn_past(i) + 4 * c(i) - c_past(i))/3;
end
c_next(1) = -1;
end