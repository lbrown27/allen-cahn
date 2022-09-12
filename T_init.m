function T = T_init(c,x,pc)
%c = (tanh((x - pc.l / 2)/(sqrt(2) * pc.ksi_c)) - 1)/2;
T = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    %if x(i) < pc.ksi_c
       % c(i) = 1 / pc.ksi_c * x(i) - 1;
        T(i) = (pc.init_T - pc.wall_T)*(tanh((x(i) - pc.l/2)/(pc.ksi_c * 2)) +1)/2 + pc.wall_T;
    %end
    
end

end