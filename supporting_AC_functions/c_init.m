function c = c_init(c,x,pc)
%c = (tanh((x - pc.l / 2)/(sqrt(2) * pc.ksi_c)) - 1)/2;
for i = 2:pc.N + 1
    %if x(i) < pc.ksi_c
       % c(i) = 1 / pc.ksi_c * x(i) - 1;
        c(i) = (tanh((x(i) - pc.l/2)/(pc.ksi_c/sqrt(200))) +1)/2;
    %end
    
end

end