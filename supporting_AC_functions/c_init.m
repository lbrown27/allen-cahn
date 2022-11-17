function c = c_init(c,x,pc)
if pc.phase_model == 'allen-cahn'
    thick_factor = sqrt(2);
else
    thick_factor = 2;
end
%c = (tanh((x - pc.l / 2)/(sqrt(2) * pc.ksi_c)) - 1)/2;
for i = 2:pc.N + 1
    %if x(i) < pc.ksi_c
    % c(i) = 1 / pc.ksi_c * x(i) - 1;
        c(i) = (tanh((x(i) - pc.x_init)/(pc.ksi_c * thick_factor)) +1)/2;
    %end
    
end

end