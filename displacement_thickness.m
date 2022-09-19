function integral_thickness_water_side = displacement_thickness(c,x_coll,pc)
%% Calculates the interface thickness assuming that water is on the right
% and ice is on the left.
[int_loc, int_idx,dist_from_left_gridpt] = find_interface_loc(c, x_coll,pc);

c =c - .5;

%% Water Side:
integral = (pc.dx - dist_from_left_gridpt) * (c(int_idx))/2;
for i = int_idx:length(c) - 1
   integral = integral + pc.dx/ 2 * (c(i)+c(i+1)); 
end

end_of_measuring = (x_coll(end)+ x_coll(end-1))/2;
integral_thickness_water_side = (end_of_measuring-int_loc) - 2 * integral;


%% Ice Side:
integral_l = dist_from_left_gridpt * (c(int_idx))/2;
for i = 1:int_idx - 1
   integral_l = integral_l + pc.dx/ 2 * (c(i)+c(i+1)); 
end

beginning_of_measuring = (x_coll(1)+ x_coll(2))/2;
integral_thickness_water_side = (end_of_measuring-beginning_of_measuring) +2 * integral_l - 2*integral;
end