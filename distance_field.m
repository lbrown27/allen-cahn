function psi = distance_field(c)
% returns a collocated field of distances from the interface location.
%int_loc = find_interface_loc(c, x_coll,pc);
small_num = 1e-16;

psi = x_small_num*log((c+small_num)/(1-c+small_num));
end