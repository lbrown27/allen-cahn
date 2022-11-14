function psi = distance_field(c,x_coll,pc)
% returns a collocated field of distances from the interface location.
int_loc = find_interface_loc(c, x_coll,pc);

psi = x_coll - int_loc;
end