function percent_lost = mass_counter(rho_new, rho_n, pc, c_new, c_n,mass_orig, x_coll)
dl = find_interface_loc(c_new, x_coll,pc) - find_interface_loc(c_n, x_coll,pc);
generated_mass = dl * (pc.rho_ice - pc.rho_water);
domain_mass_change = (sum(rho_new(2:end - 1)) - sum(rho_n(2:end - 1))) * pc.dx;
mass_lost =  generated_mass - domain_mass_change;
percent_lost = mass_lost / mass_orig;

end