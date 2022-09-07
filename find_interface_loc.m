function loc = find_interface_loc(c, x_coll,pc)
count = 1;
while c(count) < -.5
    count = count + 1;
end
loc = (-.5 - c(count - 1)) * pc.dx / (c(count) - c(count - 1)) + x_coll(count - 1);

end