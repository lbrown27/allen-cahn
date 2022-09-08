function loc = find_interface_loc(c, x_coll,pc)
found = 0;
for i = 2:length(c)
    if (c(i) > .5) && found == 0
        loc = (-.5 - c(i - 1)) * pc.dx / (c(i) - c(i - 1)) + x_coll(i - 1);
        found = 1;
    end
end
if(found == 0)
    loc = inf;
end
end