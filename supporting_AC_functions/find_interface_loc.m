function [loc,idx,dist_from_left_gridpoint] = find_interface_loc(c, x_coll,pc)
%% Finds the two gridpoints at which c=.5 is located and then performs a
% linear interpolation to find interface location.
% The index returned is the index of the point to the RIGHT of the
% interface!

found = 0;
for i = 2:length(c)
    if (c(i) > .5) && found == 0
        loc = (.5 - c(i - 1)) * pc.dx / (c(i) - c(i - 1)) + x_coll(i - 1);
        found = 1;
        idx = i;
        dist_from_left_gridpoint = loc - x_coll(i - 1);
    end
end
if(found == 0)
    loc = inf;
    idx = inf;
    dist_from_left_gridpoint = inf;
end
end