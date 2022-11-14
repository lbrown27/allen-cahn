function stag = stagger(coll,num_pts)
%% Staggers a collocated field. stag_(i + 1/2) in real life = stag(i) in matlab.

% Applies an outflow boundary condition.
stag = zeros(num_pts,1);
for i = 1: num_pts - 1
   stag(i) = (coll(i+1) + coll(i))/2;
end
stag(num_pts) = stag(num_pts - 1); % outflow bc.
end