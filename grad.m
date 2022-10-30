function g = grad(coll,pc)
%% Returns a STAGGERED quantity which is the gradient of a collocated field 
% variable. g_(i + 1/2) in real life = g(i) in matlab.
g = ones(pc.N+2,1) * 9999;
for i = 1: pc.N + 1
   g(i) = (coll(i+1) - coll(i))/(pc.dx);
end

g(end) = 0;
% grad_P(pc.N+2) = 0; 

end