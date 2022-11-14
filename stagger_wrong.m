function stag = stagger_wrong(coll,pc)
%% Please delete me
%% Staggers a collocated field. stag_(i + 1/2) in real life = stag(i) in matlab.
% Applies an outflow boundary condition.
stag = zeros(pc.N+2,1);
for i = 2: pc.N + 2
   stag(i) = (coll(i) + coll(i-1))/2;
end
% stag(1) won't be used anyway, because new bcs for c will be applied.
end