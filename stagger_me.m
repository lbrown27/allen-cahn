function stag = stagger_me(coll,pc)
%% Staggers a collocated field. stag_(i + 1/2) in real life = stag(i) in matlab.
% wrt a collocated point i, i_plus would be i and i_minus would be i - 1, 
% which is consistent with the rest of the code.
% Applies an outflow boundary condition.
stag = zeros(pc.N+2,1);
for i = 1: pc.N + 1
   stag(i) = (coll(i+1) + coll(i))/2;
end
stag(pc.N + 2) = stag(pc.N + 1); % outflow bc.
end