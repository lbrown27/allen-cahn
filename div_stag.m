function d = div_stag(s,pc)
%% Divergence of collocated and staggered variable. c & s are the collocated
% and staggered variables, respectively.
d = zeros(pc.N + 2,1);
for i = 2:pc.N + 1
    i_plus = i;
    i_minus = i - 1;
    d(i) = (s(i_plus) - s(i_minus))/(pc.dx); % ADVECTION TERM
end
end