function d = div_stag(s,num_pts,dx)
%% Divergence of collocated and staggered variable. c & s are the collocated
% and staggered variables, respectively.
d = zeros(num_pts,1);
for i = 2:num_pts - 1
    i_plus = i;
    i_minus = i - 1;
    d(i) = (s(i_plus) - s(i_minus))/(dx); 
end
end