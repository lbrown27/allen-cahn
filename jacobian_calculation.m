function A = jacobian_calculation(u_n, u_new)
%% This function calculates the jacobian so that can use in the crank-nicolson
%  method for implicit time advancement
diff_u = u_new - u_n;
epsilon_factor = .01;
A = diff_u;



end