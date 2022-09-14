function E = E_function(k_s, T_s, alpha, L, k_l, T_l)
%% This function calculates the error when determing what alpha is. This 
% lets the stefan_solver script converge on an alpha value. :)
E = 1/(sqrt(pi) * L)*(sqrt(k_s)*T_s * exp(-alpha^2/k_s)/erfc(alpha/sqrt(k_s)) + ...
    sqrt(k_l)*T_l * exp(-alpha^2/k_l)/(2-erfc(alpha/sqrt(k_l))))-alpha;
end