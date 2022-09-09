function dE_dalpha = E_deriv(k_s, T_s, alpha, L, k_l, T_l)
%% This function returns the dError_dalpha value (determined analytically
% with help from wolfram alpha) for use in a newton iteration loop in the
% script stefan_solver.
% See https://www.wolframalpha.com/input?i=derivative+of+exp%28-alpha%5E2%2Fk%29%2F%282-erfc%28alpha%2Fsqrt%28k%29%29%29+with+respect+to+alpha
% https://www.wolframalpha.com/input?i=derivative+of+exp%28-x%5E2%2FB%5E2%29%2Ferfc%28x%2FB%29
t_1 = sqrt(k_s)*T_s/(sqrt(pi)*L)*...
    (2 * exp(-2 * alpha^2/k_s)/ (sqrt(k_s *pi)*erfc(alpha / sqrt(k_s))^2));
t_1 = (1 - isnan(t_1))*t_1;
t_2 =  sqrt(k_s)*T_s/(sqrt(pi)*L)*...
    (- 2*alpha*exp(-alpha^2 / k_s)/(k_s * erfc(alpha / sqrt(k_s))));
t_1 = (1 - isnan(t_2))*t_2;
t_3 = sqrt(k_l)*T_l/(sqrt(pi)*L)*...
    (-2*alpha*exp(-alpha^2/k_l)/ (k_l * (2 - erfc(alpha / sqrt(k_l)))^2));
t_1 = (1 - isnan(t_3))*t_3;
t_4 =  sqrt(k_l)*T_l/(sqrt(pi)*L)*...
    (- 2*exp(-2 * alpha^2 / k_l)/(sqrt(k_l*pi) *(2-erfc(alpha / sqrt(k_l)))));
t_1 = (1 - isnan(t_4))*t_4;
t_5 = -1;

dE_dalpha = t_1 + t_2 + t_3 + t_4 + t_5;
end