function [x0,A,B] = better_root_fit(x1,x2,x3,y1,y2,y3,x0_guess)

fcn = @(x0) (y1*sqrt(x1-x0)-y1*sqrt(x3-x0)-y3*sqrt(x1-x0)+y1*sqrt(x3-x0))/(x1-x0-sqrt(x3-x0)*sqrt(x1-x0)) * sqrt(x2-x0) + ...
    (y3*sqrt(x1-x0)-y1*sqrt(x3-x0))/(sqrt(x1-x0)-sqrt(x3-x0))-y2;

x0 = fzero(fcn, x0_guess);
B = (y3 * sqrt(x1-x0)-y1*sqrt(x3-x0))/(sqrt(x1-x0)-sqrt(x3-x0));
A = (y1 - B)/sqrt(x1-x0);

end