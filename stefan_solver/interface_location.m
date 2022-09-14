function loc = interface_location(start_loc, alpha,t,L_domain)
%% returns the interface location for stefan problem 
% of a certain alpha.

% PARAMS:
% alpha: alpha for the given stefan problem, see stefan_solver.m
% start_loc: x location where the ice starts (t = 0)
% t: current time
loc = L_domain - start_loc -  2 * alpha * sqrt(t);
end