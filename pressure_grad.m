function grad_P = pressure_grad(P,pc)

grad_P = ones(1,pc.N) * 9999;
for i = 2: pc.N - 1
   grad_P(i) = (P(i + 1) - P(i - 1))/(2 * pc.dx);
end

grad_P(1) = 0;
grad_P(pc.N) = 0; % this would normally not be the case. However, we don't
% care what grad P is here because later we will enforce that u(N) = 0. As 
% a reminder, the reason why we are calculating grad P is so that we can 
% find the corrected velocity.

end