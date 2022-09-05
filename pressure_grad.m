function grad_P = pressure_grad(P,pc)

grad_P = ones(pc.N+2,1) * 9999;
for i = 2: pc.N + 1
   grad_P(i) = (P(i+1) - P(i ))/(pc.dx);
end

grad_P(1) = 0;
% grad_P(pc.N+2) = 0; 

end