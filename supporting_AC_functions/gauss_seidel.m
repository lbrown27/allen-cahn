function P = gauss_seidel(P_old,RHS,pc)
% performs gauss seidel to solve the pressure poisson equation
error = 1e16;
P = zeros(1,pc.N);
conv_val = 1e-16 * pc.N;
num_its = 0;
while error > conv_val
    num_its = num_its + 1;
    error = 0;
    for j = 2:pc.N - 1
        P(j) = - .5 * (pc.dx ^2 * RHS(j) - P_old(j + 1) - P(j - 1));
        error = error + abs(P(j) - P_old(j));
    end
    P(1) = P(2);
    P(pc.N) = 0;
    P_old = P;
end
fprintf("converged in %d iterations \n", num_its);

end