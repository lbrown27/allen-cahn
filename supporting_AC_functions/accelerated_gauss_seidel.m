function P = accelerated_gauss_seidel(P_old,RHS,pc)
% performs gauss seidel to solve the pressure poisson equation
error = 1e16;
P = zeros(1,pc.N);
conv_val = 1e-16 * pc.N;
count = 0;
while error > conv_val
    error = 0;
    if (mod(count, 2) == 0)
        for j = 2:pc.N - 1
            P(j) = - .5 * (pc.dx ^2 * RHS(j) - P_old(j + 1) - P(j - 1));
            error = error + abs(P(j) - P_old(j));
        end
        P(1) = P(2);
        P(pc.N) = 0;
        %fprintf("frog");
    else
        for j = 2:pc.N - 1
            P_old(j) = - .5 * (pc.dx ^2 * RHS(j) - P(j + 1) - P_old(j - 1));
            error = error + abs(P_old(j) - P(j));
        end
        P_old(1) = P_old(2);
        P_old(pc.N) = 0; % now, P_old is ahead of P.
    end
    count = count + 1;
end
if mod(count, 2) == 0
   P = P_old; % if count is even, this means the latest soln is in P_old, not P. 
end
fprintf("accelerated converged in %d iterations \n", count);
end