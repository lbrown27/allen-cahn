function P = gauss_seidel_competitor(P_old,RHS,pc)
A = gallery('tridiag', pc.N,1,-2,1) / pc.dx^2;
A(1,1) = -1;
A(end,end) = 0;
RHS = transpose(RHS);
P = inv(A)*RHS;
end