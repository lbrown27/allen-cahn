%% Plot result
figure(3);
plot(t_vec,loc_num);
hold on;
plot(t_vec, loc_ana);
title("Numerical vs Analytical Stefan Problem");
legend("Numerical","Analytical");
hold off;

figure(4);
plot(t_vec, (loc_num - loc_ana)/abs(loc_num(1) - loc_num(end)));
title("Error between Numerical & Analytical Result");
hold on;
yline(0);
hold off;
plotfixer;

error = abs(loc_num(end) - loc_ana(end))/abs(loc_num(1) - loc_num(end));

frog = (loc_num(end) - loc_ana(end))/ pc.ksi_c;