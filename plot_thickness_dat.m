clf(f3);
t = pc.dt:pc.dt:physical_time;
t = linspace(pc.dt,physical_time,length(loc_vec));
start_idx = 10000;
[A,B,t_0] =rootfit(t,loc_vec, start_idx);

f3 = figure(3);

plot(t(start_idx:end),A*sqrt(t(start_idx:end) - t_0) + B - loc_vec(start_idx));
hold on;

plot(t(start_idx:end),loc_vec(start_idx:end)-loc_vec(start_idx));
title("Ice interface location")
xlabel('time');
ylabel('location')
legend('analytical fit','numerical solution')
plotfixer;