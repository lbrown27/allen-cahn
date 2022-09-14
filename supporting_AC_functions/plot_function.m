function plot_function(x_coll, x_stag,pc,c_new, T_new,t_vec, loc_ana,loc_num,count,physical_time)

fprintf("Time step: %d \n Iteration: %d \n",pc.dt,count);
fprintf("Physical time: %d \n",physical_time);
f= figure(1);
n_plots = 3;
plot_count = 1;

subplot(n_plots,1,plot_count);
plot(x_coll,c_new); % plot phase field after one step
title("phase field");
xlim([x_coll(1), x_stag(end)])
ylim([0 1])
plot_count = plot_count + 1;

subplot(n_plots,1,plot_count);
plot(x_coll,T_new);
title('Temp');
xlim([x_coll(1), x_stag(end)])
ylim([pc.wall_T - 2 pc.init_T])
plot_count = plot_count + 1;



%         subplot(n_plots,1,3);
%         plot(x_stag,u_new);
%         title("velocity");
%         xlim([x_coll(1), x_stag(end)])
%         %ylim([0 1e-4])
%
%         subplot(n_plots,1,4);
%         plot(x_coll,P_new);
%         title('Pressure');
%         xlim([x_coll(1), x_stag(end)])
%         %ylim([-1 10])
%
%
subplot(n_plots,1,plot_count);
plot(t_vec,loc_num);
hold on;
plot(t_vec, loc_ana);
hold off;
legend('num','ana');
title('Interface Location');
%
%
%
%         subplot(n_plots,1,6);
%         hist = 4999;
%         plot(t_vec(count - hist:end),loc_vec(count - hist:end));
%         hold on;
%         [A,B] = rootfit(t_vec(count - hist:end),loc_vec(count - hist:end));
%         plot(t_vec(count - hist:end), A*sqrt(t_vec(count - hist:end)) + B);
%         hold off;
%         legend('num','ana');
%         title('Fit to last 1000 pts');

%xlim([x_coll(1), x_stag(end)])
%ylim([0 pc.l])
%       subplot(n_plots,1,5);
%         plot(x_stag,u_star);
%         title('U Star');
%         xlim([x_coll(1), x_stag(end)])
%
%         subplot(n_plots,1,6);
%         plot(x_coll, c_diff);
%         title('c diff');
%         xlim([x_coll(1), x_stag(end)])
%max(abs(u_diff))




%         subplot(n_plots,1,7);
%         plot(x_stag,RU_n_func(rho_new, eta_new, pc, c_new,u_new));
%         title('RUn function');
%         xlim([x_coll(1), x_stag(end)])
drawnow();
fprintf("graphs updated. \n");

end