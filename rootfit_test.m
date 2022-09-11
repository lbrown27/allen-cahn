%% Rootfit test
 
 t = linspace(0,100,1000);
 dat = sqrt(t);
 start_idx = 1;
 
 [A,B,t_0] = rootfit(t,dat, start_idx);
 
 figure();
 fit = A*sqrt(t - t_0) +B;
 plot(t,fit);
 
 hold on;
 plot(t,dat);
 
 diff = dat - fit;
