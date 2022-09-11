%% Better root fit test
clc;
A_actual = 2;
x0_actual = -3;
B_actual = -3;
  t = transpose(linspace(x0_actual,100+x0_actual,1000));
   dat = A_actual*sqrt(t-x0_actual)+ B_actual;

 start_idx = 1;
 
%  x1 = t(start_idx);
%  y1 = dat(start_idx);
%  x2 = t(floor(N/2));
%  y2 = dat(floor(N/2));
%  x3 = t(N);
%  y3 = dat(N);
%  x0_guess = -3.0;
 %[x0,A,B] = better_root_fit(x1,x2,x3,y1,y2,y3,x0_guess)
 
 %clf(4);
 %figure(4);
 
 %plot(t,dat);
 
 
 %% Actual problem
 f3 = figure(3);
clf(f3);

N = 15000;

scale = 10^3;
t = transpose(linspace(pc.dt,physical_time,length(loc_vec)));
dat = transpose(loc_vec)*scale;
%dat = sqrt(t);
plot(t(start_idx:end),dat(start_idx:end)-dat(start_idx));
title("Ice interface location")
xlabel('time');
ylabel('location')
hold on;

 f=fit(dat,t,'poly2');
 scale = 1000;
 p1_real = f.p1 / scale;
 p2_real = f.p2 / scale;
 p3_real = f.p3/ scale;
A = sqrt(1/p1_real);
B = -.5 * A^2*p2_real;
x0 = p3_real - B^2/A^2;

x_series = transpose(linspace(x0,t(end),1000));
ana_series = A*sqrt(x_series - x0); % REMOVED THE +B
plot(x_series, ana_series);
fprintf("A: %d,B:%d,x0:%d\n",A,B,x0);
 legend('numerical','analytical solution')

 
plotfixer;
 
