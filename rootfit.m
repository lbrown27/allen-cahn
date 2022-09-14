function [A,B,t_0] =rootfit(t,dat, start_idx)
t_0 = t(start_idx);
t = t-t(start_idx);
B = (dat(end) - dat(start_idx + 1)*sqrt(t(end)/t(start_idx+1)))/(1-sqrt(t(end)/t(start_idx+1)));
A = (dat(start_idx+1)-B)/sqrt(t(start_idx+1));
end