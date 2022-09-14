function T_ana = stefan_temp_field(T_l, k_l,alpha, start_loc,t,T_s, k_s,x,L_domain)
int_loc = interface_location(start_loc, alpha,t,L_domain);

t_1 = -T_l * erfc(alpha / sqrt(k_l)) / (2-erfc(alpha/sqrt(k_l)));
t_2 = T_l * erfc(((L_domain - x) - start_loc)/(2*sqrt(k_l * t)))/(2 - erfc(alpha/sqrt(k_l)));
t_3 = T_s;
t_4 = -T_s * erfc(((L_domain - x) - start_loc)/(2*sqrt(k_s * t)))/(erfc(alpha/sqrt(k_s)));

T_ana = (t_1 + t_2).*((L_domain - x) < int_loc) + (t_3 + t_4).*(x >= int_loc);

end