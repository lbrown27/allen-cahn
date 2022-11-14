function experimental_function(c,pc)
stag_c = stagger(c,pc);
S = div_stag(stag_c.^2 .* (1 - stag_c).^2,pc); % CDI term. ADDED A FACTOR OF FOUR.
S_e = (4 * stag_c.^3 -6*stag_c.^2+2*stag_c).*div_stag(stag_c,pc);
diff = S - S_e;

Sharpening = div_stag(pc.gas_pedal * -1 * 4 * stag_c.^2 .* (1 - stag_c).^2,pc); % CDI term. ADDED A FACTOR OF FOUR.
Sharpening_experimental = pc.gas_pedal * (-4) * (4 * stag_c.^3 -6*stag_c.^2+2*stag_c).*div_stag(stag_c,pc);
diff_2 = Sharpening - Sharpening_experimental;
end