function print_coefficients(pc)
Data = [  pc.gamma*pc.mu/pc.ksi_c^2  pc.gamma*pc.mu/pc.ksi_c^2     pc.mu / (3*sqrt(2)*pc.ksi_c)
          pc.gas_pedal / pc.ksi_c        pc.gas_pedal / pc.ksi_c pc.gas_pedal / (pc.gamma*3*sqrt(2))];
VarNames = {'Diffusion', 'Sharpening', 'Phase Change'};
T = table(Data(:,1),Data(:,2),Data(:,3), 'VariableNames',VarNames)

end