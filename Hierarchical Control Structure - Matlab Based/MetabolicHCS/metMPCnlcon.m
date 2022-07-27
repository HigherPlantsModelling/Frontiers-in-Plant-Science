function [c,ceq] = metMPCnlcon(x, u_C4b, Pn, X0_C4b, Y, PARA_C4b, plant_ind, U)
inputs_C4b
c = [];
y = zeros(Par.Np,1);
if plant_ind == 1
    d = 0;
else
    d = X0_C4b(PARA_C4b.i_O2_C4b)-Y(2,plant_ind-1);
end
y(1) = X0_C4b(PARA_C4b.i_O2_C4b);
y(2) = y(1) + x(1)*(u_C4b(in_gas_o2_C4b)-y(1))*Time(1)/3600/PARA_C4b.Vchamber+Pn(1)*U1(PARA_C4b.index,5)*Time(1)/3600/PARA_C4b.Vchamber/32;
end