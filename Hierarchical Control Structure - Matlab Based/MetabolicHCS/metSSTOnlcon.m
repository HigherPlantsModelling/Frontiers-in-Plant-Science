function [c,ceq] = metSSTOnlcon(x, u_C4b, Pn, X0_C4b, Y, PARA_C4b, plant_ind, U)
inputs_C4b
c = [];
if plant_ind == 1
    d = 0;
else
    d = X0_C4b(PARA_C4b.i_O2_C4b)-Y(2,plant_ind-1);
end

ceq = [x(1)*U(PARA_C4b.index,5)*(u_C4b(in_gas_o2_C4b)-x(2))+Pn(1)/32; x(2)-U(plant_ind,4)/100/PARA_C4b.Air_Con];
end