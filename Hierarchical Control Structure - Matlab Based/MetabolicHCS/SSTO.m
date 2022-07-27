function [Us, Ys, SSTOflag] = SSTO(U, Pn, u_C4b, plant_ind, Us, Ys, X0_C4b, Y_MPC, PARA_C4b, SSTOflag)
% Order variables to otimize: u, y
% Units: u = m3/h (but enters L/min and leaves L/min)
%        y = mol/m3
inputs_C4b
lb = [0 18/100*40.88];
ub = [1 24/100*40.88];

if plant_ind == 1
    d = 0;
else
    d = -X0_C4b(PARA_C4b.i_O2_C4b)+Y_MPC(2,plant_ind-1);
end

A = [];
b = [];
Aeq = [];
beq = [];

objfun = @(x) (0)^2;
nlcon = @(x)metSSTOnlcon(x, u_C4b, Pn, X0_C4b, Y_MPC, PARA_C4b, plant_ind, U);

x0 = [u_C4b(in_gas_flow_C4b) U(PARA_C4b.index,4)/100/PARA_C4b.Air_Con];

[x,fval,out,text] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,nlcon);
Us(plant_ind,1) = x(1);
Ys(plant_ind,1) = x(2);
SSTOflag(plant_ind) = out;

end