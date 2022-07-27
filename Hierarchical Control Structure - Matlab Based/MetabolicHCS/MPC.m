function [U, Y, flow] = MPC(U1,Pn, u_C4b, plant_ind, U, Y, PARA_C4b, X0_C4b, us, ys, Par)
% Order variables to otimize: u, y
% Units: u = m3/h (but enters L/min and leaves L/min)
%        y = mol/m3
if plant_ind == 1
    x0 = zeros(Par.Nc,1);
else
    x0 = [U(:,plant_ind-1)];
end

inputs_C4b
lb = zeros(Par.Nc,1);
ub = ones(Par.Nc,1);

A = [];
b = [];
Aeq = [];
beq = [];

options = optimoptions('fmincon','OptimalityTolerance',1e-12,'algorithm','active-set');
nlcon  = @(x) metMPCnlcon(x, u_C4b, Pn, PARA_C4b, lb, ub, us, ys, X0_C4b, Par, plant_ind, Y, U1)
objfun = @(x) metMPCobj(x, u_C4b, Pn, PARA_C4b, lb, ub, us, ys, X0_C4b, Par, plant_ind, Y, U1);
                        
[x, fval, exitflag, output] = fmincon(objfun,x0,A,b,Aeq,beq,lb,ub,[],options);
[~,y] = metMPCobj(x, u_C4b, Pn, PARA_C4b, lb, ub, us, ys, X0_C4b, Par, plant_ind, Y, U1);

U(:,plant_ind) = x(1);
Y(:,plant_ind) = y;
flow = x(1);

end