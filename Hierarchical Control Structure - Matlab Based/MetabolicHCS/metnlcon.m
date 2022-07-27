function [c,ceq] = metnlcon(x, u_C4b, Pn)
inputs_C4b
c = [];
ceq = x(1)*(u_C4b(in_gas_o2_C4b)-x(2)+Pn(1)/32);
end