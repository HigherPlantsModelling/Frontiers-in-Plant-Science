function [total_Err,y] = metMPCobj(x, u_C4b, Pn, PARA_C4b, lb, ub, us, ys, X0_C4b, Par, plant_ind, Y, U1) 
inputs_C4b;
y = zeros(Par.Np,1);
y(1) = X0_C4b(PARA_C4b.i_O2_C4b);
Err = 0;
xk_prev = u_C4b(in_gas_flow_C4b);
Time = [PARA_C4b.sample_time PARA_C4b.sample_time*2 PARA_C4b.sample_time*5 PARA_C4b.sample_time*4 PARA_C4b.sample_time*6 PARA_C4b.sample_time*6];
if plant_ind == 1
    d = 0;
else
    d = X0_C4b(PARA_C4b.i_O2_C4b)-Y(2,plant_ind-1);
end
for k=1:Par.Np
    spanU = 0.2/1000*60;
    spanY = 0.1/100/PARA_C4b.Air_Con;
    y(k+1) = y(k) + x(k)*U1(PARA_C4b.index,5)*(u_C4b(in_gas_o2_C4b)-y(k))*Time(k)/3600/PARA_C4b.Vchamber+Pn(1)*Time(k)/3600/PARA_C4b.Vchamber/32;
    %y(k+1) = (y(k)+Time(k)/3600/PARA_C4b.Vchamber*(x(k)*u_C4b(in_gas_o2_C4b)+Pn(1)*U1(PARA_C4b.index,5)/32))/(1+x(k)*Time(k)/3600/PARA_C4b.Vchamber);
    sortida = y(k+1);
    tmp1 =((x(k)-us(plant_ind,1))/spanU)^2;
    tmp2 =((sortida-ys(plant_ind,1))/spanY)^2;
    tmp3 =((x(k)-xk_prev)/spanU)^2;
    xk_prev = x(k);
    %Err=Err+*tmp1+100000*tmp2+1000*tmp3;   
    Err=Err+10*tmp2+1*tmp3;   

end    
%Err = Err + ((y(k+1)-ys(plant_ind,1))/spanY)^2;
%Y = [];
total_Err=Err;
% if nargout == 2
%     y = y;
% end