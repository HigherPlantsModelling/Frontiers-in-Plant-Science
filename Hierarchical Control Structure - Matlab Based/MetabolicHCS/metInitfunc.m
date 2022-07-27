function u0 = metInitfunc(Par,CIVa,O2Tank,info, C4b_O2)
t=info.step;
u=[];
u0 = NaN;

for i=1:Par.Np
    lb = [0 18/100*40.88]';
    ub = [1 24/100*40.88]';
    Aeq=[1 1];
    beq=info.O2_dem(i)-C4b_O2(i);
    %f=[0;0];
    %[x,~,EXITFLAG] = linprog(f,[],[],Aeq,beq,lb,ub,options); 
    % CCM: I HAD TO REPLACE LINPROG BY FMINCON BECAUSE OF AN EXPIRING
    % LICENSE UPDATE
    
    f = @(x) 0;
    [x,~,EXITFLAG] = fmincon(f,[0,0],[],[],Aeq,beq,lb,ub)
    
    if EXITFLAG < 0
       return
    end

    S_tmp = S(i,:)-x(2)/O2Tank.TC*Par.Ts1;
    u=[u;x];
    S(i+1,:)= S_tmp;

end
u0=u;
