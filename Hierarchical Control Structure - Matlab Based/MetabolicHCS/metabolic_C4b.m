function [info, C4b_O2_Pred_ave, ref] = metabolic_C4b(X0_C4b, u_C4b, C4b_matrices, plant_ind, PARA_C4b, Cint, Oint, delta, tmod, info, model_lettuce, reactionFormulas, reactionNames, SubSystems, ref)
inputs_C4b
met_O2_C4b_Day   = zeros(size(u_C4b(in_groups_C4b),1));
met_O2_C4b_Night = zeros(size(u_C4b(in_groups_C4b),1));

for i = 1:u_C4b(in_groups_C4b)
    Tleaf = X0_C4b((i-1)*9+3);
    Pb = u_C4b(in_pressure_C4b);
    I0 = u_C4b(in_I_C4b);
    Scanopy = C4b_matrices.Scanopy(plant_ind,i);
    LAI = Scanopy*(u_C4b(in_groups_C4b)/PARA_C4b.area);
    CO2 = X0_C4b(PARA_C4b.i_CO2_C4b);
    [met_O2_C4b_Day(i), met_O2_C4b_Night(i)] = metabolic_BL(Tleaf, Pb, I0, PARA_C4b, i, LAI, u_C4b, Scanopy, CO2, Cint, Oint, delta, model_lettuce, reactionFormulas, reactionNames, SubSystems);
end

%if tmod < 16
    info.met_O2_C4b_Day(plant_ind)   = sum(met_O2_C4b_Day)*3600*32;
    info.met_O2_C4b_Night(plant_ind)  = sum(met_O2_C4b_Night)*3600*32;
    ref = plant_ind;
    C4b_O2_Pred_ave = zeros(1,6);
    for t = 1:6
        t_real = rem(t+tmod,24);
        if t_real < 16
            C4b_O2_Pred_ave(1,t) = info.met_O2_C4b_Day(plant_ind); 
        else
            C4b_O2_Pred_ave(1,t) = info.met_O2_C4b_Day(plant_ind);      % Eliminar quan s'incorpori la nit i descommentar la de més a sota
            %C4b_O2_Pred_ave(1,t) = info.met_O2_C4b_Night(plant_ind); 
        end
    end
% else
%     info.met_O2_C4b_Day(plant_ind)   = info.met_O2_C4b_Day(ref);
%     info.met_O2_C4b_Night(plant_ind) = info.met_O2_C4b_Night(ref);
%     C4b_O2_Pred_ave = zeros(1,6);
%     for t = 1:6
%         t_real = rem(t+tmod,24);
%         if t_real < 16
%             C4b_O2_Pred_ave(1,t) = info.met_O2_C4b_Day(plant_ind); 
%         else
%             C4b_O2_Pred_ave(1,t) = info.met_O2_C4b_Night(plant_ind); 
%         end
%     end
% end

end