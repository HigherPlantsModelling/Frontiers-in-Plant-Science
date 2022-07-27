function [Pn_Day, Pn_Night] = metabolic_BL(Tbl, Pb, I0, PARA_C4b, i, LAI, u_C4b, Scanopy, CO2, Cint, Oint, delta, model_lettuce, reactionFormulas, reactionNames, SubSystems)
inputs_C4b
gbl_CO2 = PARA_C4b.DCO2*Pb/(PARA_C4b.R*Tbl*delta(1,i));
GCO2 = gbl_CO2*PARA_C4b.gs_CO2/(gbl_CO2+PARA_C4b.gs_CO2);
Ex_CO2 = PARA_C4b.DCO2 /delta(1,i) * LAI* (CO2 - Cint(1,i));

[Pn_Day, Pn_Night] = metabolic_CanopyChloroplastModel(I0, Tbl, Pb, PARA_C4b, i, LAI, Ex_CO2, Cint, Oint, model_lettuce, reactionFormulas, reactionNames, SubSystems);
Pn_Day = Pn_Day*(PARA_C4b.area/u_C4b(in_groups_C4b));
Pn_Night = Pn_Night*(PARA_C4b.area/u_C4b(in_groups_C4b));

end