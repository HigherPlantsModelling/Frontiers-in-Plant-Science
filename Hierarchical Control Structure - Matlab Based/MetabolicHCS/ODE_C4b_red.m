%% Description: MAIN FILE C4b. Call of list of ODEs
% Author: Carles Ciurans, UAB, 2020
% Based on: UCA gas-transfer models
% Version: Adapted to incorporate Metabolic-HCS
addpath('..\')
%% Load Initialized Data
load ('Workspace_2h.mat')
k0 = k;
PARA_C4b.index = PARA_C4b.index-1;
plant_ind = plant_ind-1;
sec_cont_ind = sec_cont_ind-1;
    for i = 1:PARA_C4b.n 
            for k = 233:PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time
                tmod=rem(k*PARA_C4b.sample_time/3600,24);
                % Update indexes and object variables
                PARA_C4b.index = PARA_C4b.index+1;
                plant_ind = plant_ind+1;
                sec_cont_ind = sec_cont_ind + 1;  
                concentration_class = hObj([]);
                delta_class = hObj([]);
                Valor_class = hObj([]);
                Carbon_Assimilation_class = hObj([]);
                Scanopy_global_class = hObj([]);
                Exitflag_class = hObj([]);
                Pgs_class = hObj([]);
                Pgd_class = hObj([]);
                Wc_class = hObj([]);
                Wj_class = hObj([]);
                concentration_class.o = [Cint Oint];
                delta_class.o = delta;
                
                % Controller 1: CO2
                [C4b_CO2Controller,C4b_CO2_pred_saved, CO2Inj] = LocalController_C4b_CO2(C4b_CO2Controller, C4b_CO2_pred_saved, sec_cont_ind, X0_HPC, PARA_C4b);
                MV(PARA_C4b.index,1) = CO2Inj;      % Output in seconds
                      
                % Controller 2: SSTO-MPC
                [info, Pn, ref] = metabolic_C4b(X0_C4b, u_C4b, C4b_matrices, plant_ind, PARA_C4b, Cint, Oint, delta, tmod, info, model_lettuce, reactionFormulas, reactionNames, SubSystems, ref);
                [Us, Ys]        = SSTO(21, Pn, u_C4b, plant_ind, Us, Ys, X0_C4b, Y_MPC, PARA_C4b);
                [U_MPC, Y_MPC, flow]   = MPC(21, Pn, u_C4b, plant_ind, U_MPC, Y_MPC, PARA_C4b, X0_C4b, Us, Ys, Par);
                if plant_ind == 1
                    error(plant_ind) = 1;
                else
                    error(plant_ind) =  X0_C4b(PARA_C4b.i_O2_C4b)-Y_MPC(2,plant_ind-1);
                end
                u_C4b(in_I_C4b)        = U(PARA_C4b.index,1);                     % umol/m2
                u_C4b(in_temp_C4b)     = U(PARA_C4b.index,2);                     % K
                u_C4b(in_RH_C4b)       = U(PARA_C4b.index,3);                     % 1/1
                u_C4b(in_gas_flow_C4b) = flow;                                    % 1/1  
                u_C4b(in_gas_co2_C4b)  = 0.4/100*40.88;                           % mol O2/m3
                
                Night_CO2 = C4b_matrices.Carbon_Assimilation(PARA_C4b.index-1,:);
                
                if u_C4b(in_I_C4b) == 0
                    Night_CO2 = C4b_matrices.Carbon_Assimilation(PARA_C4b.contador_day,:)*0.34;
                else
                    PARA_C4b.contador_day = PARA_C4b.index;
                end   
                    
                options = odeset('RelTol',1e-6, 'AbsTol', 1e-6);
                tic
                [T,y]=ode45(@(t,y) dyn_C4b_test(t, y, u_C4b, PARA_C4b, CO2Inj, concentration_class, delta_class, nicolet, u_C3, Night_CO2, Valor_class, Carbon_Assimilation_class, ...
                Scanopy_global_class, Exitflag_class, Pgs_class, Pgd_class, Wc_class, Wj_class),[0,PARA_C4b.sample_time],X0_C4b', options); 
                Cint = concentration_class.o(:,1:u_C4b(in_groups_C4b));
                Oint = concentration_class.o(:,u_C4b(in_groups_C4b)+1:2*u_C4b(in_groups_C4b));
                delta = delta_class.o(1,:);                
                toc
                
                C4b_matrices.time(PARA_C4b.index,1) = T(end)+C4b_matrices.time(PARA_C4b.index-1);
                time = C4b_matrices.time(PARA_C4b.index,1);
                Y(PARA_C4b.index,:) = y(end,:);
                X0_C4b = Y(PARA_C4b.index,:);
                X0_HPC(i_CO2_HPC) = X0_C4b(PARA_C4b.i_CO2_C4b);
                X0_HPC(i_O2_HPC) = X0_C4b(PARA_C4b.i_O2_C4b);
                Y_overall(plant_ind,:) = Y(PARA_C4b.index,:);
                C4b_matrices.Limitation(PARA_C4b.index,:) = Valor_class.o(1,:);
                C4b_matrices.Scanopy(PARA_C4b.index,:) = Scanopy_global_class.o(1,:);
                C4b_matrices.Carbon_Assimilation(PARA_C4b.index,:) = Carbon_Assimilation_class.o(1,:); 
                C4b_matrices.delta(PARA_C4b.index,:) = delta(1,:);
                C4b_matrices.exitflag(PARA_C4b.index,:) = Exitflag_class.o(1,:);
                C4b_matrices.Cint(PARA_C4b.index,:) = Cint(1,:);
                C4b_matrices.Wc(PARA_C4b.index,:) = Wc_class.o(1,:);
                C4b_matrices.Wj(PARA_C4b.index,:) = Wj_class.o(1,:);
                C4b_matrices.Pgs(PARA_C4b.index,:) = Pgs_class.o(1,:);
                C4b_matrices.Pgd(PARA_C4b.index,:) = Pgd_class.o(1,:);
                C4b_matrices.RootLength(PARA_C4b.index,:) = RL;
                C4b_matrices.LeafAreaIndex(PARA_C4b.index,:) = LAI;
                C4b_matrices.RootArea(PARA_C4b.index,:) = RA;
                
                Scanopy_overall(plant_ind,:) = Scanopy_global_class.o(1,:);
                Carbon_Assimilation_overall(plant_ind,:) = Carbon_Assimilation_class.o(1,:); 
                
                U_C4b(PARA_C4b.index,:) = u_C4b;
%                 set(p1,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',MV(1:PARA_C4b.index-1,1));
%                 set(p2,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',PARA_C4b.SP*ones(size(C4b_matrices.time(1:PARA_C4b.index-1),1),1)*100*PARA_C4b.Air_Con);
%                 set(p3,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',Y(1:PARA_C4b.index-1,PARA_C4b.i_CO2_C4b)*100*PARA_C4b.Air_Con);
%                 set(p4,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',MV(1:PARA_C4b.index-1,2));
%                 set(p5,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',PARA_C4b.SPO2*ones(size(C4b_matrices.time(1:PARA_C4b.index-1),1),1)*100*PARA_C4b.Air_Con);
%                 set(p6,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',Y(1:PARA_C4b.index-1,PARA_C4b.i_O2_C4b)*100*PARA_C4b.Air_Con);
%                 set(p7,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',MV(1:PARA_C4b.index-1,3));
%                 drawnow;
%                 %Choose only 1 hr data
               
            end
            harvest_time = rem(C4b_matrices.time(PARA_C4b.index,1)/3600/24,PARA_C4b.batch_shift);
            [Y0_C4b, C4b_matrices] = C4b_FilteredData(C4b_matrices, Y, PARA_C4b, Y01, u_C4b, U_C4b, Par);      
            X0_C4b = Y0_C4b;
            PARA_C4b.harvest_counter = PARA_C4b.harvest_counter + 1;
            PARA_C4b.index = 1;
            Y(1,:) = Y0_C4b(1,:);
            C4b_matrices.LeafAreaIndex(1,:) = LAI;     
            C4b_matrices.RootLength(1,:) = RL;
            C4b_matrices.RootArea(1,:) = RA;
            C4b_matrices.RootLength(1,1) = C4b_matrices.RL0(1,1);
            C4b_matrices.RootArea(1,1) = C4b_matrices.RA0(1,1);
            C4b_matrices.LeafAreaIndex(1,1) = C4b_matrices.LAI0(1,1); 
            PARA_C4b_matrix.Limitation(1,:) = C4b_matrices.Limitation(end,:);
            RL = C4b_matrices.RootLength(1,:);
            LAI = C4b_matrices.LeafAreaIndex(1,:);
            RA = C4b_matrices.RootArea(1,:);

        end
    
         
        