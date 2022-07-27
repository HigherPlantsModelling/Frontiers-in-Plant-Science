%% Description: MAIN FILE C4b. Call of list of ODEs
% Author: Carles Ciurans, UAB, 2020
% Based on: UCA gas-transfer models
% Version: Adapted to incorporate Metabolic-HCS
addpath('..\')
%% Load Initialized Data
load X0_C3_init_v5;
load X0_C5_init_v5;
load u_C5C4b_init_v5;
load u_C3C4b_init_v5;
load u_C4b_init_v5;
load X0_HPC_init_v5;
load X0_C4b_init_v5;
load X0_Res_init_v5;
load u_C3_init_v5;
load model_lettuce;
Lettuce_Diel_Model_30thJuly_Categories
%% Initialize Parameters
inputs_C4b;             % index of Inputs
load_compounds_HPC
load_compounds_Res
u_C4b(in_groups_C4b) = 4;
Par = systemparameters;
PARA_C4b = parameters_C4b(u_C4b, Par);
U = input_matrix(PARA_C4b);
U_MPC = zeros(Par.Np,1);
Y_MPC = zeros(Par.Np+1,1);
ref = 1;
C4b_matrices = C4b_matrices_func(PARA_C4b,u_C4b, Par);
C4bfunction;
u_C4b = C4bInputs(X0_C3, X0_C5, u_C5C4b, u_C3C4b, PARA_C4b.Air_Con);
[C4b_CO2Controller, C4b_HNO3Controller, C4b_O2Controller, C4b_O2ControllerConc] = localcontroller_C4b(PARA_C4b, Par);
num_groups = u_C4b(in_groups_C4b);
Cint = 1.012e-2*ones(2,u_C4b(in_groups_C4b));
Oint = 8.4*ones(2,u_C4b(in_groups_C4b));
%% Initialize Matrix
Par.Nsim = 24*30;                      % Simulation interval [step]
C4b_CO2_pred_saved          = zeros(Par.Nsim*(1/Par.Ts4),1);
C4b_O2_pred_saved           = zeros(Par.Nsim*(1/Par.Ts4),2);
C4b_CO2_pred_saved(1,1)     = info.X0_C4b(1,PARA_C4b.i_CO2_C4b)/100/PARA_C4b.Air_Con;
C4b_O2_pred_saved(1,1)      = info.X0_C4b(1,PARA_C4b.i_O2_C4b)/100/PARA_C4b.Air_Con;
C4b_O2_pred_saved(1,2)      = info.X0_C4b(1,PARA_C4b.i_O2_C4b)/100/PARA_C4b.Air_Con;
Y_overall                   = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,PARA_C4b.i_HNO3_C4b);
Scanopy_overall             = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,u_C4b(in_groups_C4b));
Carbon_Assimilation_overall = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,u_C4b(in_groups_C4b));
info.met_O2_C4b_Day         = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,1);
info.met_O2_C4b_Night       = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,1);
Us                          = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,1);
Ys                          = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,1);
SSTOflag                    = zeros(PARA_C4b.n*PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time,1);
SSTOflag(1,1)               = 1;
Fresh_medium = u_C4b(in_Mc_C4b:end)'; 
delta = ones(u_C4b(in_groups_C4b),1)';

%% Initialize States    
Y01 = X0_C4b(1:PARA_C4b.n_compounds_C4b);
Yin = X0_C4b(end-3:end);    
y_C4b(1,:) = X0_C4b; 
Y0_C4b(1,:) = X0_C4b;
Y(1,:) = Y0_C4b(1,:);
%% Initialize Storing Matrix
for i = 1:u_C4b(in_groups_C4b)
    C4b_matrices.Scanopy(i) = PARA_C4b.K1*X0_C4b(PARA_C4b.i_Mc_C4b+PARA_C4b.n_compounds_C4b*(i-1)) /PARA_C4b.DM;   
    C4b_matrices.L(i)=2*sqrt(C4b_matrices.Scanopy(i)/pi);                                                                               % Leaf Characteristic Length m (aka Diameter)
    C4b_matrices.Lstem(i)=PARA_C4b.K2/PARA_C4b.DM*X0_C4b(PARA_C4b.i_Mc_C4b+PARA_C4b.n_compounds_C4b*(i-1));                                           % Stem Length m
    C4b_matrices.Nvessel(i)=PARA_C4b.K3/PARA_C4b.DM*X0_C4b(PARA_C4b.i_Mc_C4b+PARA_C4b.n_compounds_C4b*(i-1));
    C4b_matrices.LAI0(i) = (PARA_C4b.plants/PARA_C4b.area)*C4b_matrices.Scanopy(i)/(PARA_C4b.plants/u_C4b(in_groups_C4b));                                                                                % Initialization Leaf Area Index
    C4b_matrices.RL0(i) = (C4b_matrices.L(i)+C4b_matrices.Lstem(i))/PARA_C4b.Ratio;                                                          % Initialization Root Length considering Shoot is L + Lstem
    C4b_matrices.RA0(i) = 2*pi*C4b_matrices.RL0(i)*PARA_C4b.r0;                                                                           % Initialization Root Area
    C4b_matrices.Lim0(i) =0;
end  
C4b_matrices.time(1) = 0;
C4b_matrices.RL = C4b_matrices.RL0';
C4b_matrices.LAI = C4b_matrices.LAI0';
C4b_matrices.RA = C4b_matrices.RA0';
RL = C4b_matrices.RL;
LAI = C4b_matrices.LAI;
RA = C4b_matrices.RA;
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
C4b_matrices.Collection = [];
C4b_matrices.Crop_Array(:,1) = 0:Par.Ts1:Par.Nsim/Par.Ts1-Par.Ts1';
C4b_matrices.Crop_Array_Mng(:,1) = 0:Par.Ts1:PARA_C4b.batch_shift*24*(u_C4b(in_groups_C4b));  
C4b_matrices.SChamber(:,1) = 0:Par.Ts1:Par.Nsim/Par.Ts1-Par.Ts1';
C4b_matrices.C4b_matrices.Crop_Cell = cell(1,u_C4b(in_groups_C4b)+(PARA_C4b.n-1));                                                     % 4 is ther basal number of crops, when n = 2 we will have 5 different crops                 
C4b_matrices.C4b_matrices.Crop_Cell_Array_Mng = cell(1,u_C4b(in_groups_C4b));
for k=1:size(C4b_matrices.C4b_matrices.Crop_Cell,2)
            C4b_matrices.C4b_matrices.Crop_Cell{k} = C4b_matrices.Crop_Array;             
end    
for k=1:size(C4b_matrices.C4b_matrices.Crop_Cell_Array_Mng,2)
            C4b_matrices.C4b_matrices.Crop_Cell_Array_Mng{k} = C4b_matrices.Crop_Array_Mng;             
end    
   
%% Initialize indexes
PARA_C4b.index = 1;
plant_ind = 0;
sec_cont_ind = 1;

        figControl = figure;
        p1 = plot(1,1,'g');
        ylabel('Flow (mL/min');
        hold on
        yyaxis right
        p2 = plot(1,1,'-r');
        p3 = plot(1,1,'r'); 
        p4 = plot(1,1,'b');
        hold off
%        ylim([0. 0.2])
        ylabel('Concentration(%)')
        legend('ControlCommand','Ref O_2','Measurement O_2', 'Pred MPC O_2','Location','northoutside');
        xlabel('Time(h)')
%        
%        figO2 = figure;
%        p4 = plot(1,1,'g');
%        ylabel('m3/s');
%        hold on
%        yyaxis right
%        p5 = plot(1,1,'-r');
%        p6 = plot(1,1,'r'); 
%        yyaxis left
%        p7 = plot(1,1,'o');
%        hold off
%        ylabel('Concentration(%)')
%        legend('Gas Flow','Setpoint O2','Measurement O2', 'O2 Conc');
%        xlabel('Time(d)')
%        
    for i = 1:PARA_C4b.n 
       
            for k = 1:PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time
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
                O2ref        = U(PARA_C4b.index,4);
                % Controller 2: SSTO-MPC
                [info, Pn, ref]           = metabolic_C4b(X0_C4b, u_C4b, C4b_matrices, plant_ind, PARA_C4b, Cint, Oint, delta, tmod, info, model_lettuce, reactionFormulas, reactionNames, SubSystems, ref);
                [Us, Ys, SSTOflag]        = SSTO(U, Pn, u_C4b, plant_ind, Us, Ys, X0_C4b, Y_MPC, PARA_C4b, SSTOflag);
                [U_MPC, Y_MPC, flow]      = MPC(U, Pn, u_C4b, plant_ind, U_MPC, Y_MPC, PARA_C4b, X0_C4b, Us, Ys, Par);
                if plant_ind == 1
                    error(plant_ind) = 1;
                    flow = 0;
                    U_MPC(:,1) = 0;
                else
                    error(plant_ind) =  X0_C4b(PARA_C4b.i_O2_C4b)-Y_MPC(2,plant_ind-1);
                end
                % Update Control Commands
                u_C4b(in_I_C4b)        = U(PARA_C4b.index,1);                     % umol/m2
                u_C4b(in_temp_C4b)     = U(PARA_C4b.index,2);                     % K
                u_C4b(in_RH_C4b)       = U(PARA_C4b.index,3);                     % 1/1
                u_C4b(in_gas_flow_C4b) = flow;                                    % 1/1  
                u_C4b(in_gas_co2_C4b)  = 0.4/100/PARA_C4b.Air_Con;                           % mol O2/m3
                
                % Update Night Respiration
                Night_CO2 = C4b_matrices.Carbon_Assimilation(PARA_C4b.index-1,:);
                if u_C4b(in_I_C4b) == 0
                    Night_CO2 = C4b_matrices.Carbon_Assimilation(PARA_C4b.contador_day,:)*0.34;
                else
                    PARA_C4b.contador_day = PARA_C4b.index;
                end 
                
                % Perturbation at time 8h
                if k == (10)*3600/(PARA_C4b.sample_time)
                    X0_C4b(1,PARA_C4b.i_O2_C4b) = 20.5/100/PARA_C4b.Air_Con;
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
                set(p1,'XData',C4b_matrices.time(2:plant_ind+1)/3600,'YData',U_MPC(1,1:plant_ind)*1000/60);
                set(p2,'XData',C4b_matrices.time(2:plant_ind+1)/3600,'YData',Ys(1:plant_ind)*100*PARA_C4b.Air_Con);
                set(p3,'XData',C4b_matrices.time(2:plant_ind+1)/3600,'YData',Y(2:plant_ind+1,37)*100*PARA_C4b.Air_Con);
                set(p4,'XData',C4b_matrices.time(2:plant_ind+1)/3600,'YData',Y_MPC(2,1:plant_ind)*100*PARA_C4b.Air_Con);
%                 set(p5,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',PARA_C4b.SPO2*ones(size(C4b_matrices.time(1:PARA_C4b.index-1),1),1)*100*PARA_C4b.Air_Con);
%                 set(p6,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',Y(1:PARA_C4b.index-1,PARA_C4b.i_O2_C4b)*100*PARA_C4b.Air_Con);
%                 set(p7,'XData',C4b_matrices.time(1:PARA_C4b.index-1)/3600/24,'YData',MV(1:PARA_C4b.index-1,3));
                drawnow;
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
    
         
        