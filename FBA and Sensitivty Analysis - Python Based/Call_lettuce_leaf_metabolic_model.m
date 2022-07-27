%% Lettuce Leaf Metabolic Model
% Author: Carles Ciurans, UAB
%         Igor Batolomé, Ivan Martinez, DTU
% Description: Call lettuce_leaf_m  etabolic_model with the following inputs:
%   - Carbon Fixation by RUBISCO  (Wc)
%   - Light Reaction (Wj)
% FILES: 
%   - Model Name: Lettuce Leaf Model
%           Model Inputs: hour (day/night shift) from 0-16 day, from 16-24
%                               night
% Aim: Obtain a matrix of flows driven by optimizing production of biomass

% MODELING SCENARIOS defined as "modelling_scenarios"
%   case 1: Scan the range of possible Wc/Wj with all exchangeable flows
%           from 0 - 1000 by FBA
%   case 2: FVA analysis. Fixing Wc at its minimum and obtain Wj that makes
%           Ratio Ex_CO2/Ex_O2 = -1.2
%   case 3: Fix Rubsico and light rates and plot Flux Balance Analysis
%           Results
close all

% Step 1: Definition of parameters for posterior calculation or used
% directly in constraint definitions
    Q = 400; % W/m2
    h700 = Q*4.7;
    Resp = 0.27;
    Ex_CO2 = -59;
    Wc = 19.75;
    Wo = 0.85;
    Wj = 42;
    shiftDay = 1;
% Step 2: Load Metabolic Model (don't pay attention on the 1 as input) and
% store it in variable model_lettuce
        [model_lettuce, reactionNames, reactionFormulas, SubSystems] = LettuceMetabolicModel(1);
% Step 3: Define upper and lower boundaries of the optimization problem
        model_lettuce = changeRxnBounds(model_lettuce,'Ex_CO2',Ex_CO2,'l');
        model_lettuce = changeRxnBounds(model_lettuce,'RuBisCo',Wc,'u');
        model_lettuce = changeRxnBounds(model_lettuce,'RuBisO',Wo,'u');
        model_lettuce = changeRxnBounds(model_lettuce,'RuBisCo_n',0,'b');
        model_lettuce = changeRxnBounds(model_lettuce,'LETC',Wj,'u');
        model_lettuce = changeRxnBounds(model_lettuce,'Ex_e',-2000,'l');
        model_lettuce = changeRxnBounds(model_lettuce,'Ex_H2S',0,'l');
        model_lettuce = changeRxnBounds(model_lettuce,'ME',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'ME_n',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEyc',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEyc_n',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEc',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEc_n',0,'u');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEm',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MEm_n',0,'u');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'ENO',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'ENO_n',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'GAPDH_n',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Ser_bio_cl_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Ser_bio_cyt_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Ser_bio_cl',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Ser_bio_cyt',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Xu5P_tra_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'E4P_trans_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'G6P_trans_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'Prot10_night',0,'u');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'LP5_night',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'GAPDHc',-93,'l');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'NonP_GAPDHyc',0.33,'u');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'NADconv',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'NADconv_n',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'NADPHoxidase_cl',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'NADPHoxidase_cyt',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'ATPase_mit',0,'b');       % Dark-inhibited
        model_lettuce = changeRxnBounds(model_lettuce,'MDH',-0.75,'l');       % Dark-inhibited
% Step 4: Add Restrictions (of the type v1 = 2·v2)
        %model_lettuce = addRatioReaction(model_lettuce, {'RuBisCo', 'RuBisO'}, [1,3]);
        model_lettuce = addRatioReaction(model_lettuce, {'RuBisCo_n', 'RuBisO_n'}, [1,1]);
        model_lettuce = addRatioReaction(model_lettuce, {'GOXper', 'PGP'}, [1,1]);
        model_lettuce = addRatioReaction(model_lettuce, {'Ex_CO2', 'Ex_O2'}, [-1.08, 1]);
        model_lettuce = addRatioReaction(model_lettuce, {'Ex_CO2', 'Ex_CO2_night'}, [-0.27/(1-0.27), 1]);
        model_lettuce = addRatioReaction(model_lettuce, {'Ex_HNO3','Ex_HNO3_night'}, [2; 3]);
% Step 5: Add Advanced Restrictions (of the type v1+2·v2+v3=4·v16)
        ATPaseList = {'ATPase_cl', 'ATPase_cyt'};
        ATPaseListVal = ones(size(ATPaseList,2),1)';
        % Maintenance Costs: All Producers
        NADPHList = {'NADPH[d]', 'NADPH[d][c]', 'NADPH[d][m]', 'NADPH[d][p]'};
        [NADPHrxnList, NADPHrxnFormulaList] = findRxnsFromMets(model_lettuce, NADPHList,'producersOnly',1);
        % Maintenance Costs:
        NADPHmaintenanceList = {'OPPPc', 'OPPP', 'ICDHym', 'ICDHyc', 'ME', 'MEc'};
        NADPHmaintenanceListVal = -ones(size(NADPHmaintenanceList,2),1)';
        % ATP MAintenance Costs:  
        ATPmaintenanceList = {'Gly1_cyt', 'Gly3_cyt', 'Gly7_cyt', 'Gly10_cyt', 'Gly1', 'Gly3', 'CC3', 'Gly10'};
        ATPmaintenanceListVal = [-1, -1, -1, 1, -1, -1, -1, 1];
        %NADPH oxidase
        NADPHoxidaseList = {'NADPHoxidase_cl', 'NADPHoxidase_mit', 'NADPHoxidase_cyt'};
        NADPHoxidaseListVal = -ones(size(NADPHoxidaseList,2),1)';
        NADPHProducers = 1; % 1 - for all producers //  2 - for biomass precursors 
        switch NADPHProducers       
            case 1 
                NADPHrxnList(ismember(NADPHrxnList,NADPHmaintenanceList)) = [];
                NADPHrxnListVal = 2*ones(size(NADPHrxnList,2),1)';
            case 2
                NADPHrxnList = {'Prot3', 'Prot19', 'Prot26', 'LP7', 'LP8', 'LP9', 'LP11'};
                NADPHrxnListVal = 2*ones(size(NADPHrxnList,2),1)';
        end      
        ATPaseList_night = {'ATPase_cl_night', 'ATPase_cyt_night'};
        ATPaseListVal_night = ones(size(ATPaseList_night,2),1)';
        % Maintenance Costs: All Producers
        NADPHList_night = {'NADPH[n]', 'NADPH[n][c]', 'NADPH[n][m]', 'NADPH[n][p]'};
        [NADPHrxnList_night, NADPHrxnFormulaList_night] = findRxnsFromMets(model_lettuce, NADPHList_night,'producersOnly',1);
        % Maintenance Costs: All Maintenance
        NADPHmaintenanceList_night = {'OPPPc_n', 'OPPP_n', 'ICDHym_n', 'ICDHyc_n', 'ME_n', 'MEc_n'};
        NADPHmaintenanceListVal_night = -ones(size(NADPHmaintenanceList_night,2),1)';
        %NADPH oxidase
        NADPHoxidaseList_night = {'NADPHoxidase_cl_night', 'NADPHoxidase_mit_night', 'NADPHoxidase_cyt_night'};
        NADPHoxidaseListVal_night = -ones(size(NADPHoxidaseList_night,2),1)';
        NADPHProducers = 2; % 1 - for all producers //  2 - for biomass precursors 
        switch NADPHProducers       
            case 1 
                NADPHrxnList_night(ismember(NADPHrxnList_night,NADPHmaintenanceList_night)) = [];
                NADPHrxnListVal_night = ones(size(NADPHrxnList_night,1),1)';
            case 2
                NADPHrxnList_night = {'Prot3_night', 'Prot19_night', 'Prot26_night', 'LP7_night', 'LP8_night', 'LP9_night', 'LP11_night'};
                NADPHrxnListVal_night = 2*ones(size(NADPHrxnList_night,2),1)';
        end   
        model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList NADPHoxidaseList], [0.59*ATPaseListVal, NADPHoxidaseListVal], 0, 'E');
        model_lettuce = addCouplingConstraint(model_lettuce, [NADPHmaintenanceList  NADPHoxidaseList], [NADPHmaintenanceListVal, NADPHoxidaseListVal*-1], 0, 'E');
        model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList_night NADPHoxidaseList_night], [0.59*ATPaseListVal_night, NADPHoxidaseListVal_night], 0, 'E');
        model_lettuce = addCouplingConstraint(model_lettuce, [NADPHmaintenanceList_night  NADPHoxidaseList_night], [NADPHmaintenanceListVal_night, NADPHoxidaseListVal_night*-1], 0, 'E');
        model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList  ATPaseList_night], [ATPaseListVal, ATPaseListVal_night*-1], 0, 'E');
% Step6: Solve Optimization problem (FVA)        
        Model_Solved = optimizeCbModel(model_lettuce,'max');
% Step 7: Store solution
        SolMatrix = [reactionNames' reactionFormulas' SubSystems' num2cell(Model_Solved.x,2)];