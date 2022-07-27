%% Gross Photosynthesis in full canopies and Metabolic Modelling
% Author: Carles Ciurans, UAB, 14th May 2020
% 
% Based on papers: 
% - Johnson et al (2010) A model A model of canopy photosynthesis incorporating protein distribution through
%   the canopy and its acclimation to light, temperature and CO2
% - Thornley et al (2002) A model of canopy photosynthesis incorporating protein distribution through
% - the canopy and its acclimation to light, temperature and CO2
% - Swathy
function [Pn_Day, Pn_Night] = metabolic_CanopyChloroplastModel(Q, Tlk, P, PARA_C4b, i, LAI, Ex_CO2, Cint, Oint, model_lettuce, reactionFormulas, reactionNames, SubSystems)
R = 8.314;
Oint_Pa = Oint(1,i)*Tlk*R;
Cint_Pa = Cint(1,i)*Tlk*R;
theta           = 0.8;                % Convexity coefficient  (BETWEEN 0.7 AND 0.9)                                                          - SPECIES DEPENDANT  -                                                        [umol/(s·m2)
a               = 0.85;               % Leaf absorption coefficient for PPFD (0.85 for green leaves Nikolov (1995))
f               = 0.045;               % Energy loss factor (fraction absorbed PPFD unavailable for photsynthesis (BETWEEN 0.05 AND 0.5)       - SPECIES DEPENDANT -
fs              = 0.7;                % Direct Solar Fraction of PPFD
phi             = a*(1-f)/2;          % Efficiency of energy conversion for electron transport
Jm25            = 100;                % Rate of Jmax at 25ºC [umol m-2 s-1]                                                                   - SPECIES DEPENDANT - 
E               = 81993;              % Activation Energy of Reaction [J/mol]
R               = 8.314;              % Universal gas constant [J/(mol·K)]
S               = 711.36;             % [J mol-1 K-1]
H               = 219814;             % [J mol-1 K-1] S and H obtained from Harley (1992) for cotton leaves
Kc25            = 27e-5;              % M.M. kinetic parameters for Co2 [Pa]
Ko25            = 41e-2;              % M.M. kinetic parameters for O2 [Pa]
Vm25            = 31.31;              % Rate of Vmax at 25ºC

Thau = Oint_Pa*(213.88e-6+8.995e-6*(Tlk-273-25)+1.772e-7*(Tlk-273-25)^2);
Jmax = Jm25*exp((3.3621e-3*Tlk-1)*E/(R*Tlk))/(1+exp((S*Tlk-H)/(R*Tlk)));
Fun_GrossPhotosynthesisDirect = @(l)GrossPhotosynthesisDirect(l);
Fun_GrossPhotosynthesisDiffuse = @(l)GrossPhotosynthesisDiffuse(l);

Pgs = integral(Fun_GrossPhotosynthesisDirect, 0, LAI);  % Gross Photosynthesis Leaf Level associated to direct irradiation (mol CO2/m2 leaf/s)
Pgd = integral(Fun_GrossPhotosynthesisDiffuse, 0, LAI); % Gross Photosynthesis Leaf Level associated to diffue irradiation (mol CO2/m2 leaf/s)

    function Js = GrossPhotosynthesisDirect(l)
        Pm = Jmax*exp(-PARA_C4b.abs.*l)*10^(-6);
        Qls = PARA_C4b.abs*Q*(fs+(1-fs)*exp(-PARA_C4b.abs.*l));
        J = ((Pm+phi.*Qls)-sqrt((Pm+phi.*Qls).^2-4*theta*phi.*Qls.*Pm))/(2*theta);
        Js = J.*exp(-PARA_C4b.abs.*l);
    end
    function Jd = GrossPhotosynthesisDiffuse(l)
        Pm = Jmax*exp(-PARA_C4b.abs.*l)*10^(-6);
        Qld = PARA_C4b.abs*Q*(1-fs)*exp(-PARA_C4b.abs.*l);
        J = ((Pm+phi.*Qld)-sqrt((Pm+phi.*Qld).^2-4*theta*phi.*Qld.*Pm))/(2*theta);
        Jd = J.*(1-exp(-PARA_C4b.abs.*l));
    end
Pg = Pgs + Pgd;
Wj = Pg*10^6;

Kc = P*Kc25*exp(32.462-80470/(R*Tlk));
Ko = P*Ko25*exp(5.854-14510/(R*Tlk));
Vcmax = Vm25*exp(46.9411-116300/(R*Tlk))/(1+exp((650*Tlk-202900)/(R*Tlk)));
Wc = Vcmax*(Cint_Pa)/(Cint_Pa+Kc*(1+Oint_Pa/Ko));


shiftDay = 1;
model_lettuce = changeRxnBounds(model_lettuce,'Ex_CO2',-Ex_CO2*10^6,'l');
model_lettuce = changeRxnBounds(model_lettuce,'RuBisCo',Wc*shiftDay,'u');
model_lettuce = changeRxnBounds(model_lettuce,'RuBisCo_n',0,'b');
model_lettuce = changeRxnBounds(model_lettuce,'LETC',Wj/2*shiftDay,'u');
model_lettuce = changeRxnBounds(model_lettuce,'Ex_e',-2000*shiftDay,'l');
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
        % Reaction Ratios
model_lettuce = addRatioReaction(model_lettuce, {'RuBisCo', 'RuBisO'}, [1,3]);
model_lettuce = addRatioReaction(model_lettuce, {'RuBisCo_n', 'RuBisO_n'}, [1,1]);
model_lettuce = addRatioReaction(model_lettuce, {'GOXper', 'PGP'}, [1,1]);
model_lettuce = addRatioReaction(model_lettuce, {'Ex_CO2', 'Ex_O2'}, [-1.08, 1]);
model_lettuce = addRatioReaction(model_lettuce, {'Ex_CO2', 'Ex_CO2_night'}, [-0.27/(1-0.27), 1]);
ATPaseList = {'ATPase_cl', 'ATPase_cyt'};
ATPaseListVal = ones(size(ATPaseList,2),1)';
% Maintenance Costs: All Producers
NADPHList = {'NADPH[d]', 'NADPH[d][c]', 'NADPH[d][m]', 'NADPH[d][p]'};
[NADPHrxnList, NADPHrxnFormulaList] = findRxnsFromMets(model_lettuce, NADPHList,'producersOnly',1);
% Maintenance Costs:
NADPHmaintenanceList = {'OPPPc', 'OPPP', 'ICDHym', 'ICDHyc', 'ME', 'MEc'};
NADPHmaintenanceListVal = -ones(size(NADPHmaintenanceList,2),1)';
ATPmaintenanceList = {'Gly1_cyt', 'Gly3_cyt', 'Gly7_cyt', 'Gly10_cyt', 'Gly1', 'Gly3', 'CC3', 'Gly10'};
ATPmaintenanceListVal = [-1, -1, -1, 1, -1, -1, -1, 1];
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
model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList NADPHoxidaseList], [0.59*ATPaseListVal, NADPHoxidaseListVal], 0, 'E');
model_lettuce = addCouplingConstraint(model_lettuce, [NADPHmaintenanceList  NADPHoxidaseList], [NADPHmaintenanceListVal, NADPHoxidaseListVal*-1], 0, 'E');
ATPaseList_night = {'ATPase_cl_night', 'ATPase_cyt_night'};
ATPaseListVal_night = ones(size(ATPaseList_night,2),1)';
NADPHList_night = {'NADPH[n]', 'NADPH[n][c]', 'NADPH[n][m]', 'NADPH[n][p]'};
[NADPHrxnList_night, NADPHrxnFormulaList_night] = findRxnsFromMets(model_lettuce, NADPHList_night,'producersOnly',1);
NADPHmaintenanceList_night = {'OPPPc_n', 'OPPP_n', 'ICDHym_n', 'ICDHyc_n', 'ME_n', 'MEc_n'};
NADPHmaintenanceListVal_night = -ones(size(NADPHmaintenanceList_night,2),1)';
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
model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList_night NADPHoxidaseList_night], [0.59*ATPaseListVal_night, NADPHoxidaseListVal_night], 0, 'E');
model_lettuce = addCouplingConstraint(model_lettuce, [NADPHmaintenanceList_night  NADPHoxidaseList_night], [NADPHmaintenanceListVal_night, NADPHoxidaseListVal_night*-1], 0, 'E');
model_lettuce = addCouplingConstraint(model_lettuce, [ATPaseList  ATPaseList_night], [ATPaseListVal, ATPaseListVal_night*-1], 0, 'E');
model_lettuce = addRatioReaction(model_lettuce, {'Ex_HNO3','Ex_HNO3_night'}, [2; 3]);
[~,~] = changeCobraSolver('glpk',[],'0'); 
Model_Solved = optimizeCbModel(model_lettuce,'max');
SolMatrix = [reactionNames' reactionFormulas' SubSystems' num2cell(Model_Solved.x,2)];


Pn_Day = Model_Solved.full(13);
Pn_Night = Model_Solved.full(319);

Pn_Day = Pn_Day/10^6;
Pn_Night = Pn_Night/10^6;

end