# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 16:41:15 2022

@author: 1243617
"""


import sys
import libsbml
from libsbml import *
import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
import numpy as np
from mlxtend.plotting import heatmap
import matplotlib as plt
from pylab import *
## Using Cobra
prova = cobra.io.sbml._get_doc_from_filename('Model2.xml')
lettuce = cobra.io.sbml._sbml_to_model(prova)
FBA1 = lettuce.optimize()
lettuce.objective = "Biomass"
lettuce.solver = 'glpk'
csfont = {'fontname':'Times New Roman'}

sensitivity_points = 8
sensitivity_fractions = np.array([(1-0.2), (1-(0.2-0.2/3*1)), (1-(0.2-0.2/3*2)), (1-(0.2-0.2/3*3)), (1+0.2-(0.2/3*2)), (1+0.2-(0.2/3*1)), (1+0.2)])
central_values = np.array([70.45/2, 1.22, 3, 3, 0.25])
LETC_fractions = np.array([42.99, 53.55, 60.25, 70.45, 78.52, 84, 87.44])/2
span_LETC = np.array([[(53.55-42.99)/53.55, (60.25-53.55)/60.25, (70.45-60.25)/70.45, (78.52-70.45)/78.52, (84-78.52)/84, (87.44-84)/87.44]])
span_general = np.full((len(central_values)-1, len(sensitivity_fractions)-1),0.1)
span = np.concatenate((span_LETC, span_general), axis=0)

cobrasummary_columns = list(range(len(sensitivity_fractions)*len(central_values)+2))
cobrasummary_index   = list(range(len(lettuce.reactions)))
cobrasummary = pd.DataFrame(columns = cobrasummary_columns, index = cobrasummary_index)
variationsummary_columns = list(range(len(central_values)+2))
variationsummary_index   = list(range(len(lettuce.reactions)))
variationsummary = pd.DataFrame(columns = variationsummary_columns, index = variationsummary_index)


contador = -1
for j in range(len(central_values)):
    
    for i in range(len(sensitivity_fractions)):
        
        prova = cobra.io.sbml._get_doc_from_filename('Model2.xml')
        lettuce = cobra.io.sbml._sbml_to_model(prova)
        FBA1 = lettuce.optimize()
        lettuce.objective = "Biomass"
        lettuce.solver = 'glpk'
        
        
        ## Define Upper and Lower Bounds
        Ex_CO2 = -59;
        Wc = 19.75;
        #Wc = 30;
        Wo = 0.85;
        lettuce.reactions.RuBisCo.upper_bound = Wc
        lettuce.reactions.RuBisO.lower_bound = Wo

        ## Define Ratios
        RubiscoCarboxylation_night = lettuce.problem.Constraint(
            lettuce.reactions.RuBisCo_n.flux_expression,
            lb = 0,
            ub = 0)
        RubiscoOxygenation_night = lettuce.problem.Constraint(
            lettuce.reactions.RuBisO_n.flux_expression,
            lb = 0,
            ub = 0)
        Photorespiration_day = lettuce.problem.Constraint(
            lettuce.reactions.SGTper.flux_expression - lettuce.reactions.GTper.flux_expression,
            lb = 0,
            ub = 0)
        Photorespiration_night = lettuce.problem.Constraint(
            lettuce.reactions.SGTper_n.flux_expression - lettuce.reactions.GTper_n.flux_expression,
            lb = 0,
            ub = 0) 
        Nitrate_Export = lettuce.problem.Constraint(
            2*lettuce.reactions.Ex_HNO3.flux_expression - 3*lettuce.reactions.Ex_HNO3_night.flux_expression,
            lb = 0,
            ub = 0)
        NightResp = lettuce.problem.Constraint(
             -1*lettuce.reactions.Ex_CO2_night.flux_expression - lettuce.reactions.Ex_O2_night.flux_expression,
             lb = 0,
             ub = 0)
        NADPHMaintenance_day = lettuce.problem.Constraint(
            lettuce.reactions.OPPP.flux_expression+lettuce.reactions.OPPPc.flux_expression+lettuce.reactions.ICDHym.flux_expression+lettuce.reactions.ICDHyc.flux_expression+lettuce.reactions.ME.flux_expression+lettuce.reactions.MEc.flux_expression
            - 1*(lettuce.reactions.NADPHoxidase_cl.flux_expression+lettuce.reactions.NADPHoxidase_mit.flux_expression+lettuce.reactions.NADPHoxidase_cyt.flux_expression),
            lb = 0,
            ub = 0)
        NADPHMaintenance_night = lettuce.problem.Constraint(
            lettuce.reactions.OPPPc_n.flux_expression+lettuce.reactions.OPPP_n.flux_expression+lettuce.reactions.ICDHym_n.flux_expression+lettuce.reactions.ICDHyc_n.flux_expression+lettuce.reactions.ME_n.flux_expression+lettuce.reactions.MEc_n.flux_expression
            - 1*(lettuce.reactions.NADPHoxidase_cl_night.flux_expression+lettuce.reactions.NADPHoxidase_mit_night.flux_expression+lettuce.reactions.NADPHoxidase_cyt_night.flux_expression),
            lb = 0,
            ub = 0)

        ## Introduce Ratios
        lettuce.add_cons_vars(RubiscoCarboxylation_night)
        lettuce.add_cons_vars(RubiscoOxygenation_night)
        lettuce.add_cons_vars(Photorespiration_day)
        lettuce.add_cons_vars(Photorespiration_night)
        lettuce.add_cons_vars(Nitrate_Export)
        lettuce.add_cons_vars(NightResp)
        
        
        
        contador = contador+1
        central_values = np.array([70.45/2, 1.22, 3, 3, 0.25])
        tempo = central_values
        if j == 0:
            tempo[j] = LETC_fractions[i]
        else:
            tempo[j] = central_values[j]*sensitivity_fractions[i]
        
        lettuce.reactions.LETC.upper_bound = tempo[0]
        PhotosynthesisRatio = lettuce.problem.Constraint(
            -tempo[1]*lettuce.reactions.Ex_CO2.flux_expression - lettuce.reactions.Ex_O2.flux_expression,
            lb = 0,
            ub = 0)
        ATPaseNADPHoxidase_day = lettuce.problem.Constraint(
            (tempo[2]*lettuce.reactions.ATPase_cl.flux_expression+lettuce.reactions.ATPase_mit.flux_expression+lettuce.reactions.ATPase_cyt.flux_expression)
            - (lettuce.reactions.NADPHoxidase_cl.flux_expression+lettuce.reactions.NADPHoxidase_mit.flux_expression+lettuce.reactions.NADPHoxidase_cyt.flux_expression),
            lb = 0,
            ub = 0)
        ATPaseNADPHoxidase_night = lettuce.problem.Constraint(
            (tempo[3]*lettuce.reactions.ATPase_cl_night.flux_expression+lettuce.reactions.ATPase_mit_night.flux_expression+lettuce.reactions.ATPase_cyt_night.flux_expression)
            - (lettuce.reactions.NADPHoxidase_cl_night.flux_expression+lettuce.reactions.NADPHoxidase_mit_night.flux_expression+lettuce.reactions.NADPHoxidase_cyt_night.flux_expression),
            lb = 0,
            ub = 0)
        # NADPHMaintenance_day = lettuce.problem.Constraint(
        #     lettuce.reactions.OPPP.flux_expression+lettuce.reactions.OPPPc.flux_expression+lettuce.reactions.ICDHym.flux_expression+lettuce.reactions.ICDHyc.flux_expression+lettuce.reactions.ME.flux_expression+lettuce.reactions.MEc.flux_expression
        #     - tempo[4]*(lettuce.reactions.NADPHoxidase_cl.flux_expression+lettuce.reactions.NADPHoxidase_mit.flux_expression+lettuce.reactions.NADPHoxidase_cyt.flux_expression),
        #     lb = 0,
        #     ub = 0)
        # NADPHMaintenance_night = lettuce.problem.Constraint(
        #     lettuce.reactions.OPPPc_n.flux_expression+lettuce.reactions.OPPP_n.flux_expression+lettuce.reactions.ICDHym_n.flux_expression+lettuce.reactions.ICDHyc_n.flux_expression+lettuce.reactions.ME_n.flux_expression+lettuce.reactions.MEc_n.flux_expression
        #     - tempo[5]*(lettuce.reactions.NADPHoxidase_cl_night.flux_expression+lettuce.reactions.NADPHoxidase_mit_night.flux_expression+lettuce.reactions.NADPHoxidase_cyt_night.flux_expression),
        #     lb = 0,
        #     ub = 0)
        # MaintenanceDayNight = lettuce.problem.Constraint(
        #     lettuce.reactions.ATPase_cl.flux_expression+lettuce.reactions.ATPase_mit.flux_expression+lettuce.reactions.ATPase_cyt.flux_expression
        #     - tempo[4]*(lettuce.reactions.ATPase_cl_night.flux_expression+lettuce.reactions.ATPase_mit_night.flux_expression+lettuce.reactions.ATPase_cyt_night.flux_expression),
        #     lb = 0,
        #     ub = 0)
        DayNightResp = lettuce.problem.Constraint(
            -tempo[4]*lettuce.reactions.Ex_CO2.flux_expression - lettuce.reactions.Ex_CO2_night.flux_expression,
            lb = 0,
            ub = 0)

        lettuce.add_cons_vars(DayNightResp)
        lettuce.add_cons_vars(PhotosynthesisRatio)
        lettuce.add_cons_vars(ATPaseNADPHoxidase_day)
        lettuce.add_cons_vars(ATPaseNADPHoxidase_night)
        lettuce.add_cons_vars(NADPHMaintenance_day)
        lettuce.add_cons_vars(NADPHMaintenance_night)
        #lettuce.add_cons_vars(MaintenanceDayNight)
        
        # Solve FBA
        FBA = lettuce.optimize()

        #lettuce.reactions.get_by_id("Exp_sucrosePhloemDay_night").lower_bound = 0
        for k in range(len(lettuce.reactions)):
            ValFBA = FBA[k]
            
            name = lettuce.reactions[k].name
            expression = lettuce.reactions[k].reaction
            
            cobrasummary.loc[k][0] = name
            cobrasummary.loc[k][1] = expression
            cobrasummary.loc[k][contador+2] = ValFBA
            
variationsummary_columns = list(range(len(central_values)+2))
variationsummary_index   = list(range(len(lettuce.reactions)))
variationsummary = pd.DataFrame(columns = variationsummary_columns, index = variationsummary_index)  
#Populate Name and Formula Reaction         
for k in range(len(lettuce.reactions)):
    name = lettuce.reactions[k].name
    expression = lettuce.reactions[k].reaction
    variationsummary.loc[k][0] = name
    variationsummary.loc[k][1] = expression 
    

# Calculation of variation as: abs(SUM(slope))(max-min)
calc = np.array([])
size = len(sensitivity_fractions)
calc = np.empty(len(sensitivity_fractions)-1, dtype=object)

for k in range(len(lettuce.reactions)):
    counter2 = -1
    for m in range(len(central_values)):
        counter2 = counter2+1
        counter3 = -1
        for n in range(len(sensitivity_fractions)-1):
            counter3 = counter3+1
            calc[n] = (abs((cobrasummary.loc[k][size*counter2+counter3+3]-cobrasummary.loc[k][size*counter2+counter3+2]))/span[m,n])/(max(cobrasummary.loc[k][size*counter2+2:size*(counter2+1)+1])-min(cobrasummary.loc[k][size*counter2+2:size*(counter2+1)+1]))
            calc[np.isnan(calc[n])]=0
            if max(cobrasummary.loc[k][size*counter2+2:size*(counter2+1)+1])-min(cobrasummary.loc[k][size*counter2+2:size*(counter2+1)+1]) < 1e-3:
                calc[n] = 0
            
        #calc = np.nan_to_num(calc)
        #sum_slope = np.sum(calc)
        variationsummary.loc[k][m+2] = np.sum(calc)
# -min(cobrasummary.loc[k][size*counter2+2:size*(counter2+1)+1]) 
        
        
# Generation of Heat Map
RatioNames = ["$I_{u}$", "$Ex_{O_{2}}[d]:Ex_{CO_{2}}[d]$", "$ATP_{maint}[d]:NADPH_{maint}[d]$", "$ATP_{maint}[n]:NADPH_{maint}[n]$", "$Ex_{CO_{2}}[d]:Ex_{CO_{2}}[n]$"]
np_variationsummary =  variationsummary.to_numpy()
Day_Names = ["MP1", "MP2", "PRK", "MDHm", "G3P_trans"]
daynames = np.empty(len(Day_Names), dtype=object)
day_variation = np.empty([len(Day_Names), len(central_values)+2], dtype=object)

for p in range(len(Day_Names)):
    daynames[p] = np.argwhere((np_variationsummary == Day_Names[p]))
    daynames[p] = daynames[[p]][0][0][0]
    day_variation[p,:] = np_variationsummary[daynames[p]]
    
Night_Names = ["MP1_n", "MP2_n", "GAPDc_n","MDHm_n", "Citrate_DayImport"]
nightnames = np.empty(len(Night_Names), dtype=object)
night_variation = np.empty([len(Night_Names), len(central_values)+2], dtype=object)
for p in range(len(Night_Names)):
    nightnames[p] = np.argwhere((np_variationsummary == Night_Names[p]))
    nightnames[p] = nightnames[[p]][0][0][0] 
    night_variation[p] = np_variationsummary[nightnames[p]]
    
minval_day = np.amin(day_variation[:,2:7][day_variation[:,2:7] != -np.inf])
maxval_day = norm_night_variation = night_variation[:,2:8].astype(int)/maxval_night
np.amax(day_variation[:,2:7][day_variation[:,2:7] != np.inf])

minval_night = np.amin(night_variation[:,2:7][night_variation[:,2:7] != -np.inf])
maxval_night = np.amax(night_variation[:,2:7][night_variation[:,2:7] != np.inf])

# Divide both matrix by maximum of day and night
norm_day_variation = day_variation[:,2:8].astype(int)/maxval_day

fig, ax = heatmap(norm_day_variation, column_names=RatioNames, row_names=Day_Names, figsize=(20, 10), cell_values=False)
plt.rcParams.update({'font.family':'Times'})
plt.title('a) Normalized Flux Variability of Day Metabolism', fontsize=16, **csfont)
plt.figure(dpi=1200)
for tick in ax.get_xticklabels():
    tick.set_fontname("Times New Roman")
for tick in ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()  
fig.tight_layout()
fig.savefig('DayMet.jpg', dpi=1200)


#plt.yaxis.set_label_position("right")

fig, ax = heatmap(norm_night_variation, column_names=RatioNames, row_names=Night_Names, figsize=(20, 10), cell_values=False)
plt.rcParams.update({'font.family':'Times'})
plt.title('b) Normalized Flux Variability of Night Metabolism', fontsize=16, **csfont)
for tick in ax.get_xticklabels():
    tick.set_fontname("Times New Roman")
for tick in ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
fig.tight_layout()
fig.savefig('NightMet.jpg', dpi=1200)
#ax.yaxis.set_label_position("right")
#ax.yaxis.tick_right()  


