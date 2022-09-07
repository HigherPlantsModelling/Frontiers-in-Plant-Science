# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 16:31:53 2022

@author: 1243617
"""

import sys
import libsbml
from libsbml import *
import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from cobra.flux_analysis.loopless import add_loopless, loopless_solution

## Using Cobra
prova = cobra.io.sbml._get_doc_from_filename('Model2.xml')
lettuce = cobra.io.sbml._sbml_to_model(prova)
FBA1 = lettuce.optimize()
lettuce.objective = "Biomass"
lettuce.solver = 'glpk'
cobrasummary_columns = ['Name', 'Reaction', 'FBA']
cobrasummary_index   = list(range(len(lettuce.reactions)))
cobrasummary = pd.DataFrame(columns = cobrasummary_columns, index = cobrasummary_index)

## Define Upper and Lower Bounds - Obtained From Matlab Hierarchical Model
# 1000 ppm CO2
Ex_CO2 = -59;
Wc = 19.75;
Wo = 0.85;
Wj = 38;

# 400 ppm CO2
Ex_CO2 = -56;
Wc = 7.82;
Wo = 1.9;
Wj = 37;


lettuce.reactions.RuBisCo.upper_bound = Wc
lettuce.reactions.RuBisO.lower_bound = Wo
lettuce.reactions.LETC.upper_bound = Wj
lettuce.reactions.GAPDHy.upper_bound = 0
lettuce.reactions.GAPDHy.lower_bound = 0
lettuce.reactions.GAPDHy_n.upper_bound = 0
lettuce.reactions.GAPDHy_n.lower_bound = 0
lettuce.reactions.MDHy_n.upper_bound = 0
lettuce.reactions.MDHy_n.lower_bound = 0
lettuce.reactions.PEPCx.upper_bound = 0
lettuce.reactions.PEPCx.lower_bound = 0
lettuce.reactions.PEPCx_n.upper_bound = 0
lettuce.reactions.PEPCx_n.lower_bound = 0
#lettuce.reactions.Ser_bio_cyt.upper_bound = 0
#lettuce.reactions.Ser_bio_cyt.lower_bound = 0
#lettuce.reactions.Ser_bio_cl.upper_bound = 0
#lettuce.reactions.Ser_bio_cl.lower_bound = 0
# lettuce.reactions.MEyc_n.upper_bound = 0
# lettuce.reactions.MEyc_n.lower_bound = 0
# lettuce.reactions.MEc_n.upper_bound = 0
# lettuce.reactions.MEc_n.lower_bound = 0
# lettuce.reactions.MEy_n.upper_bound = 0
# lettuce.reactions.MEy_n.lower_bound = 0
# lettuce.reactions.ME_n.upper_bound = 0
# lettuce.reactions.ME_n.lower_bound = 0
# lettuce.reactions.MEm_n.upper_bound = 0
# lettuce.reactions.MEm_n.lower_bound = 0

# Removing NAD reactions
#lettuce.reactions.NADconv_mit.lower_bound = 0
#lettuce.reactions.NADconv_mit.upper_bound = 0
#lettuce.reactions.NADconv_mit_n.lower_bound = 0
#lettuce.reactions.NADconv_mit_n.upper_bound = 0
#lettuce.reactions.NADconv.lower_bound = 0
#lettuce.reactions.NADconv.upper_bound = 0
#lettuce.reactions.NADconv_n.lower_bound = 0
#lettuce.reactions.NADconv.upper_bound = 0
#lettuce.reactions.NADconv_cl.lower_bound = 0
#lettuce.reactions.NADconv_cl.upper_bound = 0
#lettuce.reactions.NADconv_cl_n.lower_bound = 0
#lettuce.reactions.NADconv_cl_n.upper_bound = 0
#lettuce.reactions.Ex_HNO3.lower_bound = -0.5

## Break Loops
# lettuce.reactions.GS.upper_bound = 0.2
# lettuce.reactions.GS.lower_bound = -0.2
# lettuce.reactions.PYKc.upper_bound = 0.2
# lettuce.reactions.PYKc.lower_bound = -0.2
# lettuce.reactions.PGM.upper_bound = 0.2
# lettuce.reactions.PGM.lower_bound = -0.2
# lettuce.reactions.ENO.upper_bound = 0.2
# lettuce.reactions.ENO.lower_bound = -0.2
# lettuce.reactions.Ser_bio_cl.upper_bound = 0.2
# lettuce.reactions.Ser_bio_cl.lower_bound = -0.2
# lettuce.reactions.CysBiosynth.upper_bound = 0.2
# lettuce.reactions.CysBiosynth.lower_bound = -0.2
# lettuce.reactions.GOGAT.upper_bound = 0.2
# lettuce.reactions.GOGAT.lower_bound = -0.2
# lettuce.reactions.Prot32.upper_bound = 0.2
# lettuce.reactions.Prot32.lower_bound = -0.2
# lettuce.reactions.Prot33.upper_bound = 0.2
# lettuce.reactions.Prot33.lower_bound = -0.2

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
PhotosynthesisRatio = lettuce.problem.Constraint(
    -1.22*lettuce.reactions.Ex_CO2.flux_expression - lettuce.reactions.Ex_O2.flux_expression,
    lb = 0,
    ub = 0)
ATPaseNADPHoxidase_day = lettuce.problem.Constraint(
    (lettuce.reactions.ATPase_cl.flux_expression+lettuce.reactions.ATPase_mit.flux_expression+lettuce.reactions.ATPase_cyt.flux_expression)
    - 3*(lettuce.reactions.NADPHoxidase_cl.flux_expression+lettuce.reactions.NADPHoxidase_mit.flux_expression+lettuce.reactions.NADPHoxidase_cyt.flux_expression),
    lb = 0,
    ub = 0)
ATPaseNADPHoxidase_night = lettuce.problem.Constraint(
    (lettuce.reactions.ATPase_cl_night.flux_expression+lettuce.reactions.ATPase_mit_night.flux_expression+lettuce.reactions.ATPase_cyt_night.flux_expression)
    - 3*(lettuce.reactions.NADPHoxidase_cl_night.flux_expression+lettuce.reactions.NADPHoxidase_mit_night.flux_expression+lettuce.reactions.NADPHoxidase_cyt_night.flux_expression),
    lb = 0,
    ub = 0)
NADPHMaintenance_day = lettuce.problem.Constraint(
    lettuce.reactions.OPPP.flux_expression+lettuce.reactions.OPPPc.flux_expression+lettuce.reactions.ICDHym.flux_expression+lettuce.reactions.ICDHyc.flux_expression+lettuce.reactions.ME.flux_expression+lettuce.reactions.MEc.flux_expression
    - (lettuce.reactions.NADPHoxidase_cl.flux_expression+lettuce.reactions.NADPHoxidase_mit.flux_expression+lettuce.reactions.NADPHoxidase_cyt.flux_expression),
    lb = 0,
    ub = 0)
NADPHMaintenance_night = lettuce.problem.Constraint(
    lettuce.reactions.OPPPc_n.flux_expression+lettuce.reactions.OPPP_n.flux_expression+lettuce.reactions.ICDHym_n.flux_expression+lettuce.reactions.ICDHyc_n.flux_expression+lettuce.reactions.ME_n.flux_expression+lettuce.reactions.MEc_n.flux_expression
    - (lettuce.reactions.NADPHoxidase_cl_night.flux_expression+lettuce.reactions.NADPHoxidase_mit_night.flux_expression+lettuce.reactions.NADPHoxidase_cyt_night.flux_expression),
    lb = 0,
    ub = 0)
MaintenanceDayNight = lettuce.problem.Constraint(
    lettuce.reactions.OPPPc_n.flux_expression+lettuce.reactions.OPPP_n.flux_expression+lettuce.reactions.ICDHym_n.flux_expression+lettuce.reactions.ICDHyc_n.flux_expression+lettuce.reactions.ME_n.flux_expression+lettuce.reactions.MEc_n.flux_expression
    - (lettuce.reactions.OPPPc.flux_expression+lettuce.reactions.OPPP.flux_expression+lettuce.reactions.ICDHym.flux_expression+lettuce.reactions.ICDHyc.flux_expression+lettuce.reactions.ME.flux_expression+lettuce.reactions.MEc.flux_expression),
    lb = 0,
    ub = 0)
DayNightResp = lettuce.problem.Constraint(
     -0.25*lettuce.reactions.Ex_CO2.flux_expression - lettuce.reactions.Ex_CO2_night.flux_expression,
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

## Introduce Ratios
lettuce.add_cons_vars(RubiscoCarboxylation_night)
lettuce.add_cons_vars(RubiscoOxygenation_night)
lettuce.add_cons_vars(Photorespiration_day)
lettuce.add_cons_vars(Photorespiration_night)
lettuce.add_cons_vars(PhotosynthesisRatio)
lettuce.add_cons_vars(ATPaseNADPHoxidase_day)
lettuce.add_cons_vars(ATPaseNADPHoxidase_night)
lettuce.add_cons_vars(NADPHMaintenance_day)
lettuce.add_cons_vars(NADPHMaintenance_night)
lettuce.add_cons_vars(MaintenanceDayNight)
lettuce.add_cons_vars(Nitrate_Export)
lettuce.add_cons_vars(DayNightResp)
lettuce.add_cons_vars(NightResp)

# Loopless
loop_reactions = [lettuce.reactions.PYKc, 
                  lettuce.reactions.PGM, 
                  lettuce.reactions.ENO, 
                  lettuce.reactions.EB1, 
                  lettuce.reactions.EB2, 
                  lettuce.reactions.ACS, 
                  lettuce.reactions.Ser_bio_cl, 
                  lettuce.reactions.GOGAT,
                  lettuce.reactions.Prot32, 
                  lettuce.reactions.OASTL,
                  lettuce.reactions.GS]

if __name__ == '__main__':
    FVA_Sol_loopless = flux_variability_analysis(lettuce,reaction_list=loop_reactions,loopless=True)
    lettuce.reactions.get_by_id("PYKc").upper_bound = FVA_Sol_loopless['maximum']['PYKc']
    lettuce.reactions.get_by_id("PYKc").lower_bound = FVA_Sol_loopless['minimum']['PYKc']
    lettuce.reactions.get_by_id("PGM").upper_bound = FVA_Sol_loopless['maximum']['PGM']
    lettuce.reactions.get_by_id("PGM").lower_bound = FVA_Sol_loopless['minimum']['PGM']
    lettuce.reactions.get_by_id("ENO").upper_bound = FVA_Sol_loopless['maximum']['ENO']
    lettuce.reactions.get_by_id("ENO").lower_bound = FVA_Sol_loopless['minimum']['ENO']
    lettuce.reactions.get_by_id("EB1").upper_bound = FVA_Sol_loopless['maximum']['EB1']
    lettuce.reactions.get_by_id("EB1").lower_bound = FVA_Sol_loopless['minimum']['EB1']
    lettuce.reactions.get_by_id("EB2").upper_bound = FVA_Sol_loopless['maximum']['EB2']
    lettuce.reactions.get_by_id("EB2").lower_bound = FVA_Sol_loopless['minimum']['EB2']
    lettuce.reactions.get_by_id("ACS").upper_bound = FVA_Sol_loopless['maximum']['ACS']
    lettuce.reactions.get_by_id("ACS").lower_bound = FVA_Sol_loopless['minimum']['ACS']
    lettuce.reactions.get_by_id("Ser_bio_cl").upper_bound = FVA_Sol_loopless['maximum']['Ser_bio_cl']
    lettuce.reactions.get_by_id("Ser_bio_cl").lower_bound = FVA_Sol_loopless['minimum']['Ser_bio_cl']
    lettuce.reactions.get_by_id("GOGAT").upper_bound = FVA_Sol_loopless['maximum']['GOGAT']
    lettuce.reactions.get_by_id("GOGAT").lower_bound = FVA_Sol_loopless['minimum']['GOGAT']
    lettuce.reactions.get_by_id("Prot32").upper_bound = FVA_Sol_loopless['maximum']['Prot32']
    lettuce.reactions.get_by_id("Prot32").lower_bound = FVA_Sol_loopless['minimum']['Prot32']
    lettuce.reactions.get_by_id("OASTL").upper_bound = FVA_Sol_loopless['maximum']['OASTL']
    lettuce.reactions.get_by_id("OASTL").lower_bound = FVA_Sol_loopless['minimum']['OASTL']
    lettuce.reactions.get_by_id("GS").upper_bound = FVA_Sol_loopless['maximum']['GS']
    lettuce.reactions.get_by_id("GS").lower_bound = FVA_Sol_loopless['minimum']['GS']
# Solve FBA
    FBA = lettuce.optimize()
#lettuce.reactions.get_by_id("Exp_sucrosePhloemDay_night").lower_bound = 0
    for i in range(len(lettuce.reactions)):
        ValFBA = FBA[i]
        name = lettuce.reactions[i].name
        expression = lettuce.reactions[i].reaction
        cobrasummary.loc[i].Name = name
        cobrasummary.loc[i].Reaction = expression
        cobrasummary.loc[i].FBA = ValFBA
    
    prova = cobrasummary.loc[645].FBA/cobrasummary.loc[654].FBA

