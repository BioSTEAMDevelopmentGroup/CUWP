# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:28:33 2021

@author: yrc2
"""

__all__ = (
    'set_CFs',
    'GWP',
    'WU',
    'FFC',
    'indicators',
)

GWP = 'GWP' # Global warming potential 100-yr
WU = 'WU' # Water consumption
FFC = 'FFC' # Fossil fuel consumption
indicators = {}

# All values in cradle-to-gate except for CH4, which is in cradle-to-grave
indicators[GWP] = { # Material GWP cradle-to-gate [kg*CO2*eq / kg]
    'cellulase': 8.0482, # GREET
    'H3PO4': 1.0805, # GREET
    'lime': 1.2824, # GREET
    'MgSO4': 3.9902e-1, # Ecoinvent 3.6, IPCC 2013
    'urea': 1.0488, # GREET
    'HCl': 1.9873, # GREET
    'NaOH': 2.0521, # GREET
    'NaOCH3': 1.5871, # Ecoinvent, TRACI, sodium methoxide
    'CH4': 0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
    'Electricity': 0.3869, # [kg*CO2*eq / kWhr] From GREET; NG-Fired Simple-Cycle Gas Turbine CHP Plant
    'gasoline': 0.8415, # GREET, Gasoline blendstock USA
    'citric acid': 1.4749, # GREET
    
    # The following ecoinvent entries are from v3.8, allocation at the point of substitution
    # Ecoinvent, market for sodium hydrogen sulfite, GLO,
    # converted to 38% solution
    'bisulfite': 1.2871,

    # Ecoinvent, market for sodium hypochlorite, without water, in 15% solution state, RoW,
    # converted to 12.5 wt% solution (15 vol%)
    'NaOCl': 2.4871,
    
    # Corn dry-grind ethanol (GREET), To well to biorefinery-gate (no distribution and storage) 
    'corn-ethanol': 1.0741,
    
    # Corn stover ethanol (GREET), To well to biorefinery-gate (no distribution and storage) 
    'cornstover-ethanol': 0.2972,
}
indicators[FFC] = { # MJ / kg
    'H3PO4': 15, 
    'cellulase': 96,
    'lime': 4.868,
    'MgSO4': 1.514, # Update
    'urea': 27, 
    'HCl': 29,
    'NaOH': 29,
    'NaOCH3': 22, # Update
    'CH4': 54,
    'Electricity': 6.303,
    'gasoline': 53,
    'citric acid': 16, # GREET
    'bisulfite': 18, # Update
    'NaOCl': 35, # Update
    'corn-ethanol': 10,
    'cornstover-ethanol': 4.533,
}
indicators[WU] = { # kg / kg
    'H3PO4': 36.7572505,
    'cellulase': 197.7,
    'lime': 4.6528090,
    'MgSO4': 1.4477260,
    'urea': 4.6678101,
    'HCl': 5.1341402,
    'NaOH': 13.6079225,
    'NaOCH3': 10.5244061, # Update
    'CH4': 0.5513673,
    'Electricity': 0.7014350,
    'gasoline': 4.2677178,
    'citric acid': 26.9968670, 
    'bisulfite': 8.5350407, # Update
    'NaOCl': 16.4925023, # Update
    'corn-ethanol': 30.6762116,
    'cornstover-ethanol': 10.7020759, # 6.127890632318501 comes from corn stover
}


default_dilutions = {
    'NaOCl': 0.125,
    'bisulfite': 0.38,
}

def set_CFs(stream, name, dilution=None):
    if dilution is None:
        dilution = default_dilutions.get(name, 1)
    for i in indicators:
        if i == 'WU':
            stream.characterization_factors[i] = indicators[i][name] * dilution + (1 - dilution)
        else:
            stream.characterization_factors[i] = indicators[i][name] * dilution
    
    
    
# from thermosteam.units_of_measure import convert
# from thermosteam import Chemical
# CH4 = Chemical('CH4')
# CO2 = Chemical('CO2')
# electricty_produced_per_kg_CH4 = - convert(0.8 * 0.85 * CH4.LHV / CH4.MW, 'kJ', 'kWhr')
# GWP_per_kg_CH4 = 0.33 + CO2.MW / CH4.MW
# GWP_per_kWhr = GWP_per_kg_CH4 / electricty_produced_per_kg_CH4