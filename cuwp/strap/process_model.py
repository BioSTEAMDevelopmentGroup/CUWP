# -*- coding: utf-8 -*-
"""

"""
import biosteam as bst
from .property_package import STRAP_chemicals_outline, create_property_package, create_property_package_MSW
from .process_settings import GWP as iGWP, FFC as iFFC, WU as iWU, load_STRAP_MSW_process_settings, load_process_settings
from .systems import (
    create_single_layer_batch_separation_system, 
    create_multilayer_batch_separation_system,
    create_STRAPMSW_system
)
from .tea import create_baseline_tea
from .data import price_distributions_2023 as dist
from .data.lca_characterization_factors import indicators, set_CFs
from chaospy import distributions as shape
from cuwp import strap
from biosteam.utils import CABBI_colors, GG_colors, colors
from scipy.optimize import minimize
import flexsolve as flx
import thermosteam as tmo
import numpy as np
import os

__all__ = (
    'BaselineSTRAPProcess',
    'define_solvent',
    'define_dissolution',
    'define_precipitation',
)

CFs = indicators[iGWP]

# %% Old STRAP-MSW compositional analysis
# import biosteam as bst

# class ChemicalData:
#     __slots__ = ('titer', 'conversion', 'theoretical_yield', 'chemical')
    
#     def __init__(self, ID, titer, conversion, theoretical_yield):
#         self.titer = titer
#         self.conversion = conversion
#         self.theoretical_yield = theoretical_yield
#         if ID is not None: self.chemical = bst.Chemical(ID, db='BioSTEAM')
        
#     def __getattr__(self, name):
#         return getattr(self.chemical, name)

# # Glucose, xylose, total
# glucan = bst.Chemical('glucan', db='BioSTEAM')
# xylan = bst.Chemical('xylan', db='BioSTEAM')
# glucose = bst.Chemical('glucose')
# xylose = bst.Chemical('xylose')
# glucan = 31.3 * (glucan.MW / glucose.MW)
# xylan = 5.2 * (xylan.MW / xylose.MW)
# glucose = ChemicalData('glucose', 34.4, 17.2 / 100, 54.9 / 100)
# xylose = ChemicalData('xylose', 6.8, 3.4 / 100, 65.8 / 100)
# sugar = ChemicalData(None, 41.2, 20.6 / 100, 50.6 / 100)
# total_material = sugar.titer / sugar.conversion # per L
# total_glucan = glucose.titer / glucose.theoretical_yield * glucan.MW / glucose.MW
# total_xylan = xylose.titer / xylose.theoretical_yield * xylan.MW / xylan.MW
# biogenic_content = 0.781 
# target_polymer_fraction = 0.766 # Of the plastic
# total_without_target_plastic = 1 - (1 - biogenic_content) * target_polymer_fraction
# biogenic_content_without_PEPP = biogenic_content / total_without_target_plastic
# total_biomaterial = total_material * biogenic_content_without_PEPP
# total_lignin_and_others = total_biomaterial - total_glucan - total_xylan

# %%

bst.System.strict_convergence = False
kg_per_ton = 907.18474
kg_per_MT = 1000
L_per_gal = 3.7854
ethanol_kg_per_gal = 2.98668849
ethanol_gal_per_kg = 1. / ethanol_kg_per_gal
ethanol_L_per_kg = ethanol_gal_per_kg * L_per_gal
ethanol_kg_per_L = 1. / ethanol_L_per_kg

# https://www.recyclingpyrolysisplant.com/FAQ/pyrolysis_plant/pyrolysis-oil-applications-68.html
pyrolysis_oil_density_kg_per_m3 = 820.5 # kg / m3
pyrolysis_oil_price_range_USD_per_L = np.array([0.455, 0.78]) # USD / L
pyrolysis_oil_price_range_USD_per_kg = 1000 * pyrolysis_oil_price_range_USD_per_L / pyrolysis_oil_density_kg_per_m3
# https://pubs.acs.org/doi/10.1021/acssuschemeng.9b04763

# https://resource-recycling.com/plastics/2024/05/15/recycled-plastic-prices-continue-to-climb-higher/
recycled_PP_price_USD_per_ton = 0.06 * 2.20462 * 907.185 # Recycled PP
recycled_HDPE_price_USD_per_ton = 0.3244 * 2.20462 * 907.185 # Recycled natural HDPE

# Bloomberg original source: https://www.statista.com/statistics/1171074/price-high-density-polyethylene-forecast-globally/
HDPE_price_range = (837., 1211.)

def update_chemicals_outline(plastic, solvent):
    if plastic not in STRAP_chemicals_outline:
        STRAP_chemicals_outline.extend([
            bst.ChemicalDraft(
                plastic, # Model generic plastic film as PET.
                formula='C10H8O4',
                search_db=False,
                phase='s',
                rho=1380, # kg / m3,
                Tm=523,
                Tb=623,
                Cp=1,
                default=True,
                LHV= 21285 * 192.16812,
            ),
            bst.ChemicalDraft(
                plastic + 'oligomer', # Model generic dissolved resin as hexene
                search_ID='1-Hexene',
                CAS=plastic + 'oligomer',
            ),
        ])
    if (solvent not in STRAP_chemicals_outline 
        and solvent not in solvent_mixture_names):
        STRAP_chemicals_outline.append(solvent)

solvent_mixtures = [
]

solvent_mixture_names = set(['DMSOWater'])

def define_solvent(
        name, chemicals, composition, wt=True
    ):
    solvent_mixtures.append(
        (name, chemicals, composition, wt)
    )
    solvent_mixture_names.add(name)

def default_plastic_solvent_pair(plastic, solvent):
    update_chemicals_outline(plastic, solvent)
    define_dissolution(plastic, solvent, override=False)
    define_precipitation(plastic, solvent, override=False)

def define_dissolution(
        plastic: str,
        solvent: str,
        capacity: float=0.05,
        solvent_content: float=0.5,
        T: float=130 + 273.15,
        tau: float=0.5,
        override: bool=True,
    ):
    name = f'{plastic}_{solvent}_dissolution'
    if not override and hasattr(strap.dissolution_steps, name): return 
    def f():
        return strap.dissolution_steps.DissolutionStep(
            plastic, plastic + 'oligomer', solvent, 
            tmo.Reaction(f'{plastic} -> {plastic}oligomer', plastic, X=1.0, basis='wt'), 
            capacity, solvent_content, T, tau,
        )
    setattr(strap.dissolution_steps, name, f)
    f.__name__ = name

def define_precipitation(
        plastic: str,
        solvent: str,
        solubility: float=0,
        precipitate_solvent_content: float=0.8,
        screw_press_solvent_content: float=0.4,
        T: float=308.15,
        tau: float=0.5,
        T_condensation: float=None, # Not actually used in new configuration
        override: bool=True,
    ):
    name = f'{plastic}_{solvent}_precipitation'
    if not override and hasattr(strap.precipitation_steps, name): return 
    def f():
        return strap.precipitation_steps.PrecipitationStep(
            solvent,
            plastic,
            plastic + 'oligomer',
            solubility,
            precipitate_solvent_content,
            screw_press_solvent_content,
            T,
            tau,
            T_condensation,
        )
    setattr(strap.precipitation_steps, name, f)
    f.__name__ = name


class BaselineSTRAPProcess(bst.ProcessModel):
    """
    Create a model for a solvent targeted precipitation and dissolution process.
    The dissolution and precipitation steps default to PE.
    
    Examples
    --------
    >>> from cuwp.strap import BaselineSTRAPProcess
    >>> pm = BaselineSTRAPProcess(simulate=False)
    >>> pm.system.diagram(kind='cluster', number=True)
    >>> pm.system.simulate()
    >>> assumptions, results = pm.baseline()
    >>> assumptions
    Natural gas          Price [USD/m3]                0.167
    Feedstock            Processing capacity [MT/yr]   5e+03
                         Price [USD/kg]                 0.01
    -                    IRR [%]                        0.15
    Solvent              Price [USD/kg]                 2.17
    Polymer              Mass fraction                   0.5
    Centrifuged plastic  Solvent content [%]              50
    Plastic              Feedstock distance [km]         500
    Solvent              Solvent loss [%]                0.1
    Dissolution          Temperature [K]                 368
    Precipitation        Temperature [K]                 308
    Dissolution          Solvent capacity [wt %]           3
    dtype: float64
    
    >>> results
    -  GWP [kg*CO2e/kg]   1.54
       FFC [MJ/kg]        1.87
       MSP [USD/kg]       2.95
    dtype: float64
    
    """
    @bst.scenario
    class Scenario:
        solvent: str|tuple[str, ...] = '# Solvent used to separate the target plastic'
        target_plastic: str|tuple[str, ...] = '# The polymer layer being dissolved'
        target_plastic_percent: float|tuple[float, ...] = '# Fraction in feedstock [%]'
        processing_capacity: float = 5000, '# Feedstock flow rate [MT-plastic/yr]'
        sell_leftover_plastic: bool = False, '# Whether the MSP will include all products'
        burn_leftover_plastic: bool = True, '# Produce heat and power from leftover plastic'
        facilities: bool = True, '# On-site heat and power generation'
        precipitation_temperature_format: str = 'constant', "# Use 'drop' for % temperature drop to solvent melting point. Use 'constant' to set in Kelvin."
        precipitation_configuration: str = 'integrated heat transfer', "# Must be either 'solvent mixing' or 'integrated heat transfer'."
        turbogenerator: bool = True, '# On-site electricity generation'
        percent_inks: float = None, '# Percent inks and solubles'
        cooling_tower: bool = True, '# On-site cooling tower'
        
        @property
        def multistep(self):
            return not isinstance(self.target_plastic, str)
        
        @property    
        def N_steps(self):
            return len(self.target_plastic) if self.multistep else 1

    @classmethod
    def get_scenarios(cls):
        return tuple(cls._scenarios.values())
    
    @classmethod
    def get_scenario(cls, scenario):
        return cls._scenarios[scenario] 
    
    @classmethod
    def set_scenario(cls, scenario):
        cls._scenarios[f'{scenario.target_plastic}/{scenario.solvent}'] = scenario
    
    _scenarios = [
        Scenario('THF', 'PC', 65, 250, False, False, False),
        Scenario('Toluene', 'PE', 50, 5000, False, True, True),
        Scenario(('Toluene', 'DMSOWater'), ('PE', 'EVOH'), (50, 3.22), 5000, False, True, True),
        Scenario('Xylene', 'PE', 90, 5e3, False, False, False),
    ]
    _scenarios = {
        f'{i.target_plastic}/{i.solvent}': i for i in _scenarios
    }
    _scenarios['pilot verification'] = Scenario(
        'Xylene', 'PE', 65, 3e3, False, False, False,
        cooling_tower=False,
    )
    _scenarios['film verification'] = Scenario(
        ('Toluene', 'DMSOWater'), ('PE', 'EVOH'), (68, 10), 3e3, True, False, False,
        percent_inks=0.05,
        cooling_tower=False,
    )
    
    @property
    def name(self):
        scenario = self.scenario
        if scenario.multistep:
            name = f"{'_'.join(scenario.solvent)}_{'_'.join(scenario.target_plastic)}"
            if scenario.sell_leftover_plastic:
                name += "_sell_leftover_plastic"
            elif scenario.burn_leftover_plastic: 
                name += "_burn_leftover_plastic"
            if scenario.facilities: 
                name += "_with_facilities"
        else:
            name = f"{scenario.solvent}_{scenario.target_plastic}"
            if scenario.sell_leftover_plastic:
                name += "_sell_leftover_plastic"
            elif scenario.burn_leftover_plastic: 
                name += "_burn_leftover_plastic"
            if scenario.facilities: 
                name += "_with_facilities"
        return name
    
    @classmethod
    def default_scenario(cls):
        return cls._scenarios['PE/Toluene']
        
    @classmethod
    def as_scenario(cls, scenario):
        if isinstance(scenario, (str, tuple)):
            return cls._scenarios[scenario]
        else:
            raise TypeError('invalid scenario type')
    
    def create_thermo(self):
        scenario = self.scenario
        if scenario.multistep:
            for i, j in zip(scenario.target_plastic, scenario.solvent):
                default_plastic_solvent_pair(i, j)
        else:
            default_plastic_solvent_pair(scenario.target_plastic, scenario.solvent)
        chemicals = create_property_package()
        for i in solvent_mixtures: chemicals.define_group(*i)
        return chemicals
    
    def create_system(self):
        scenario = self.scenario
        chemicals = self.chemicals
        load_process_settings()
        if scenario.multistep:
            dissolution_step = tuple([
                getattr(strap.dissolution_steps, f'{i}_{j}_dissolution')()
                for i, j in zip(scenario.target_plastic, scenario.solvent)
            ])
            precipitation_step = tuple([
                getattr(strap.precipitation_steps, f'{i}_{j}_precipitation')()
                for i, j in zip(scenario.target_plastic, scenario.solvent)
            ])
        else:
            dissolution_step = getattr(strap.dissolution_steps, f'{scenario.target_plastic }_{scenario.solvent}_dissolution')()
            precipitation_step = getattr(strap.precipitation_steps, f'{scenario.target_plastic }_{scenario.solvent}_precipitation')()
        facilities = scenario.facilities
        target_plastic_percent = scenario.target_plastic_percent
        if scenario.multistep:
            bulk_plastic_percent = 100 - sum(target_plastic_percent) # PET
            chemicals.define_group('Plastic', [*scenario.target_plastic, 'BulkPlastic'], [*target_plastic_percent, bulk_plastic_percent], wt=True)
        else:
            bulk_plastic_percent = 100 - target_plastic_percent # PET
            chemicals.define_group('Plastic', [scenario.target_plastic, 'BulkPlastic'], [target_plastic_percent, bulk_plastic_percent], wt=True)
        chemicals.define_group('Solutes', ['Minerals', 'Solubles'], [0.8, 0.2], wt=True)
        if scenario.multistep:
            system = create_multilayer_batch_separation_system(
                dissolution_steps=dissolution_step, 
                precipitation_steps=precipitation_step,
                shred=True,
                facilities=scenario.facilities,
                core_facilities=True,
                turbogenerator=scenario.turbogenerator,
                precipitation_configuration=scenario.precipitation_configuration,
                cooling_tower=scenario.cooling_tower,
            )
            for i, step in zip(system.subsystems, dissolution_step):
                i.ID = f"{step.plastic}/{step.solvent}"
        else:
            system = create_single_layer_batch_separation_system(
                dissolution_step=dissolution_step,
                precipitation_step=precipitation_step,
                facilities=facilities,
                relative_molar_tolerance=1e-6, 
                molar_tolerance=1e-2,
                method='fixed-point',
                burn_leftover_plastic=scenario.burn_leftover_plastic,
                core_facilities=True,
                precipitation_configuration=scenario.precipitation_configuration,
                turbogenerator=scenario.turbogenerator,
                cooling_tower=scenario.cooling_tower,
            )
            for i, step in zip(system.subsystems, [dissolution_step]):
                i.ID = f"{step.plastic}/{step.solvent}"
        if scenario.multistep:
            self.dissolution_steps = dissolution_step
            self.precipitation_steps = precipitation_step
        else:
            self.dissolution_step = dissolution_step
            self.precipitation_step = precipitation_step
        self.tea = create_baseline_tea(system)
        self.direct_nonbiogenic_emissions = lambda: self.emissions.imass['CO2'] * system.operating_hours
        system.define_process_impact(
            key=iGWP,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=self.direct_nonbiogenic_emissions,
            CF=1.,
        )
        system.set_tolerance(mol=1e-6, rmol=1e-9, T=1e-6, rT=1e-9, subsystems=True, maxiter=200)      
        return system
        
    def create_model(self):
        scenario = self.scenario
        if not scenario.turbogenerator: self.BT = self.B
        facilities = scenario.facilities
        processing_capacity = scenario.processing_capacity
        target_plastic_percent = scenario.target_plastic_percent
        system = self.system
        if scenario.sell_leftover_plastic:
            products = system.outs
        else:
            _, *products = system.outs
        
        self.products = products
        model = bst.Model(system)
        parameter = model.parameter
        indicator = model.indicator
        self.tea_parameters = []
        def tea_param(f):
            self.tea_parameters.append(f)
            return f
        
        self.lca_parameters = []
        def lca_param(f):
            self.lca_parameters.append(f)
            return f
        
        self.general_parameters = []
        def gen_param(f):
            self.tea_parameters.append(f)
            self.lca_parameters.append(f)
            self.general_parameters.append(f)
            return f
        
        #global warming potential
        @indicator(units='kg*CO2e/kg')
        def GWP():
            if scenario.burn_leftover_plastic:
                GWP_material = system.get_total_feeds_impact(iGWP)
                GWP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                GWP_emissions = system.get_process_impact(iGWP) # kg CO2 eq. / y
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production # kg CO2 eq. / y
                return GWP_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                GWP = system.get_property_allocated_impact(
                    key=iGWP, name='mass', basis='kg',
                    products=products
                ) # kg-CO2e / kg
                if GWP < 0: breakpoint()
            return GWP
        
        #fossil fuel consumption
        @indicator(units='MJ/kg')
        def FFC():
            if scenario.burn_leftover_plastic:
                FFC_material = system.get_total_feeds_impact('FFC')
                FFC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                FFC_total = FFC_material  - FFC_electricity_production # kg CO2 eq. / y
                return FFC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                FFC = system.get_property_allocated_impact(
                    key='FFC', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if FFC < 0: breakpoint()
            return FFC
        
        #water usage
        @indicator(units='m3/kg')
        def WU():
            if scenario.burn_leftover_plastic:
                WU_material = system.get_total_feeds_impact('WU')
                WU_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
                #FFC_emissions = system.get_process_impact('FFC') # kg CO2 eq. / y
                WU_total = WU_material  - WU_electricity_production # kg CO2 eq. / y
                return WU_total / (self.PE_resin.F_mass * self.tea.operating_hours)
            else:
                WU = system.get_property_allocated_impact(
                    key='WU', name='mass', basis='kg',
                    products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
                ) # MJ / kg
                if WU < 0: breakpoint()
            return WU
        
        # Unverified indicators
        # #human toxicity - CANCER
        # @indicator(units='CTUh/kg')
        # def HTC():
        #     if scenario.burn_leftover_plastic:
        #         HTC_material = system.get_total_feeds_impact('HTC')
        #         HTC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         HTC_total = HTC_material  - HTC_electricity_production # kg CO2 eq. / y
        #         return HTC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         HTC = system.get_property_allocated_impact(
        #             key='HTC', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if HTC < 0: breakpoint()
        #     return HTC
        
        # #human toxicity - non CANCER
        # @indicator(units='CTUh/kg')
        # def HTNC():
        #     if scenario.burn_leftover_plastic:
        #         HTNC_material = system.get_total_feeds_impact('HTNC')
        #         HTNC_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         HTNC_total = HTNC_material  - HTNC_electricity_production # kg CO2 eq. / y
        #         return HTNC_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         HTNC = system.get_property_allocated_impact(
        #             key='HTNC', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if HTNC < 0: breakpoint()
        #     return HTNC
        
        # #acidification
        # @indicator(units='MOL H+ eq/kg')
        # def ACD():
        #     if scenario.burn_leftover_plastic:
        #         ACD_material = system.get_total_feeds_impact('ACD')
        #         ACD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         ACD_total = ACD_material  - ACD_electricity_production # kg CO2 eq. / y
        #         return ACD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         ACD = system.get_property_allocated_impact(
        #             key='ACD', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if ACD < 0: breakpoint()
        #     return ACD
        
        # #ecotoxicity
        # @indicator(units='CTU eq/kg')
        # def ETOX():
        #     if scenario.burn_leftover_plastic:
        #         ETOX_material = system.get_total_feeds_impact('ETOX')
        #         ETOX_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         ETOX_total = ETOX_material  - ETOX_electricity_production # kg CO2 eq. / y
        #         return ETOX_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         ETOX = system.get_property_allocated_impact(
        #             key='ETOX', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if ETOX < 0: breakpoint()
        #     return ETOX
        
        # #ozone depletion
        # @indicator(units='kg CFC11 eq/kg')
        # def OZD():
        #     if scenario.burn_leftover_plastic:
        #         OZD_material = system.get_total_feeds_impact('OZD')
        #         OZD_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         OZD_total = OZD_material  - OZD_electricity_production # kg CO2 eq. / y
        #         return OZD_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         OZD = system.get_property_allocated_impact(
        #             key='OZD', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if OZD < 0: breakpoint()
        #     return OZD
        
        # #photochemical ozone creation potential
        # @indicator(units='kg CFC11 eq/kg')
        # def POCP():
        #     if scenario.burn_leftover_plastic:
        #         POCP_material = system.get_total_feeds_impact('POCP')
        #         POCP_electricity_production = CFs['Electricity'] * (self.system.get_electricity_production() - self.system.get_electricity_consumption())
        #         POCP_total = POCP_material  - POCP_electricity_production # kg CO2 eq. / y
        #         return POCP_total / (self.PE_resin.F_mass * self.tea.operating_hours)
        #     else:
        #         POCP = system.get_property_allocated_impact(
        #             key='POCP', name='mass', basis='kg',
        #             products=[i for i in products if 'resin' in i.ID or 'product' in i.ID]
        #         ) # MJ / kg
        #         if POCP < 0: breakpoint()
        #     return POCP
        
        @indicator(units='USD/kg')
        def MSP():
            return self.tea.solve_price(products)
        
        V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]
        
        # https://www.eia.gov/energyexplained/natural-gas/prices.php
        # @parameter(analysis='Sobol', group='tea')
        if facilities:
            @tea_param
            @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                       baseline=4.73 * 35.3146667/1e3)
            def set_natural_gas_price(price): 
                self.BT.natural_gas_price = price * V_ng
    
        
        # Processing capacity is entirely arbitrary for now
        @parameter(
            bounds=(processing_capacity * 0.5, processing_capacity * 2),
            element='Feedstock',
            units='MT/yr',
            baseline=processing_capacity,
        )
        def set_processing_capacity(processing_capacity):
            self.feedstock.F_mass = processing_capacity * 1000 / self.tea.operating_hours
        
        # Feedstock price will be equal to transportation cost.
        # TODO: Base estimate to transportation cost on availability of 
        # post-industrial plastic per area.
        
        @tea_param
        @parameter(
            baseline=0.01,
            element='Feedstock',
            units='USD/kg',
            distribution=shape.Uniform(0, 0.02)
        )
        def set_feedstock_price(price):
            self.feedstock.price = price
        
        @gen_param
        @parameter(
            baseline = 500, #km
            element = 'feedstock',
            units='km',
            distribution=shape.Uniform(20, 2000)
        )
        def set_feedstock_distance(distance):
            GWP = 'GWP'
            FFC = 'FFC'
            WU = 'WU'
            HTC = 'HTC'
            HTNC = 'HTNC'
            ETOX = 'ETOX'
            ACD = 'ACD'
            OZD = 'OZD'
            POCP = 'POCP'
            #from Articulated lorry transport, Total weight 12-14 t, mix Euro 0-5, consumption mix, to consumer, diesel driven, Euro 0 - 5 mix, cargo, 12 - 14t gross weight / 9.3t payload capacity - RNA
            self.plastic.set_CF(GWP, distance * 0.083 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(FFC, distance * 1.08828 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(WU, distance * 0.00799 * 1/1000, ) #distance (km) * .08 kg co2 / (t*km) * t/1000kg EF 2.0
            self.plastic.set_CF(HTC, distance * 1.25436e-9 * 1/1000, )
            self.plastic.set_CF(HTNC, distance * 2.884e-9 * 1/1000, )
            self.plastic.set_CF(ETOX, distance * 0.02214 * 1/1000, )
            self.plastic.set_CF(ACD, distance * 0.00065 * 1/1000, )
            self.plastic.set_CF(OZD, distance * 2.70507e-13 * 1/1000, )
            self.plastic.set_CF(POCP, distance * 0.0006 * 1/1000, )
        
        @tea_param
        @parameter(
            element='Cashflow analysis',
            baseline=0.10,
            units='%',
            distribution=shape.Uniform(0.1, 0.2)
        )
        def set_IRR(IRR):
            self.tea.IRR = IRR
        
        if scenario.multistep:
            baseline = 100 * self.dissolution_steps[0].solvent_content
            @gen_param
            @parameter(
                distribution=shape.Uniform(baseline - 10, baseline + 10),
                element='centrifuged plastic', units='%',
            )
            def set_centrifuged_plastic_solvent_content(solvent_content):
                for i in self.dissolution_steps:
                    i.solvent_content = solvent_content / 100
            
            @gen_param
            @parameter(
                baseline=0.005,
                element='Solvent',
                units='%',
                distribution=shape.Uniform(0.001 * 100, 0.01 * 100)
            )
            def set_solvent_loss(solvent_loss):
                for i in self.dissolution_steps:
                    getattr(self, i.solvent + '_loss').split[:] = solvent_loss / 100
            
            first_step, *other_steps = self.dissolution_steps
            self.target_plastics_ratio = target_plastics_ratio = {i.plastic: 0 for i in other_steps}    
            
            def create_parameters(dissolution_step, precipitation_step):
                
                def step_parameter(*args, element=None, **kwargs):
                    return lambda f: _step_parameter(f, *args, element=element, **kwargs)
                    
                def _step_parameter(f, *args, element, **kwargs):
                    element = f'{dissolution_step.plastic} step-{element}'
                    f.__name__ = f.__name__.replace('set_', f'set_{dissolution_step.plastic}_')
                    return parameter(f, *args, element=element, **kwargs)
                
                baseline = 2.17
                @tea_param
                @step_parameter(
                    baseline=baseline,
                    element='solvent',
                    units='USD/kg',
                    distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline),
                )
                def set_solvent_price(price):
                    getattr(self, dissolution_step.solvent).price = price
                
                if dissolution_step is not first_step:
                    @gen_param
                    @step_parameter(
                        baseline=target_plastic_percent[self.dissolution_steps.index(dissolution_step)] / target_plastic_percent[0],
                        element='polymer',
                        distribution=shape.Uniform(0.1, 1),
                    )
                    def set_polymer_ratio(ratio):
                        self.target_plastics_ratio[dissolution_step.plastic] = ratio
                
                @gen_param
                @step_parameter(
                    element='dissolution', units='wt %',
                    distribution=shape.Uniform(1, 5),
                    baseline=dissolution_step.capacity * 100,
                )
                def set_dissolution_capacity(solvent_capacity):
                    dissolution_step.capacity = solvent_capacity / 100
                
                chemicals = bst.settings.chemicals
                solvent = chemicals[dissolution_step.solvent]
                if isinstance(solvent, list):
                    T = dissolution_step.T
                    Tmax = min([i.Tb for i in solvent]) - 5
                    Tmin = max(
                        max([i.Tm for i in solvent]) + 5, 265
                    )
                    if T > Tmax: T = Tmax - 1
                else:
                    T = dissolution_step.T
                    Tmax = solvent.Tb - 5
                    Tmin = max(solvent.Tm + 5, 265)
                    if T > Tmax: T = Tmax - 1
                
                @step_parameter(
                    element='dissolution', units='K', 
                    distribution=shape.Triangle(Tmin, T, Tmax)
                )
                def set_dissolution_temperature(temperature):
                    dissolution_step.T = temperature
                
                if scenario.precipitation_temperature_format == 'drop':
                    @gen_param
                    @step_parameter(
                        element='precipitation', units='%',
                        distribution=shape.Uniform(50, 100),
                    )
                    def set_precipitation_temperature_drop(temperature_drop):
                        T = dissolution_step.T
                        precipitation_step.T = (
                            T - temperature_drop / 100 * (T - Tmin)
                        )
                elif scenario.precipitation_temperature_format == 'constant':
                    T_precipitation = precipitation_step.T
                    @gen_param
                    @step_parameter(
                        element='precipitation', units='K',
                        baseline=T_precipitation,
                        distribution=shape.Uniform(T_precipitation - 5, T_precipitation + 5),
                    )
                    def set_precipitation_temperature(temperature):
                        precipitation_step.T = temperature
                
            
            for i in range(scenario.N_steps):
                create_parameters(
                    dissolution_step=self.dissolution_steps[i],
                    precipitation_step=self.precipitation_steps[i]
                )
                
            @gen_param
            @parameter(
                baseline=sum(target_plastic_percent) / 100.,
                element='target polymer',
                distribution=shape.Uniform(0.3, 0.9)
            )
            def set_polymer_mass_fraction(mass_fraction):
                self.target_polymer_mass_fraction = mass_fraction
        else:
            @gen_param
            @parameter(
                baseline=target_plastic_percent / 100.,
                element='target polymer',
                distribution=shape.Uniform(0.3, 0.9)
            )
            def set_polymer_mass_fraction(mass_fraction):
                s = self.feedstock
                F_mass = s.F_mass
                plastic = self.dissolution_step.plastic
                s.imass[plastic] = 0
                other_composition = s.mass / s.F_mass
                s.mass = other_composition * F_mass * (1 - mass_fraction)
                s.imass[plastic] = mass_fraction * F_mass
            
            def solvent_content(baseline, *args, **kwargs):
                bounds = (baseline - 10, baseline + 10)
                return parameter(*args, bounds=bounds, baseline=baseline,
                                 distribution=shape.Uniform(*bounds), 
                                 units='%', **kwargs)
            
            @gen_param
            @solvent_content(
                100 * self.dissolution_step.solvent_content,
                element='centrifuged plastic'
            )
            def set_centrifuged_plastic_solvent_content(solvent_content):
                self.dissolution_step.solvent_content = solvent_content / 100
            
            baseline = 2.17
            @tea_param
            @parameter(
                baseline=baseline,
                element='Solvent',
                units='USD/kg',
                distribution=shape.Uniform(0.5 * baseline, 1.5 * baseline)
            )
            def set_solvent_price(price):
                self.solvent.price = price
            
            @gen_param
            @parameter(
                baseline=0.001 * 100,
                element='solvent',
                units='%',
                distribution=shape.Uniform(0.0001 * 100, 0.002 * 100)
            )
            def set_solvent_loss(solvent_loss):
                self.solvent_loss.split[:] = solvent_loss / 100.
            
            @gen_param
            @parameter(
                element='dissolution', units='wt %',
                distribution=shape.Uniform(1, 5),
                baseline=self.dissolution_step.capacity * 100,
            )
            def set_dissolution_capacity(solvent_capacity):
                self.dissolution_step.capacity = solvent_capacity / 100
            
            chemicals = bst.settings.chemicals
            solvent = chemicals[self.dissolution_step.solvent]
            T = self.dissolution_step.T
            Tmax = solvent.Tb - 5
            Tmin = max(solvent.Tm + 5, 265)
            if T > Tmax: T = Tmax - 1
            @parameter(
                element='dissolution', units='K',
                distribution=shape.Triangle(Tmin, T, Tmax),
                baseline=T,
            )
            def set_dissolution_temperature(temperature):
                self.dissolution_step.T = temperature
            
            if scenario.precipitation_temperature_format == 'drop':
                @gen_param
                @parameter(
                    element='precipitation', units='%',
                    distribution=shape.Uniform(50, 100),
                )
                def set_precipitation_temperature_drop(temperature_drop):
                    T = self.dissolution_step.T
                    self.precipitation_step.T = (
                        T - temperature_drop / 100 * (T - Tmin)
                    )
            elif scenario.precipitation_temperature_format == 'constant':
                T_precipitation = self.precipitation_step.T
                @gen_param
                @parameter(
                    element='precipitation', units='K',
                    baseline=T_precipitation,
                    distribution=shape.Uniform(T_precipitation - 5, T_precipitation + 5),
                )
                def set_precipitation_temperature(temperature):
                    self.precipitation_step.T = temperature
        
        if scenario.percent_inks:
            inks = scenario.percent_inks
            @parameter(
                bounds=(0.80 * inks, 1.2 * inks),
                baseline=inks,
                units='wt % plastic',
            )
            def set_solute_content(solute_content):
                solute_content = solute_content / 100
                plastics = self.feedstock.imass['Plastic']
                self.feedstock.imass['Solutes'] = 0
                self.feedstock.imass['Plastic'] = plastics * (1 - solute_content)
                self.feedstock.imass['Solutes'] = plastics * solute_content
        
        # @gen_param
        # @solvent_content(
        #     precipitation_step.centrifuge_solvent_content,
        #     element='centrifuged precipitate'
        # )
        # def set_centrifuged_precipitate_solvent_content(solvent_content):
        #     precipitation_step.centrifuge_solvent_content = solvent_content
        
        # @gen_param
        # @solvent_content(
        #     precipitation_step.screw_press_solvent_content,
        #     element='screw pressed precipitate'
        # )
        # def set_screw_press_solvent_content(solvent_content):
        #     precipitation_step.screw_press_solvent_content = solvent_content
        
        # chemicals = bst.settings.chemicals
        # solvent = chemicals[dissolution_step.solvent]
        # solvent.Psat.method = 'BOILING_CRITICAL'
        # # This line resets the extrapolation coefficients
        # solvent.Psat.extrapolation_coeffs.clear()
        
        # @gen_param
        # @parameter(
        #     baseline=solvent.Tb,
        #     element='Solvent', units='K',
        #     distribution=shape.Uniform(solvent.Tb - 25, solvent.Tb + 25),
        # )
        # def set_boiling_point(normal_boiling_point):
        #     solvent.Tb = normal_boiling_point
        #     # This line resets the extrapolation coefficients
        #     solvent.Psat.extrapolation_coeffs.clear()
        if scenario.multistep:
            @system.add_specification(simulate=True)
            def adjust_composition():
                s = self.feedstock
                F_mass = s.F_mass
                IDs = [i.plastic for i in self.dissolution_steps]
                s.imass[IDs] = 0
                other_composition = s.mass / s.F_mass
                s.mass = other_composition * F_mass * (1 - self.target_polymer_mass_fraction)
                composition = np.array([1, *self.target_plastics_ratio.values()], dtype=float)
                composition /= composition.sum()
                s.imass[IDs] = self.target_polymer_mass_fraction * F_mass * composition
            
        self.load_model(model)
        for i in ('emissions', 'natural_gas', 'makeup_water', 'cooling_tower_makeup_water'):
            if not hasattr(self, i): setattr(self, i, bst.Stream(i))
        if facilities:
            self.natural_gas.set_CF(
                iGWP,
                0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
            )
            self.natural_gas.set_CF(
                'FFC',
                51, # [MJ / kg NG] From Open-LCA Environmental Footprint 2.0
            )
            
            #TODO: ADD NG WATER USE
        # TODO: Adjust solvent CF accordingly
        if scenario.multistep:
            for i in self.dissolution_steps:
                solvent = getattr(self, i.solvent)
                solvent.set_CF(
                    iGWP,
                    0.8199, # GREET; Mixed xylenes production from catalytic reforming of naphtha
                )
                solvent.set_CF(
                    'FFC',
                    54, # GREET; Mixed xylenes production from catalytic reforming of naphtha
                )
        else:
            self.solvent.set_CF(
                iGWP,
                0.8199, # GREET; Mixed xylenes production from catalytic reforming of naphtha
            )
            self.solvent.set_CF(
                'FFC',
                54, # GREET; Mixed xylenes production from catalytic reforming of naphtha
            )
        return model
