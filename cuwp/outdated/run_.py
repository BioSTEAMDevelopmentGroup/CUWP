from biosteam import units
import json
import biosteam as bst
from biosteam import settings
import numpy as np
import cond

"""
Import file Setup.py with all common settings for the
simulation (chemicals, utilities, prices, etc.)

Additional settings must be added there instead of here
"""

# General commments:
# I tried to generalize the solvent variable, only when it is "" we need to specify the solvent
# For this script, we are using dodecane

#this function is used to run the process two times 
def run_simulation(input_data, price_solvent):
    #impact factors from openLCA 
    #electricity: Electricity grid mix 1kV-60kV , consumption mix, to consumer, AC, technology mix, 1kV - 60kV
    impact_elec = 0.152997222 #kg CO2-eq/MJ 
    #steam: steam production, in chemical industry
    impact_steam = 0.123945455 #kg CO2-eq/MJ 
    #cooling water: tap water production, artificial recharged wells
    impact_water = 0.42 #kg CO2-eq/m3


    input_data["Price"]["solvent"] = price_solvent

    chemicals = bst.Chemicals(
        ["Dodecane", "DMSO", bst.Chemical(
            "C2H4", phase="s"), "Water", "glucose", "N2", "O2", "CH4", "CO2", "C5H8O2"]

    )
    bst.settings.set_thermo(chemicals)

    bst.CE = input_data['CE']

    bst.PowerUtility.price = input_data["power_utility_price"]
    HeatUtility = bst.HeatUtility
    steam_utility1 = HeatUtility.get_agent("high_pressure_steam")
    steam_utility2 = HeatUtility.get_agent("low_pressure_steam")
    steam_utility3 = HeatUtility.get_agent("medium_pressure_steam")
    # HeatUtility.heating_agents([steam_utility1,steam_utility2])
    steam_utility1.heat_transfer_efficiency = 0.85
    steam_utility1.regeneration_price = 0.08064
    steam_utility1.T = 529.2
    steam_utility1.P = 44e5
    steam_utility2.heat_transfer_efficiency = 0.85
    steam_utility2.regeneration_price = 0.06768
    steam_utility2.T = 428.6
    steam_utility2.P = 55e4
    steam_utility3.heat_transfer_efficiency = 0.85
    steam_utility3.regeneration_price = 0.07974
    steam_utility3.T = 480.3
    steam_utility3.P = 18e5

    product_stream = bst.Stream(
        "PROD",
        price=2.4023,
    )

    # Solvent Feed tank
    V100 = units.MixTank(
        "V100", tau=input_data["tau_solv_feed_tank"], vessel_material="Carbon steel"
    )

    # new mixer for output solvent streams of dryer
    M101 = units.Mixer("M101")

    # Feeder
    FE101 = units.StorageTank(
        "FE101", tau=input_data["tau_feeder"], vessel_material="Carbon steel"
    )
    # process_units["FE101"] = FE101

    # Dissolution vessel
    V101 = units.MixTank(
        "V101", vessel_material="Carbon steel", tau=input_data["tau_dissolution_vessel"]
    )
    # process_units["V101"] = V101

    # Heater
    E101 = units.HXutility(
        "E101", T=input_data["heating_temp_after_dissolution"])
    # process_units["E101"] = E101

    # Hot filtration
    FHOT = units.RVF("FHOT", moisture_content=0,
                     split=input_data["splits_hot_filter"])
    # process_units["FHOT"] = FHOT

    # Pump 1
    P101 = units.Pump("P101")
    # process_units["P101"] = P101

    # Precipitation vessel
    V102 = units.MixTank(
        "V102",
        vessel_material="Carbon steel",
        tau=input_data["tau_precipitation_vessel"],
    )
    # process_units["V102"] = V102

    # Cooler
    E102 = units.HXutility(
        "E102", T=input_data["cooling_temp_after_precipitation"])
    # process_units["E102"] = E102

    # Pump 2
    P102 = units.Pump("P102")
    # process_units["P102"] = P102

    # Filtration #tol_44%
    F101 = units.RVF("F101", moisture_content=0,
                     split=input_data["splits_filter"])
    # process_units["F101"] = F101

    # Dryer
    H101 = units.DrumDryer(
        ID="H101",
        ins=("Mix_s6"),
        outs=("dry_P1", "solvent_out_1"),
        thermo=None,
        split=input_data["splits_dryer"],
        R=1.4,
        H=20,
        length_to_diameter=25,
        T=input_data["temp_dryer"],
        moisture_content=0,
        utility_agent="Natural gas",
    )
    # process_units["H101"] = H101

    # # Condenser
    # E103 = units.HXutility("E103", T=input_data["temp_condenser"], rigorous=True)
    # process_units["E103"] = E103

    # # Condenser
    E103 = units.HXutility("E103", T=cond.Temp, rigorous=True)
    # process_units["E103"] = E103

    # Splitter
    S101 = units.PhaseSplitter("S101", outs=("gas", "liquid_solvent"))
    # process_units["S101"] = S101
    # top gas

    @S101.add_specification(run=True)
    def phases():
        S101.ins[0].phases = "gl"

    # Pump 6
    P106 = units.Pump("P106")
    # process_units["P106"] = P106

    # Mixer
    M100 = units.Mixer("M100")
    # process_units["M100"] = M100

    # Dryer
    H1012 = units.DrumDryer(
        ID="H1012",
        ins=("Mix_s6"),
        outs=("dry_P2", "solvent_out_2"),
        thermo=None,
        split=input_data["splits_dryer_EVOH"],
        R=1.4,
        H=20,
        length_to_diameter=25,
        T=input_data["temp_dryer"],
        moisture_content=0,
        utility_agent="Natural gas",
    )

    # Mixer
    # we need this mixer for the outputs and then estimate the MSP with product_stream-"PROD"
    M1003 = units.Mixer("M1003", outs=product_stream)
    # process_units["M1003"] = M1003

    inlet_pe = input_data["total_PE"]

    inlet_solvent = input_data["total_solvent"]

    N_runs = 10

    class STRAPTEA(bst.TEA):
        """
        Create a STRAPTEA object for techno-economic analysis.

        Parameters
        ----------
        system : System
            Should contain feed and product streams.
        IRR : float
            Internal rate of return (fraction).
        duration : tuple[int, int]
            Start and end year of venture (e.g. (2018, 2038)).
        depreciation : str
            'MACRS' + number of years (e.g. 'MACRS7').
        operating_days : float
            Number of operating days per year.
        income_tax : float
            Combined federal and state income tax rate (fraction).
        lang_factor : float
            Lang factor for getting fixed capital investment from
            total purchase cost. If no lang factor, estimate capital investment
            using bare module factors.
        startup_schedule : tuple[float]
            Startup investment fractions per year
            (e.g. (0.5, 0.5) for 50% capital investment in the first year and 50%
            investment in the second).
        WC_over_FCI : float
            Working capital as a fraction of fixed capital investment.
        labor_cost : float
            Total labor cost (USD/yr).
        fringe_benefits : float
            Cost of fringe benefits as a fraction of labor cost.
        property_tax : float
            Fee as a fraction of fixed capital investment.
        property_insurance : float
            Fee as a fraction of fixed capital investment.
        supplies : float
            Yearly fee as a fraction of labor cost.
        maintenance : float
            Yearly fee as a fraction of fixed capital investment.
        administration : float
            Yearly fee as a fraction of fixed capital investment.

        References
        ----------
        .. [1] Huang, H., Long, S., & Singh, V. (2016). Techno-economic analysis of biodiesel
            and ethanol co-production from lipid-producing sugarcane. Biofuels, Bioproducts
            and Biorefining, 10(3), 299–315. https://doi.org/10.1002/bbb.1640
        """

        __slots__ = (
            "labor_cost",
            "fringe_benefits",
            "maintenance",
            "property_tax",
            "property_insurance",
            "_FCI_cached",
            "supplies",
            "maintanance",
            "administration",
        )

        def __init__(
            self,
            system,
            IRR,
            duration,
            depreciation,
            income_tax,
            operating_days,
            lang_factor,
            construction_schedule,
            WC_over_FCI,
            labor_cost,
            fringe_benefits,
            property_tax,
            property_insurance,
            supplies,
            maintenance,
            administration,
        ):
            super().__init__(
                system,
                IRR,
                duration,
                depreciation,
                income_tax,
                operating_days,
                lang_factor,
                construction_schedule,
                startup_months=0,
                startup_FOCfrac=0,
                startup_VOCfrac=0,
                startup_salesfrac=0,
                finance_interest=0,
                finance_years=0,
                finance_fraction=0,
                WC_over_FCI=WC_over_FCI,
            )
            self.labor_cost = labor_cost
            self.fringe_benefits = fringe_benefits
            self.property_tax = property_tax
            self.property_insurance = property_insurance
            self.supplies = supplies
            self.maintenance = maintenance
            self.administration = administration

        # This returns the ISBL

        def _ENG(self, installed_equipment_cost):
            return (
                installed_equipment_cost
                + installed_equipment_cost * input_data["OSBL_factor"]
            ) * input_data["Engineering"]

        def _CON(self, installed_equipment_cost):
            return (
                installed_equipment_cost
                + installed_equipment_cost * input_data["OSBL_factor"]
            ) * input_data["Contingency"]

        def _DPI(self, installed_equipment_cost):
            return (
                installed_equipment_cost
                + installed_equipment_cost * input_data["OSBL_factor"]
                + (
                    installed_equipment_cost
                    + installed_equipment_cost * input_data["OSBL_factor"]
                )
                * input_data["Engineering"]
                + (
                    installed_equipment_cost
                    + installed_equipment_cost * input_data["OSBL_factor"]
                )
                * input_data["Contingency"]
            )

        def _TDC(self, DPI):
            return DPI

        def _FCI(self, TDC):
            self._FCI_cached = TDC
            return TDC

        def _FOC(self, FCI):
            return FCI * (
                self.property_tax
                + self.property_insurance
                + self.maintenance
                + self.administration
            ) + self.labor_cost + (self.fringe_benefits * self.labor_cost + self.supplies)

    # this stream include the two polymers (film)
    film = bst.Stream("film", C2H4=inlet_pe, units="kg/hr")

    solvent_flow = bst.Stream(
        "solvent_in",
        Dodecane=inlet_solvent,
        units="kg/hr",
        price=input_data["Price"]["solvent"],
    )

    # P1 sep
    (film) - FE101
    (solvent_flow) - V100
    (FE101 - 0, V100 - 0) - V101
    (V101 - 0) - E101
    (E101 - 0) - FHOT
    (FHOT - 0) - P101
    (P101 - 0) - V102
    (V102 - 0) - E102
    (E102 - 0) - P102
    (P102 - 0) - F101
    (F101 - 0) - H101
    (S101 - 1, F101 - 1) - M100
    (M100 - 0) - P106

    # P2 SEP
    (FHOT - 1) - H1012
    (H101 - 1, H1012 - 1) - M101
    (M101 - 0) - E103
    (E103 - 0) - S101
    (H101 - 0, H1012 - 0) - M1003

    pp_sep_sys = bst.System(
        "PP_Recovery_norecycle",
        path=(
            [
                FE101,
                V100,
                V101,
                E101,
                FHOT,
                P101,
                V102,
                E102,
                P102,
                F101,
                H101,
                M101,
                E103,
                S101,
                H101,
                M100,
                P106,
                H1012,
                M1003,
            ]
        ),
        N_runs=N_runs,
    )

    pp_sep_sys.simulate()
    pp_sep_sys.save_report('results_n.xlsx')
    bst.main_flowsheet.diagram(file='FlowDiagram', format='png', number=True)

    pe_recovered_mass_flow = H1012.outs[0].get_flow("kg/hr", "C2H4")
    evoh_recovered_mass_flow = H101.outs[0].get_flow("kg/hr", "C2H4")

    solvent_dissoltion_tank = V101.outs[0].get_flow("kg/hr", "Dodecane")
    plastic_dissolution_tank = V101.outs[0].get_flow("kg/hr", "C2H4")

    solvent_plastic_ratio = solvent_dissoltion_tank/(plastic_dissolution_tank)

    print("===============================")
    print("Process completed.")
    print("Solvent/plastic ratio at dissolution tank: ", solvent_plastic_ratio)

    print("Recovered pe: ", pe_recovered_mass_flow, "kg/hr")
    print("Recovered evoh: ", evoh_recovered_mass_flow, "kg/hr")

    print("===============================")

    strap_tea = STRAPTEA(
        system=pp_sep_sys,
        IRR=input_data["IRR"],
        duration=input_data["duration"],
        depreciation="MACRS7",
        income_tax=input_data["Tax_rate"],
        operating_days=input_data["year_operation_hours"]/24,
        lang_factor=input_data["Lang_factor"],
        construction_schedule=(0.4, 0.6),
        WC_over_FCI=0.10,
        labor_cost=input_data["Shift_operators"] *
        input_data["Salariy_per_operator"],
        fringe_benefits=input_data["Benefits_and_overhead"],
        property_tax=0,
        property_insurance=input_data["Insurance"],
        supplies=0,
        maintenance=input_data["Maintenance"],
        administration=0.00,
    )

    feed = bst.main_flowsheet("PROD")
    msp = strap_tea.solve_price(feed)  # USD/kg
    print("Minimum Selling price: ", int(msp * 1000), "USD/TON")
    print("Total purchase cost: ", int(strap_tea.purchase_cost), "USD")
    print("Total utility cost: ", int(strap_tea.utility_cost), "USD/y")

    # cost of makeup solvent
    solvent_cost = (
        input_data["Price"]["solvent"]
        * inlet_solvent
        * input_data["year_operation_hours"]
    )

    installed_equipment_cost = strap_tea.installed_equipment_cost
    DPI = strap_tea._DPI(installed_equipment_cost)
    TDC = strap_tea._TDC(DPI)
    FCI = strap_tea._FCI(TDC)
    FOC = strap_tea._FOC(FCI)
    # AOC = strap_tea._AOC(AOC)
    ISBL = installed_equipment_cost
    OSBL = installed_equipment_cost * input_data["OSBL_factor"]
    ENG = strap_tea._ENG(installed_equipment_cost)
    CON = strap_tea._CON(installed_equipment_cost)

    electricity = pp_sep_sys.get_electricity_consumption()
    # print("Total electricity consumption (kW-h/y): ", electricity)
    cooling = pp_sep_sys.get_cooling_duty()
    # print("Total cooling duty (kJ/y): ", cooling)
    heating = pp_sep_sys.get_heating_duty()
    # print("Total heating duty (kJ/y): ", heating)

    results = {}
    results["Processing capacity (ton/y)"] = (
        (inlet_pe) * input_data["year_operation_hours"] / 1000
    )
    results['Solvent-plastic ratio at dissolution'] = solvent_plastic_ratio
    # results['Solvent-plastic ratio at dissolution_2'] = solvent_plastic_ratio_2
    # results["Recovered plastic PE (kg/hr)"] = round(PE, 4)
    # results["Recovered plastic PE (kg/hr)"] = round(PE2, 4)
    results["Solvent cost (USD/y)"] = round(solvent_cost, 2)
    # results["Solvent recovery (%)"] = 100-solvent_p_loss
    results["Minimum selling price (USD/kg)"] = round(msp, 4)
    results["ISBL"] = ISBL
    results["installed_equipment_cost"] = installed_equipment_cost
    results["OSBL"] = OSBL
    results["Engineering"] = ENG
    results["Contingency"] = CON
    results["Fixed capital investment"] = FCI
    results["Fixed operating cost"] = FOC
    # results["operating cost"] = AOC

    # changes
    # variable_cost = results["Solvent cost (USD/y)"] + results["Solvent cost_2 (USD/y)"] + results["Solvent cost_3 (USD/y)"] + strap_tea.utility_cost
    variable_cost = results["Solvent cost (USD/y)"] + strap_tea.utility_cost

    results["Variable operating cost"] = variable_cost

    results["Total electricity consumption (kW-h/y)"] = electricity
    results["Total cooling duty (kJ/y)"] = cooling
    results["Total heating duty (kJ/y)"] = heating

    print("Fixed capital investment: {:.2f} USD".format(FCI))
    print("Fixed operating cost: {:.2f} USD/y".format(FOC))
    print("Variable operating cost : {:.2f} USD/y".format(variable_cost))
    print(
        "Contribution of solvent to variable cost: ",
        round(100 * (results["Solvent cost (USD/y)"]) / variable_cost, 2),
        "%",
    )

    makeup_stream = inlet_solvent - P106.outs[0].get_flow("kg/hr", "Dodecane")
    results['Solvent makeup flow (kg/hr)'] = makeup_stream
    print('Solvent makeup flow (kg/hr)', makeup_stream)

    adjusted_solvent_price = makeup_stream * \
        input_data["Price"]["solvent"]/inlet_solvent
    # print('Adjusted solvent price ($/kg)', adjusted_solvent_price)
    

    #######LCA
    functional_unit = (
        inlet_pe * input_data["year_operation_hours"]
    )
    #Total electricity in MJ/y
    El = electricity * 3.6
    #Generated emissions due to electricity 
    Em_electricity = impact_elec * El / functional_unit
    #kg CO2-eq/MJ * MJ/y / kg CO2/y
    results['Climate change impact-electricity (kg CO2-eq/kg polymer)'] = Em_electricity
    print ('Climate change impact-electricity (kg CO2-eq/kg polymer)', Em_electricity)
    #Total steam in MJ/y
    St = heating / 1000
    #Generated emissions due to steam 
    Em_steam = impact_steam * St / functional_unit
    #kg CO2-eq/MJ * MJ/y / kg CO2/y
    results['Climate change impact-steam (kg CO2-eq/kg polymer)'] = Em_steam
    print ('Climate change impact-steam (kg CO2-eq/kg polymer)', Em_steam)
    # flow rate of chilled water in kmol/yr
    flow_water = pp_sep_sys.get_utility_flow("chilled_water")
    # flow rate of chilled brine in kmol/yr
    flow_brine = pp_sep_sys.get_utility_flow("chilled_brine")
    #Total water in m3/y
    Wa = (flow_water+flow_brine) * 18 / 1000
    #Generated emissions due to water consumption (assuming 20% make up)
    Em_wa = impact_water * Wa * 0.02 / functional_unit  
    results['Climate change impact-water (kg CO2-eq/kg polymer)'] = Em_wa
    print ('Climate change impact-water (kg CO2-eq/kg polymer)', Em_wa)
    #Total climate change impact in kg CO2-eq/kg film
    Emissions = Em_electricity + Em_steam + Em_wa
    results['Climate change impact (kg CO2-eq/kg polymer)'] = Emissions
    print ('Climate change impact (kg CO2-eq/kg polymer)', Emissions)

    with open("_results_simulation.json", "w") as fp:
        json.dump(results, fp, indent=4)

    
    return adjusted_solvent_price


def main_simulation(input_data, real_price):
    # calculate adjusted price using real price
    adjusted_price = run_simulation(input_data, real_price)

    # re-calculate costs with adjusted price
    run_simulation(input_data, adjusted_price)

# ====================================================== run simulation


# Import input_data.json with main information about the process.
with open('input_data.json', 'r') as f:
    input_data = json.load(f)

main_simulation(input_data, input_data["Price"]["solvent"])

