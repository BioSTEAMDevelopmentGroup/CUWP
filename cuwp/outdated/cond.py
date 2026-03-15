# Import the relevant libraries 
import biosteam as bst
from biosteam import units
import thermosteam as tmo
from biosteam import Stream, settings
import numpy as np

# Create a function to simulate this process
def process(Temp):
    # Set the thermodynamics of the inlet stream and specify the phase
    chemicals = bst.Chemicals([bst.Chemical("C2H4", phase="s"), "Dodecane", "Water", "N2", "O2", "CH4", "CO2"], cache=True)  ## Change
    bst.settings.set_thermo(chemicals)

    # Specify the operating temperatures of the unit operations
    cooling_temp_after_precipitation = Temp
    temp_dryer = 373.15

    # Create the streams and provide the relevant parameters 
    feed = bst.Stream('feed', Dodecane=800, C2H4=200, units="kg/hr")                                                        ## Change
    H101 = units.DrumDryer(ID='H101',thermo=None, split={'C2H4':0, 'Dodecane':1}, R= 1.4, H=20,length_to_diameter=25,          ## Change
                        T=temp_dryer, moisture_content=0, utility_agent="Natural gas")
    E102 = units.HXutility(ID='E102', T=cooling_temp_after_precipitation, rigorous=True)
    S101 = units.PhaseSplitter(ID='S101', outs=('Vapour_solvent_A', 'Liquid_solvent_A'))

    # connect the streams using the pipe notation
    (feed) - H101
    (H101 - 1) - E102
    (E102) - S101

    # Create the system and specify the path 
    pp_cond_sys = bst.System("pp_condenser", path=([H101, E102, S101]))

    # simulate the process
    pp_cond_sys.simulate()
    rec_solvent = pp_cond_sys.outs[3].get_flow('kmol/hr', 'Dodecane')                                         ## Change
    feed_solv = pp_cond_sys.ins[0].get_flow('kmol/hr', 'Dodecane')                                               ## Change
    return rec_solvent, feed_solv

# Get the required condenser temperature
def condenser_temp (recovery_percentage):
    Temperature = np.arange(70+273.15, -100+273.15, -10)
    outlets = np.zeros(len(Temperature))
    inlets = np.zeros(len(Temperature))
    for i,j in enumerate(Temperature):
        outlets[i] = (process(j)[0])
        inlets[i] = (process(j)[1])
        percentage_recovery = (outlets[i]/inlets[i]) * 100
        if percentage_recovery >= recovery_percentage:
            temp_condenser = j
            break
    return temp_condenser

Temp = float(condenser_temp(99.99))