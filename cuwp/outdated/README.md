Only two files are need now: input_data and run

### Input data

Process data must be input with `input_data.json`. 

### Estimate the T of condensor

This is done with the file `cond.py` . The file `run_.py` calls `cond.py`

### Run simulation and TEA

Run the script `run_.py` to simulate the process and find the minimum selling price of the recovered plastic. Results are saved in `_results_simulation.json`. `FlowDiagram.png` shows the process flow diagram generated and `results_n.xlsx`. shows the detailed results for the streams and costs.  

We are running the process without recycling. To account for the cost of only the makeup solvent (assuming the solvent output stream from Pump P106 would be recycled), I estimated a "new" price for the solvent. I called this "adjusted_solvent_price" and it was estimated considering the makeup solvent (obtained running once the simulation), the "real" price for the solvent (which can be found in input_data) and the inlet solvent. 