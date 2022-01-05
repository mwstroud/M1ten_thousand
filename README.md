# M1ten_thousand

Matt Stroud, Tyler Banks, and Feng Feng
A ten thousand cell network of the primary motor cortex with 4 unique cell types in layer 5A and 5B: CP, CS, FSI, and LTS. Utilizes distance dependent connections for more realism.


There are 3 experiments to run:
* baseline        (No 50 Hz input)
* short burst     (50 Hz input in 125 ms bursts)
* long burst      (50 Hz input in 1000 ms bursts) 

They are defined in the "simulation_config_... .json" files

To switch which experiment is run, edit the "config.json file" to reference the correct simulation .json experiment.

### Building Network

```
sbatch build_batch.sh
```
### Building Input

```
python build_input.py
```

### Running
Compress the whole folder to a zip folder. Run on NSG:
* NEURON on EXPANSE
* Run script is run_network.py

OR...
From command line:
```
sbatch batchfile_newserver.sh
```

