from neuron import h, gui   # Standard "import" of the NEURON library into Python...
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import Real_neuron as cell 
import Electrical_stim as stimEl
import Light_stim as stimLight

h.load_file('stdrun.hoc')

sup=cell.Simulate_Real_Neuron()   ## create empty cell
sup.createMorphology()       ## create cell morphology (now soma+axon)
sup.InjectDynamics('WBCN',opsin_bool=True)   ## create cell dynamic 


#stimulation=stimEl.Electrical(sup)
#stimulation.Electrical_Train_Impulse_Soma(T=1000)

stimulation=stimLight.Light(sup)
stimulation.Ligh_Train_Impulse_Soma(T=30000)


##### EXAMPLE OF ACCESS TO VARIABLES: Print voltage vs time recorded ######
plt.plot(stimulation.Record_pack[0][5], stimulation.Record_pack[0][0],label="soma")
#plt.plot(stimulation.Record_pack[3][5], stimulation.Record_pack[3][0],label="End axon")
plt.xlabel('time (ms)')
plt.title('Delay record along the axon')
plt.ylabel('Voltage (mV)')
plt.legend(loc="upper right")
plt.show()


##### EXAMPLE OF ACCESS TO VARIABLES: Print latency vs spike time recorded ######
plt.plot(stimulation.spike_t, stimulation.lat, '.', color='darkred', label="soma")
plt.xlabel('time (s)')
plt.ylabel('Latency (ms)')
plt.title(' Soma Latency response phases freq=40 Hz')
plt.legend(loc="lower right")
plt.show()
