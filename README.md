# Neuron_Simulations


Repository for Neuron simulation code. 


#### Run instruction for real neuron python code

In order to run simulation I recommend using Python3. You can check your current python version typing in your terminal:

```
python --version
```

Clone the repository locally typing in your terminal:

```
git clone https://github.com/GiuliaFranco/Neuron_Simulations.git

```

Enter in the repository folder. First we need to install all the required packages for the program.Type in your terminal (in Neuron_Simulation folder):

```
pip3 install -r requirements.txt

```
Now it's time to compile the .mod files needed for the real neuron dynamics.Type in your terminal (in Neuron_Simulation folder):

```
nrnivmodl mods
```
You should notice a new folder created called x86_64.
Finally run the code with:

```
python3 main.py
```

#### Real Neuron - Electrical & Light Simulations

So far only electrical simulation on real neuron are implemented. Inside Real_neuron.py you can create the neuron structure 'soma+axon' specifying the biological constants as you wish. Moreover, the code is implemented for accepting 3 kind of biophysical dynamics:

- 'WB': Wang Buzsaki model --> see WB.mod file in mods folder for implementation details
- 'WBS': Wang Buzsaki model + slow sodium inactivation --> see WBS.mod file in mods folder for implementation details
- 'WBCN': Wang Buzsaki model + slow sodium inactivation + Guler stochasticity --> see WBCN.mod file in mods folder for implementation details

The model is specified by the user as input of a function.

##### Slow Sodium Note: In this version the slow sodium only one time scale is present (parametes in the s differential equation in .mod file are constants). You can eventually change them if you wish by modifying WBS.mod file or WBCN.mod file according to the model implemented.REMEMBER to run 'nrnivmodl mods' each time you change the .mod files, in order to make you changes effective in simulations.

Simulations are implemented injecting an electrical/light input train (variables can be defined by the user, some defaults one are already present in the code). Recording are performed all along the neuron (code is commented and highlighted just for soma in the light case). Voltage and gates variables are recorderd along the simulations and results are stored in Gates.csv file. Latency is calculated as well with spike time along the simulations and results are stored in Latency.csv file.

##### Simulation duration Note: as you can notice in the code, the duration of the entire simulations is the same of the duration input electrical train, which means that stimulations persist all along simualtions period. We can change this parameter if you prefer.
