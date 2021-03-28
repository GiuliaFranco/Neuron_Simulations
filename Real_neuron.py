from neuron import h, gui   # Standard "import" of the NEURON library into Python...
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Simulate_Real_Neuron:
    '''
    This class is responsible for the creation of neuron's morphology
    and dynamics injection.
    '''

    def __init__(self):
        h.load_file('stdrun.hoc')
        self.soma=''     ### soma variable
        self.axon=''     ### axon variable
        self.V=self.m=self.h=self.n=self.s=self.t=[] ### Voltage, m,n,h,s gate values, time recording. 
        self.model_name=''   ### dynamic name choosen 
        self.spike_t=self.lat=[]    ### spike time [s], latency [ms]
        
    def createSoma(self,cm=0.5,diam=500,L=100,Ra=35.4):
        '''
        Define soma biological structure.
        '''
        self.soma = h.Section(name='soma')
        self.soma.cm = cm
        self.soma.diam = diam
        self.soma.L=L
        self.soma.Ra=Ra
        h.psection(sec=self.soma)
        
    def createAxon(self,nseg=100,diam=2,L=200,Ra=123):
        '''
        Define axon biological structure.
        '''
        self.axon=h.Section(name='axon')
        self.axon.nseg=nseg
        self.axon.Ra=Ra
        self.axon.L=L
        self.axon.diam=diam
        h.psection(sec=self.axon)
        
    def createMorphology(self):
        '''
        Auxiliary function for creating morphology
        '''
        self.createSoma()
        self.createAxon()
        
    def AdjustNumberChannels(self,sec):
        '''
        Auxiliary function for adjusting number of channels 
        in neuron's compartment
        '''
        if(self.model_name=='WBCN'):
            for i in sec:
                i.WBCN.NK*=0.1
                i.WBCN.NNa*=0.1
            if(sec==self.axon):
                self.axon(0.005).WBCN.NNa*=100
                self.axon(0.005).WBCN.NK*=100

    def AuxiliaryInjectDynamics(self,flag):
        '''
        Auxiliary function for injecting dynamics
        '''
        for sec in h.allsec():
            sec.insert(flag)
            self.AdjustNumberChannels(sec)
        
    def InjectDynamics(self,flag):
        '''
        Inject the desired dynamics inside neuron's compartment:
            - "WB": Basic Wang Buzsaki model
            - "WBS": Wang Buzsaki model + slow sodium inactivation
            - "WBCN": Wang Buzsaki model + slow sodium inactivation + Guler stochasticity
        '''
        try:
            self.model_name=flag
            self.AuxiliaryInjectDynamics(flag)
            self.axon.connect(self.soma,1)
            h.topology()
        except:
            print('No dynamic module named as: '+flag)

    def Electrical_Train_Impulse_Soma(self,T=5000,freq=40,Intensity=20.9,dur=0.27303):
        '''
        Create and Run simulation with electrical impulse train in center Soma.
        Setting parameters:
            - T: duration of the electrical impulse train [ms]
            - freq: frequency of electrical impulse train
            - Intensity: the intensity of the current square pulse [micro Ampere]
            - dur: duration of the single the square pulse waveform [ms]
        '''
        int_sec=1000/40 # distance in ms between pulses
        n=int(T/int_sec) # number of pulses
        iclamp=[]   # impulse train
        t_in=[]
        counter=0

        for i in range(0,n):
            imp=h.IClamp(self.soma(0.5))    
            imp.delay=0+counter        
            imp.amp=Intensity
            imp.dur=dur         
            iclamp.append(imp)   
            t_in.append(imp.delay)
            counter=imp.delay+int_sec

        self.Record_Simulation_Electrical(T)
        self.Latency(n,int_sec)
        self.WriteResultsCSV()

    def Record_Simulation_aux(self):
        '''
        Auxiliary function for recording dynamics
        '''

        if(self.model_name=='WBCN'):
            self.m.record(self.soma(0.5).WBCN._ref_m) 
            self.h.record(self.soma(0.5).WBCN._ref_h) 
            self.n.record(self.soma(0.5).WBCN._ref_n) 
            self.s.record(self.soma(0.5).WBCN._ref_n) 

        if(self.model_name=='WBS'):

            self.m.record(self.soma(0.5).WBS._ref_m) 
            self.h.record(self.soma(0.5).WBS._ref_h) 
            self.n.record(self.soma(0.5).WBS._ref_n) 
            self.s.record(self.soma(0.5).WBS._ref_s)

        if(self.model_name=='WB'):
            self.m.record(self.soma(0.5).WB._ref_m) 
            self.h.record(self.soma(0.5).WB._ref_h) 
            self.n.record(self.soma(0.5).WB._ref_n) 

    def Record_Simulation_Electrical(self,T,dt=0.01):
        '''
        Run & Recording dynamics
        '''
        h.dt=dt
        h.psection()
        self.V = h.Vector()             # Membrane potential vector is created here
        self.t = h.Vector() 
        self.m = h.Vector() 
        self.n = h.Vector()
        self.h = h.Vector()  
        self.s = h.Vector()
        self.V.record(self.soma(0.5)._ref_v)
        self.Record_Simulation_aux()
        self.t.record(h._ref_t)   
        h.v_init = -60                # Let's set the initial condition of the membrane potential
        h.t     =   0.0               # Let's reset the initial time of the simulation to 0 ms
        h.tstop = float(T)   
        h.run()


    def Latency(self,n,int_sec):
        v=np.array(self.V).tolist()
        tem=np.array(self.t)
        t1=tem.tolist()
        t_out=[]
        self.lat=[]
        self.spike_t=[]
        for i in range(0,n):
            indx1=np.where(tem > i*int_sec)[0][0]
            indx2=np.where(tem <(i*int_sec+int_sec))[0][-1]
            aux=max(v[indx1:indx2])
            if(aux>5):
                aux_t=int(t1[v.index(aux)]/int_sec)*int_sec
                self.lat.append(t1[v.index(aux)]-aux_t)
                self.spike_t.append(aux_t/1000)

    def WriteResultsCSV(self):
        df = pd.DataFrame({'time':self.t,'V':self.V,'m':self.m,'n':self.n,'h':self.h})             
        if(len(self.s)!=0):
            df['s']=self.s
        df.to_csv('Gates.csv') #### Writes voltage & gates results to file

        df_1 = pd.DataFrame({'Spike_time':self.spike_t,'Latency':self.lat})             
        df_1.to_csv('Latency.csv')  #### Writes spike time & Latency results to file

sup=Simulate_Real_Neuron()   ## create empty cell
sup.createMorphology()       ## create cell morphology (now soma+axon)
sup.InjectDynamics('WBCN')   ## create cell dynamic 
sup.Electrical_Train_Impulse_Soma(20000)    ## Run simulation recording V,n,m,s (now with electrical impulse in soma, recording soma) 


##### EXAMPLE OF ACCESS TO VARIABLES: Print voltage vs time recorded ######
plt.plot(sup.t, sup.V,label="soma")
plt.xlabel('time (ms)')
plt.title('Delay record along the axon')
plt.ylabel('Voltage (mV)')
plt.legend(loc="upper right")
plt.show()

##### EXAMPLE OF ACCESS TO VARIABLES: Print latency vs spike time recorded ######
plt.plot(sup.spike_t, sup.lat, '.', color='darkred', label="soma")
plt.xlabel('time (s)')
plt.ylabel('Latency (ms)')
plt.title(' Soma Latency response phases freq=40 Hz')
plt.legend(loc="lower right")
plt.show()
