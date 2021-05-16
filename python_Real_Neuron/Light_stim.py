from neuron import h, gui   # Standard "import" of the NEURON library into Python...
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Light():

	def __init__(self,cell):
		'''
		Record pack is a 4x6 matrix of lists containing the recording outputs from the neuron. 
		row 0: recordings from soma (0.5)
		row 1: recordings from AIS (axon(0.005))
		row 2: recordings from middle axon (0.50)
		row 3: recordings from end axon (0.99)

		column 0: recordings of V
		column 1: recordings of m
		column 2: recordings of h
		column 3: recordings of n
		column 4: recordings of s
		column 5: recordings of time

		Uncomment the following section for recordings all along the neuron, only soma code is highlighted
		'''
		'''
		self.Record_pack=[[[] for i in range(0,6)] for j in range(0,4)]
		self.spike_t=self.spike_t_AIS=self.spike_t_MA=self.spike_t_EA=[]
		self.lat=self.lat_AIS=self.lat_MA=self.lat_EA=[]    ### spike time [s], latency [ms]
		self.model_name=cell.model_name   ### dynamic name choosen
		self.soma=cell.soma
		self.axon=cell.axon
		self.stimulation_point=''

		'''
		self.Record_pack=[[[] for i in range(0,6)] for j in range(0,1)]
		self.spike_t=[]
		self.lat=[]    ### spike time [s], latency [ms]
		self.model_name=cell.model_name   ### dynamic name choosen
		self.soma=cell.soma
		self.stimulation_point=''

	def StimulationPoint(self,point):
		'''
		stimulation point definition. Commented lines for all the neuron, only soma code is highlighted
		'''
		'''
		if(point=='AIS'):
			self.stimulation_point=self.axon(0.005)
		if(point=='middleAxon'):
			self.stimulation_point=self.axon(0.5050000000000002)
		if(point=='EndAxon'):
			self.stimulation_point=self.axon(0.9950000000000007)
		else:
			self.stimulation_point=self.soma(0.5) 
		'''
		self.stimulation_point=self.soma(0.5) 

	def Inject_Pulse(self,tot,int_sec,wavelength,light_intensity,pulse_width):
		'''
		Definition of the injected train light pulse specifics.Commented lines for all the neuron, only soma code is highlighted
		'''
		self.stimulation_point.ChR2.wavelength=wavelength
		self.stimulation_point.ChR2.pulse_width=pulse_width
		self.stimulation_point.ChR2.light_delay=int_sec
		self.stimulation_point.ChR2.n=int(tot/self.stimulation_point.ChR2.light_delay)
		self.stimulation_point.ChR2.light_intensity=light_intensity
		n=int(self.stimulation_point.ChR2.n)
		#self.axon.connect(self.soma,1) ### axon connection to soma
		h.psection()
		h.topology()
		return n

	def Ligh_Train_Impulse_Soma(self,point='soma',T=5000,freq=25,wavelength=470,pulse_width=0.85,light_intensity=20):
		'''
		Simulation and recordings of optogenetic behavior.Commented lines for all the neuron, only soma code is highlighted
		'''
		self.StimulationPoint(point)
		int_sec=1000/freq # distance in ms between pulses
		n=self.Inject_Pulse(T,int_sec,wavelength,light_intensity,pulse_width)
		self.Record_Simulation_Light(T)
		self.Latency(n,int_sec,self.Record_pack[0][0],self.lat,self.spike_t)
		#self.Latency(n,int_sec,self.Record_pack[1][0],self.lat_AIS,self.spike_t_AIS)
		#self.Latency(n,int_sec,self.Record_pack[2][0],self.lat_MA,self.spike_t_MA)
		#self.Latency(n,int_sec,self.Record_pack[3][0],self.lat_EA,self.spike_t_EA)
		self.WriteResultsCSV()
		
	def Record_Simulation_aux(self):
		'''
        Auxiliary function for recording dynamics
		'''
		#points=[self.soma(0.5),self.axon(0.005),self.axon(0.5050000000000002),self.axon(0.9950000000000007)]
		points=[self.soma(0.5)]
		if(self.model_name=='WBCN'):
			for i in range(0,len(self.Record_pack)):
				self.Record_pack[i]=[h.Vector() for j in self.Record_pack[i]]
				self.Record_pack[i][0].record(points[i]._ref_v)
				self.Record_pack[i][1].record(points[i].WBCN._ref_m) 
				self.Record_pack[i][2].record(points[i].WBCN._ref_h) 
				self.Record_pack[i][3].record(points[i].WBCN._ref_n) 
				self.Record_pack[i][4].record(points[i].WBCN._ref_s) 
				self.Record_pack[i][5].record(h._ref_t)

		if(self.model_name=='WBS'):
			for i in range(0,len(self.Record_pack)):
				self.Record_pack[i]=[h.Vector() for j in self.Record_pack[i]]
				self.Record_pack[i][0].record(points[i]._ref_v)
				self.Record_pack[i][1].record(points[i].WBS._ref_m) 
				self.Record_pack[i][2].record(points[i].WBS._ref_h) 
				self.Record_pack[i][3].record(points[i].WBS._ref_n) 
				self.Record_pack[i][4].record(points[i].WBS._ref_s)
				self.Record_pack[i][5].record(h._ref_t)

		if(self.model_name=='WB'):
			for i in range(0,len(self.Record_pack)):
				self.Record_pack[i]=[h.Vector() for j in self.Record_pack[i]]
				self.Record_pack[i][0].record(points[i]._ref_v)
				self.Record_pack[i][1].record(points[i].WB._ref_m) 
				self.Record_pack[i][2].record(points[i].WB._ref_h) 
				self.Record_pack[i][3].record(points[i].WB._ref_n) 
				self.Record_pack[i][5].record(h._ref_t)

	def Record_Simulation_Light(self,T,dt=0.01):
		'''
        Run & Recording dynamics
		'''
		h.dt=dt
		self.Record_Simulation_aux()
		h.v_init = -60                # Let's set the initial condition of the membrane potential
		h.t     =   0.0               # Let's reset the initial time of the simulation to 0 ms
		h.tstop = float(T)   
		h.run()


	def Latency(self,n,int_sec,V,lat_out,t_out):
		'''
		Latency calculations
		'''
		v=np.array(V).tolist()
		tem=np.array(self.Record_pack[0][5])
		t1=tem.tolist()
		for i in range(0,n):
			indx1=np.where(tem > i*int_sec)[0][0]
			indx2=np.where(tem <(i*int_sec+int_sec))[0][-1]
			aux=max(v[indx1:indx2])
			if(aux>5):
				aux_t=int(t1[v.index(aux)]/int_sec)*int_sec
				lat_out.append(t1[v.index(aux)]-aux_t)
				t_out.append(aux_t/1000)

	def WriteResultsCSV(self):
		'''
		output results: Write just for the soma.Commented lines for all the neuron, only soma code is highlighted
		'''
		names=['soma']
		for i in range(0,len(names)):
			df = pd.DataFrame({'time':self.Record_pack[i][5],'V':self.Record_pack[i][0],'m':self.Record_pack[i][1],'n':self.Record_pack[i][3],'h':self.Record_pack[i][2]})             
			if(len(self.Record_pack[i][4])!=0):
				df['s']=self.Record_pack[i][4]
			df.to_csv('Gates_'+names[i]+'.csv') #### Writes voltage & gates results to file

		'''df_1 = pd.DataFrame({'Spike_time_Soma':self.spike_t,'Latency_Soma':self.lat,
						        'Spike_time_AIS':self.spike_t_AIS,'Latency_AIS':self.lat_AIS,
						         'Spike_time_MA':self.spike_t_MA,'Latency_MA':self.lat_MA,
						         'Spike_time_EA':self.spike_t_EA,'Latency_EA':self.lat_EA})'''   
		df_1 = pd.DataFrame({'Spike_time_Soma':self.spike_t,'Latency_Soma':self.lat})          
		df_1.to_csv('Latency.csv')  #### Writes spike time & Latency results to file
