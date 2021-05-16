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
        self.soma=''     ### soma variable
        self.axon=''     ### axon variable
        self.model_name=''   ### dynamic name choosen 
        
    def createSoma(self,cm=0.5,diam=500,L=100,Ra=35.4):
        '''
        Define soma biological structure.
        '''
        self.soma = h.Section(name='soma')
        self.soma.cm = cm
        self.soma.diam = diam
        self.soma.L=L
        self.soma.Ra=Ra
        #h.psection(sec=self.soma)
        
    def createAxon(self,nseg=100,diam=2,L=200,Ra=123,cm=0.5):
        '''
        Define axon biological structure.
        '''
        self.axon=h.Section(name='axon')
        self.axon.nseg=nseg
        self.axon.Ra=Ra
        self.axon.L=L
        self.axon.diam=diam
        self.axon.cm=cm
        #h.psection(sec=self.axon)
        
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

    def Adjust_S_dynamics(self,sec,flag,As,Bs):

    	if(flag=='WBS' and (As!=0 or Bs!=0)):
    		for i in sec:
    			i.WBS.Phi_as=As
    			i.WBS.Phi_bs=Bs
    			print(i.WBS.Phi_as)

    	if(flag=='WBCN' and (As!=0 or Bs!=0)):
    		for i in sec:
    			i.WBCN.Phi_as=As
    			i.WBCN.Phi_bs=Bs


    def AuxiliaryInjectDynamics(self,flag,opsin_bool=False,As=0,Bs=0):
        '''
        Auxiliary function for injecting dynamics
        '''
        for sec in h.allsec():
            sec.insert(flag)
            if(opsin_bool):
            	sec.insert('ChR2')
            self.AdjustNumberChannels(sec)
            self.Adjust_S_dynamics(sec,flag,As,Bs)
        
    def InjectDynamics(self,flag,opsin_bool=False,As=0,Bs=0):
        '''
        Inject the desired dynamics as flag inside neuron's compartment:
            - "WB": Basic Wang Buzsaki model
            - "WBS": Wang Buzsaki model + slow sodium inactivation
            - "WBCN": Wang Buzsaki model + slow sodium inactivation + Guler stochasticity

        Set Opsin_bool as True in order to include ChR2 dynamic in the neuron. Dafault value is False
        '''
        try:
        	self.model_name=flag
        	self.AuxiliaryInjectDynamics(flag,opsin_bool,As,Bs)
        except:
     	   print('No dynamic module named as: '+flag)

