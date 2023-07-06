# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 15:03:24 2023

@author: junhyeok
"""

import numpy as np
import pandas as pd

# Define the aggrergator class

class aggregator:
    def __init__(self, agg, char_ess, model_dict):
        self.name = agg['name']
        self.code = agg['code']
        self.ngen_dict = agg['gen']
        self.nTotalGen = 0
        
        self.wt_profile = agg['profile'][0]
        self.pv_profile = agg['profile'][1]
        
        self.char_ess = char_ess
        
        self.wt_list = []
        self.pv_list = []
        self.ess_list = []
        
        self.model_dict = model_dict
        self.nTimeslot = self.model_dict['nTimeslot']
        
        try:
            self.uncertainty_dict = self.model_dict['uncertainty']
            self.wt_uncert = self.uncertainty_dict['wt']
            self.pv_uncert = self.uncertainty_dict['pv']
        except Exception as e:
            print("Error")
            print(e)
            print("Class Aggregator")
            print("No Uncertainty Sets in this cases")
            print("No Uncertainty Sets in this cases")
            print("No Uncertainty Sets in this cases")
        self.initialize_res()
    
    def initialize_res(self):
        
        self.nWT = self.ngen_dict['WT']
        self.nPV = self.ngen_dict['PV']
        self.nESS = self.ngen_dict['ESS']
        self.nTotalGen = self.nPV + self.nWT + self.nESS
        
        count = 1
        
        for i in range(self.nWT):
            self.wt_list.append(WT(self.name, self.code, count, self.wt_profile))
            count = count + 1
        
        for i in range(self.nPV):
            self.pv_list.append(PV(self.name, self.code, count, self.pv_profile))
            count = count + 1
            
        for i in range(self.nESS):
            self.ess_list.append(ESS(self.name, self.code,count, self.char_ess))
            count = count + 1
            
    def set_wt_power(self, max_power_list):
        
        for i in range(len(self.wt_list)):
            self.wt_list[i].set_power(max_power_list[i])
        
    def set_pv_power(self, max_power_list):
        for i in range(len(self.pv_list)):
            self.pv_list[i].set_power(max_power_list[i])
    
    def set_ess_power(self, max_power_list):
        for i in range(len(self.ess_list)):
            self.ess_list[i].set_power(max_power_list[i])
            
    def set_ess_capacity(self,max_capacity_list):
        for i in range(len(self.ess_list)):
            self.ess_list[i].set_capacity(max_capacity_list[i])

    def set_der_power(self,max_list):
        self.set_wt_power(max_list[0])
        self.set_pv_power(max_list[1])
        self.set_ess_power(max_list[2])
        self.set_ess_capacity(max_list[3])
        self.set_res_table()
    
    def set_res_table(self):
        
        data_list = []
        self.total_max_power = np.zeros(self.nTimeslot) 
        self.total_min_power = np.zeros(self.nTimeslot)
        res_list = [self.wt_list, self.pv_list, self.ess_list]
        
        try:
            uncertainty_list = [self.wt_uncert, self.pv_uncert, np.zeros(self.nTimeslot)]
            for i in range(len(res_list)):
                for j in range(len(res_list[i])):
                    data_list.append(res_list[i][j].get_res_data())
                    for step in range(self.nTimeslot):
                        self.total_max_power[step] += res_list[i][j].max_power * (1+uncertainty_list[i][step])
                        self.total_min_power[step] += res_list[i][j].min_power * (1-uncertainty_list[i][step])
                    
        except Exception as e:
            print("Error")
            print(e)
            
            print("Aggregator set_res_table method")
            print("Uncertainty does not exist")
            print("Uncertainty does not exist")
            print("Uncertainty does not exist")
            for i in range(len(res_list)):
                for j in range(len(res_list[i])):
                    data_list.append(res_list[i][j].get_res_data())
                    for step in range(self.nTimeslot):
                        self.total_max_power[step] += res_list[i][j].max_power 
                        self.total_min_power[step] += res_list[i][j].min_power           
            
        if self.ess_list:
            self.res_table = pd.DataFrame(data_list,columns=['name', 'type', 'number', 'min_power', 'max_power','capacity'])
        else:
            self.res_table = pd.DataFrame(data_list,columns=['name', 'type', 'number', 'min_power', 'max_power'])
    def get_res_table(self):
        
        try:
            return self.res_table
        except:
            
            self.set_res_table()
            return self.res_table
            
    
class WT:
    def __init__(self, name,code, count, wt_profile):
        self.name = f'WT{count}_{name}'
        self.type = 'WT'
        self.cvpp_name = name
        self.cvpp_code = code
        self.busNumber = count
        self.min_power = 0
        self.max_power = 0
        self.profile = wt_profile
        
    def set_power(self, max_power):
        # Unit [kW]
        self.max_power = max_power
    
    def get_res_data(self):
        self.res_data = [
            self.name,
            self.type,
            self.busNumber,
            self.min_power,
            self.max_power
        ]  
        return self.res_data 
            
                
        
class PV:
    def __init__(self, name,code, count, pv_profile):
        self.name = f'PV{count}_{name}'
        self.type = 'PV'
        self.cvpp_name = name
        self.cvpp_code = code
        self.busNumber = count
        self.min_power = 0
        self.max_power = 0
        self.profile = pv_profile
        
    def set_power(self, max_power):
        # Unit [kW]
        self.max_power = max_power
        
    def get_res_data(self):
        self.res_data = [
            self.name,
            self.type,
            self.busNumber,
            self.min_power,
            self.max_power
        ]  
        return self.res_data 
        
class ESS:
    def __init__(self, name,code, count, ess):
        self.name = f'ESS{count}_{name}'
        self.type = 'ESS'
        self.cvpp_name = name
        self.cvpp_code = code
        self.busNumber = count
        self.min_power = 0
        self.max_power = 0
        self.initSOC = ess['initSOC']
        self.termSOC = ess['termSOC']
        self.minSOC = ess['minSOC']
        self.maxSOC = ess['maxSOC']
        self.max_capacity = 0
        self.efficiency = ess['efficiency']

    def set_power(self, max_power):
        # Unit [kW]
        self.min_power = - max_power
        self.max_power = max_power
        
    def set_capacity(self, capacity):
        # Unit [kWh]
        self.max_capacity = capacity
                
    def get_res_data(self):
        self.res_data = [
            self.name,
            self.type,
            self.busNumber,
            self.min_power,
            self.max_power,
            self.max_capacity
        ]  
        return self.res_data 