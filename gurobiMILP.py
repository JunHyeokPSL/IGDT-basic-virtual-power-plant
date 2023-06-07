# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 15:03:33 2023

@author: junhyeok
"""
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import time

class gurobi_MILP:
    def __init__(self, NAME, vpp, model_dict, case_dict):
        
        self.m = gp.Model(name=NAME)
        self.vpp = vpp
        
        self.model_dict = model_dict
        self.dayahead_smp = self.model_dict['da_smp']
        
        
        self.wt_list = self.vpp.wt_list
        self.pv_list = self.vpp.pv_list
        self.ess_list = self.vpp.ess_list
        
        self.nWT = self.vpp.nWT
        self.nPV = self.vpp.nPV
        self.nESS = self.vpp.nESS      
        
        self.case_dict = case_dict

        self.UNIT_TIME = self.case_dict['UNIT_TIME'] 
        self.nTimeslot = int (24 / self.UNIT_TIME)
        
    def add_Parameters(self):
        
        if self.case_dict['res_var'] == False:
            self.P_wt = np.zeros([self.nWT, self.nTimeslot])
            self.P_pv = np.zeros([self.nPV, self.nTimeslot])

            for j in range(self.nTimeslot):    
                for i in range(self.nWT):
                    self.P_wt[i,j] = self.wt_list[i].max_power * self.wt_list[i].profile[j]
                for i in range(self.nPV):
                    self.P_pv[i,j] = self.pv_list[i].max_power * self.pv_list[i].profile[j]
    
    def add_Variables(self):
        
        vpp = self.vpp

        self.Pbid = self.m.addVars(self.nTimeslot, vtype =GRB.CONTINUOUS,
                          lb = vpp.total_min_power, ub= vpp.total_max_power, name='Pbid')
        
        # self.Pbid = self.m.addVars(self.nTimeslot, vtype =GRB.CONTINUOUS,
        #                   lb = -10000000, ub= 10000000, name='Pbid')
        
        
        
        self.uncertainty_dict = self.model_dict['uncertainty']
        self.wt_uncertainty = self.uncertainty_dict['wt']
        self.pv_uncertainty = self.uncertainty_dict['pv']
        self.smp_uncertainty = self.uncertainty_dict['smp']
        
        self.P_wt_uncertainty = self.m.addVars(vpp.nWT, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                               lb = -self.wt_uncertainty , ub = self.wt_uncertainty ,
                                               name = 'P_wt_uncertainty')
        self.P_pv_uncertainty = self.m.addVars(vpp.nPV, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                               lb = - self.pv_uncertainty, ub = self.pv_uncertainty,
                                               name = 'P_pv_uncertainty')
        
        self.Smp_uncertainty = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS,
                                               lb = - self.smp_uncertainty, ub = self.smp_uncertainty,
                                               name = 'Smp_uncertainty')
       
            
        
        
        
        if self.case_dict['res_var'] == True:
            
            if self.wt_list:
                self.P_wt = self.m.addVars(vpp.nWT, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.wt_list[i].min_power * (1-self.wt_uncertainty)
                                               for _ in range(self.nTimeslot)] for i in range(self.nWT)],
                                          ub=[[self.wt_list[i].max_power * (1+self.wt_uncertainty)
                                               for _ in range(self.nTimeslot)] for i in range(self.nWT)],
                                          name='P_wt'
                                          )  
            if self.pv_list:
                self.P_pv = self.m.addVars(vpp.nPV, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.pv_list[i].min_power * (1-self.pv_uncertainty)
                                               for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                          ub=[[self.pv_list[i].max_power * (1+self.pv_uncertainty)
                                               for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                          name='P_pv'
                                          )
                
            except Exception as e:
                print("")
                print(e)
                print("gurobi_MILP add Variables")
                print("No Uncertainty Sets in this case")
                print("")
                if self.wt_list:
                    self.P_wt = self.m.addVars(vpp.nWT, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                              lb=[[self.wt_list[i].min_power for _ in range(self.nTimeslot)] for i in range(self.nWT)],
                                              ub=[[self.wt_list[i].max_power for _ in range(self.nTimeslot)] for i in range(self.nWT)],
                                              name='P_wt'
                                              )  
                if self.pv_list:
                    self.P_pv = self.m.addVars(vpp.nPV, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                              lb=[[self.pv_list[i].min_power for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                              ub=[[self.pv_list[i].max_power for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                              name='P_pv'
                                              )
    
        if self.ess_list:
            self.P_essChg = self.m.addVars(vpp.nESS, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                     lb= 0,
                                     ub=[[-self.ess_list[i].min_power for _ in range(self.nTimeslot)] for i in range(self.nESS)],
                                     name='P_essChg'
                                     )
            self.P_essDis = self.m.addVars(vpp.nESS, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                     lb= 0,
                                     ub=[[self.ess_list[i].max_power for _ in range(self.nTimeslot)] for i in range(self.nESS)],
                                     name='P_essDis'
                                     )
            self.U_essChg = self.m.addVars(vpp.nESS, self.nTimeslot, vtype =GRB.BINARY, name='U_essChg')
            self.U_essDis = self.m.addVars(vpp.nESS, self.nTimeslot, vtype =GRB.BINARY, name='U_essDis')
            
        
            
            
    def add_bid_constraints(self):
        for j in range(self.nTimeslot):
            self.m.addConstr(self.Pbid[j] == gp.quicksum(self.P_wt[i,j] for i in range(self.nWT))
                             + gp.quicksum(self.P_pv[i,j] for i in range(self.nPV))
                             + gp.quicksum(self.P_essDis[i,j] for i in range(self.nESS))
                            - gp.quicksum(self.P_essChg[i,j] for i in range(self.nESS))
                            , name=f'const_bid{j}')
                             
    def add_res_constraints(self):
        
        if self.case_dict['res_var'] == True:
            
            for j in range(self.nTimeslot):
                
                if self.case_dict['case'] == 2:

                    if self.case_dict['bid_type'] == 'risk_averse':

                        for i in range(self.nWT):
                            self.m.addConstr(self.P_wt[i,j] == (1 - self.P_wt_uncertainty[i,j] ) *                                       
                                             self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                            name=f'const_wt{i}_{j}_profile_uncertainty')
        
                        for i in range(self.nPV):
                            self.m.addConstr(self.P_pv[i,j] ==(1 - self.P_pv_uncertainty[i,j] ) *                                           
                                             self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                            name=f'const_pv{i}_{j}_profile_uncertainty')
                        
                    else:
                        for i in range(self.nWT):
                            self.m.addConstr(self.P_wt[i,j] == (1 + self.P_wt_uncertainty[i,j] ) *                                       
                                             self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                            name=f'const_wt{i}_{j}_profile_uncertainty')
        
                        for i in range(self.nPV):
                            self.m.addConstr(self.P_pv[i,j] ==(1 + self.P_pv_uncertainty[i,j] ) *                                           
                                             self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                            name=f'const_pv{i}_{j}_profile_uncertainty')
                        
                else:
                    for i in range(self.nWT):
                        self.m.addConstr(self.P_wt[i,j] == self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                        name=f'const_wt{i}_{j}_profile')
    
                    for i in range(self.nPV):
                        self.m.addConstr(self.P_pv[i,j] == self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                        name=f'const_pv{i}_{j}_profile')
        
    def add_ess_contraints(self):
        # ESS
        ess_list = self.ess_list
        if ess_list:
            for i in range(self.nESS):

                for j in range(self.nTimeslot):    
                    
                    # Ess minmax with preventing Discharging and Charging simultaneously
            
                    self.m.addConstr(self.P_essDis[i, j] <= self.U_essDis[i, j] * ess_list[i].max_power,
                                    name=f'const_ess{i}_{j}_power_max')
                    self.m.addConstr(self.P_essChg[i, j] <= self.U_essChg[i, j] * - ess_list[i].min_power,
                                    name=f'const_ess{i}_{j}_power_min')
                    self.m.addConstr(self.U_essDis[i, j] + self.U_essChg[i, j] <= 1,
                                    name=f'const_ess{i}_{j}_on/off')
                    
                    # SoC min max constraint
                    self.m.addConstr(ess_list[i].initSOC 
                                     - sum(self.P_essDis[i, k] * self.UNIT_TIME 
                                           / ess_list[i].efficiency / ess_list[i].max_capacity
                                                for k in range(j + 1)) 
                                     + sum(self.P_essChg[i, k] * self.UNIT_TIME 
                                           * ess_list[i].efficiency / ess_list[i].max_capacity
                                                for k in range(j + 1)) <= ess_list[i].maxSOC, 
                                          name=f'const_ess{i}_{j}_soc_max')
                    
                    self.m.addConstr(ess_list[i].initSOC
                                          - sum(self.P_essDis[i, k] * self.UNIT_TIME 
                                                / ess_list[i].efficiency / ess_list[i].max_capacity
                                                for k in range(j + 1))
                                          + sum(self.P_essChg[i, k] * self.UNIT_TIME 
                                                * ess_list[i].efficiency / ess_list[i].max_capacity
                                                for k in range(j + 1)) >= ess_list[i].minSOC,
                                          name=f'const_ess{i}_{j}_soc_min')
                # SoC terminal constraint
                self.m.addConstr(ess_list[i].initSOC
                                      - sum(self.P_essDis[i, k] * self.UNIT_TIME 
                                            / ess_list[i].efficiency / ess_list[i].max_capacity
                                            for k in range(self.nTimeslot))
                                      + sum(self.P_essChg[i, k] * self.UNIT_TIME 
                                            * ess_list[i].efficiency / ess_list[i].max_capacity
                                            for k in range(self.nTimeslot)) == ess_list[i].termSOC,
                                     name=f'const_ess{i}_term')  
    def set_Objectives(self):
        
        self.obj = gp.quicksum(self.dayahead_smp[j]*self.Pbid[j] for j in range(self.nTimeslot))
        #self.obj = gp.quicksum(self.dayahead_smp[j] for j in range(self.nTimeslot))

        self.set_obj = self.m.setObjective(self.obj, GRB.MAXIMIZE)
        #self.set_obj = self.m.setObjective(self.obj, GRB.MINIMIZE)
    
        return self.obj
    
    def solve(self, tol, timelimit=None):
        
        mip_gap = tol[0]
        feas_tol = tol[1]
        self.m.setParam(GRB.Param.MIPGap, mip_gap)
        self.m.setParam(GRB.Param.FeasibilityTol, feas_tol)
        
        time_start_op = time.time()
        sol = self.m.optimize()
        time_end_op = time.time()
        
        print("Optimization Duration Time:", time_end_op - time_start_op)
        
        if self.m.status == GRB.OPTIMAL:
            print("Optimal Solution:")
        elif self.m.status == GRB.INFEASIBLE:
            print("Model is infeasible")
            
            self.m.computeIIS()
            if self.m.IISMinimal:
                print("Model is infeasible.")
                infeasible_constrs = self.m.getConstrs()
                infeasible_vars = self.m.getVars()
                for constr in infeasible_constrs:
                    if constr.IISConstr:
                        print(f"Infeasible constraint: {constr.ConstrName}")
                for var in infeasible_vars:
                    if var.IISLB > 0 or var.IISUB > 0:
                        print(f"Infeasible variable: {var.VarName}")

        return sol, self.obj
        
    def get_sol(self):
        
        P_BidSol = np.zeros([self.nTimeslot])
        
        P_essDisSol = np.zeros([self.nESS, self.nTimeslot])
        P_essChgSol = np.zeros([self.nESS, self.nTimeslot])
        U_essDisSol = np.zeros([self.nESS, self.nTimeslot])
        U_essChgSol = np.zeros([self.nESS, self.nTimeslot])
        if self.case_dict['res_var'] == True:
            P_wtSol = np.zeros([self.nWT, self.nTimeslot])
            P_pvSol = np.zeros([self.nPV, self.nTimeslot])
        
        for j in range(self.nTimeslot):
            P_BidSol[j] = self.m.getVarByName(f"Pbid[{j}]").X
            for i in range(self.nESS):
                P_essDisSol[i,j] = self.m.getVarByName(f"P_essDis[{i},{j}]").X
                P_essChgSol[i,j] = self.m.getVarByName(f"P_essChg[{i},{j}]").X
                U_essDisSol[i,j] = self.m.getVarByName(f"U_essDis[{i},{j}]").X
                U_essChgSol[i,j] = self.m.getVarByName(f"U_essChg[{i},{j}]").X
                
                if self.case_dict['res_var'] == True:
                    for i in range(self.nWT):
                        P_wtSol[i,j] = self.m.getVarByName(f"P_wt[{i},{j}]").X
                    for i in range(self.nPV):
                        P_pvSol[i,j] = self.m.getVarByName(f"P_pv[{i},{j}]").X
        
        P_dict = {'bid': P_BidSol, 'essDis': P_essDisSol, 'essChg': P_essChgSol}
        U_dict = {'essDis': U_essDisSol, 'essChg': U_essChgSol}
        
        if self.case_dict['res_var'] == True:
            P_dict['wt'] = P_wtSol
            P_dict['pv'] = P_pvSol
            
        return P_dict, U_dict
