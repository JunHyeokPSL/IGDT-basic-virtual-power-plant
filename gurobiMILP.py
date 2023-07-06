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
        
        self.is_case1 = self.case_dict['case'] == 1
        self.is_case2 = self.case_dict['case'] == 2
        self.is_case3 = self.case_dict['case'] == 3
        self.is_case4 = self.case_dict['case'] == 4
        

        if self.is_case1:
            self.is_res_var = True
            self.is_uncertainty = False
        elif self.is_case2:
            self.is_res_var = True
            self.is_uncertainty = True
        elif self.is_case3:
            self.is_res_var = False
            self.is_uncertainty = False
            self.delta = self.model_dict['delta']
        elif self.is_case4:
            self.is_res_var = True
            self.is_uncertainty = False
            self.delta = self.model_dict['delta']            
            
        else:
            raise Exception("No Considered Case at init is_res_var")
    
        
        
        self.is_case_risk_averse = self.case_dict['bid_type'] == 'risk_averse'
        self.is_igdt_risk_averse = self.case_dict['bid_type'] == 'igdt_risk_averse'
        self.is_igdt_risk_seeking = self.case_dict['bid_type'] == 'igdt_risk_seeking'
        
        
        self.UNIT_TIME = self.case_dict['UNIT_TIME'] 
        self.nTimeslot = int (24 / self.UNIT_TIME)
        
        self.base_obj = 0
        self.delta = 0
        
        self.check_set = {}
        
    def add_Variables(self):
        
        vpp = self.vpp

        if self.is_case1 or self.is_case2: 
            self.Pbid = self.m.addVars(self.nTimeslot, vtype =GRB.CONTINUOUS,
                              lb = [vpp.total_min_power[i] for i in range(self.nTimeslot)],
                              ub= [vpp.total_max_power[i] for i in range(self.nTimeslot)],
                              name='Pbid')
            self.da_smp = np.zeros(self.nTimeslot)
            
        elif self.is_case3:
            self.Pbid = np.zeros(self.nTimeslot)
            self.da_smp = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS, name = 'da_smp')
            self.revenue = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS, name = 'revenue')
            print("does not considered the lb, ub of da_smp")
            
        elif self.is_case4:
            self.Pbid = self.m.addVars(self.nTimeslot, vtype =GRB.CONTINUOUS,
                              name='Pbid')         
            self.da_smp = np.zeros(self.nTimeslot)
            self.revenue = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS, name = 'revenue')
            print("case 4 add Variables")
            
            
        else:
            raise Exception("No Considered Case at add_Variables")
                
            
        if self.is_case1:
            print("")
            print("case 1 add Variables")
            print("gurobi_MILP add Variables")
            print("No Uncertainty Sets in this case")
            print("")
            if self.wt_list:
                self.P_wt = self.m.addVars(vpp.nWT, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.wt_list[i].min_power for step in range(self.nTimeslot)] for i in range(self.nWT)],
                                          ub=[[self.wt_list[i].max_power for step in range(self.nTimeslot)] for i in range(self.nWT)],
                                          name='P_wt'
                                          )  
                print("Assign P_wt in No uncertainty")
                
            if self.pv_list:
                self.P_pv = self.m.addVars(vpp.nPV, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.pv_list[i].min_power for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                          ub=[[self.pv_list[i].max_power for _ in range(self.nTimeslot)] for i in range(self.nPV)],
                                          name='P_pv'
                                          )    
            
                print("Assign P_pv in No uncertainty")
        
        elif self.is_case2:
            self.uncertainty_dict = self.model_dict['uncertainty']
            
            self.wt_uncertainty = self.uncertainty_dict['wt']
            self.pv_uncertainty = self.uncertainty_dict['pv']
            self.smp_uncertainty = self.uncertainty_dict['smp']
            
            if self.is_case_risk_averse:
                self.P_wt_uncertainty = self.m.addVars(vpp.nWT, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [[-self.wt_uncertainty[i] 
                                                              for _ in range(vpp.nWT) for i in range(self.nTimeslot)]],
                                                       ub = [[-self.wt_uncertainty[i] 
                                                             for _ in range(vpp.nWT) for i in range(self.nTimeslot)]],
                                                       name = 'P_wt_uncertainty')
                self.P_pv_uncertainty = self.m.addVars(vpp.nPV, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [[-self.pv_uncertainty[i] 
                                                              for _ in range(vpp.nPV) for i in range(self.nTimeslot)]],
                                                       ub = [[-self.pv_uncertainty[i] 
                                                             for _ in range(vpp.nPV) for i in range(self.nTimeslot)]],
                                                       name = 'P_pv_uncertainty')
                
                self.smp_uncertainty = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [-self.smp_uncertainty[step] for step in range(self.nTimeslot)],
                                                       ub = [self.smp_uncertainty[step] for step in range(self.nTimeslot)],
                                                       name = 'Smp_uncertainty')

            else:
                
                self.P_wt_uncertainty = self.m.addVars(vpp.nWT, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [[-self.wt_uncertainty[i] 
                                                              for _ in range(vpp.nWT) for i in range(self.nTimeslot)]],
                                                       ub = [[self.wt_uncertainty[i] 
                                                             for _ in range(vpp.nWT) for i in range(self.nTimeslot)]],
                                                       name = 'P_wt_uncertainty')
                self.P_pv_uncertainty = self.m.addVars(vpp.nPV, self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [[-self.pv_uncertainty[i] 
                                                              for _ in range(vpp.nPV) for i in range(self.nTimeslot)]],
                                                       ub = [[self.pv_uncertainty[i] 
                                                             for _ in range(vpp.nPV) for i in range(self.nTimeslot)]],
                                                       name = 'P_pv_uncertainty')
                
                self.smp_uncertainty = self.m.addVars(self.nTimeslot, vtype = GRB.CONTINUOUS,
                                                       lb = [-self.smp_uncertainty[step] for step in range(self.nTimeslot)],
                                                       ub = [self.smp_uncertainty[step] for step in range(self.nTimeslot)],
                                                       name = 'Smp_uncertainty')                   
                
            if self.wt_list:
                self.P_wt = self.m.addVars(vpp.nWT, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.wt_list[i].min_power * (1-self.wt_uncertainty[step])
                                               for step in range(self.nTimeslot)] for i in range(self.nWT)],
                                          ub=[[self.wt_list[i].max_power * (1+self.wt_uncertainty[step])
                                               for step in range(self.nTimeslot)] for i in range(self.nWT)],
                                          name='P_wt'
                                          )  
                print("Assign the P_wt")
                
            if self.pv_list:
                self.P_pv = self.m.addVars(vpp.nPV, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          lb=[[self.pv_list[i].min_power * (1-self.pv_uncertainty[step])
                                               for step in range(self.nTimeslot)] for i in range(self.nPV)],
                                          ub=[[self.pv_list[i].max_power * (1+self.pv_uncertainty[step])
                                               for step in range(self.nTimeslot)] for i in range(self.nPV)],
                                          name='P_pv'
                                          )                
                print("Assign the P_pv")
            
        elif self.is_case3:
            print("case 3: Add Parameters")
            self.add_Parameters()
        elif self.is_case4:
            print("")
            print("case 4 add Variables")
            print("gurobi_MILP add Variables")
            print("")
            if self.wt_list:
                self.P_wt = self.m.addVars(vpp.nWT, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          name='P_wt'
                                          )  
                print("Assign P_wt in case 4")
                
            if self.pv_list:
                self.P_pv = self.m.addVars(vpp.nPV, self.nTimeslot, vtype =GRB.CONTINUOUS,
                                          name='P_pv'
                                          )    
                print("Assign P_pv in case 4")
       
        else:
            raise Exception("No Considered Case at add_Variables")
    
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
               
        
        if self.is_igdt_risk_averse or self.is_igdt_risk_seeking:
            
            if self.is_case3 or self.is_case4:
                print("case 3 or 4: alpha is constant")
                self.alpha = self.m.addVar(lb = 0, vtype = GRB.CONTINUOUS, name= 'alpha')
            else:                
                raise Exception("No Considered Case at add_Variables for alpha")
        else:
            print("Does not Cosidered alpha")
                
    
    def add_Parameters(self):
        
        if not self.is_res_var:
            print("add_Parameter res_var False")
            print("add_Parameter res_var False")
            self.P_wt = np.zeros([self.nWT, self.nTimeslot])
            self.P_pv = np.zeros([self.nPV, self.nTimeslot])

            for j in range(self.nTimeslot):    
                for i in range(self.nWT):
                    self.P_wt[i,j] = self.wt_list[i].max_power * self.wt_list[i].profile[j]
                for i in range(self.nPV):
                    self.P_pv[i,j] = self.pv_list[i].max_power * self.pv_list[i].profile[j]
    
    
    def add_bid_constraints(self):
        
        if self.is_case1 or self.is_case2 or self.is_case4: 
            for j in range(self.nTimeslot):
                self.m.addConstr(self.Pbid[j] == gp.quicksum(self.P_wt[i,j] for i in range(self.nWT))
                                 + gp.quicksum(self.P_pv[i,j] for i in range(self.nPV))
                                 + gp.quicksum(self.P_essDis[i,j] for i in range(self.nESS))
                                - gp.quicksum(self.P_essChg[i,j] for i in range(self.nESS))
                                , name=f'const_bid{j}')
        elif self.is_case3:
            for j in range(self.nTimeslot):
                self.Pbid[j] = (
                    sum(self.P_wt[i,j] for i in range(self.nWT))
                    + sum(self.P_pv[i,j] for i in range(self.nPV))
                    )
        else:
            raise Exception("No Considered Case at Bid Constraints")
    
    def add_smp_constraints(self):
        
        if self.is_case1 or self.is_case2 or self.is_case4:
            for j in range(self.nTimeslot):
                self.da_smp[j] = self.dayahead_smp[j]
            print("ADD smp constraints as fixed")
            
        elif self.is_case3:
            for j in range(self.nTimeslot):
                if self.is_igdt_risk_averse:
                    self.m.addConstr(self.da_smp[j] == (1-self.alpha) * self.dayahead_smp[j],
                                     name= f'const_da_smp{j}_uncertainty_igdt_risk_averse') 
                elif self.is_igdt_risk_seeking:
                    self.m.addConstr(self.da_smp[j] == (1+self.alpha) * self.dayahead_smp[j],
                                     name= f'const_da_smp{j}_uncertainty_igdt_risk_seeking')
        else: 
            raise Exception("No Considered Case at smp constraints Constraints")           
                         
    def add_res_constraints(self):
        
        if self.is_res_var:
            
            for j in range(self.nTimeslot):
                
                if self.is_case1:
                    for i in range(self.nWT):
                        self.m.addConstr(self.P_wt[i,j] == self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                        name=f'const_wt{i}_{j}_profile')
    
                    for i in range(self.nPV):
                        self.m.addConstr(self.P_pv[i,j] == self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                        name=f'const_pv{i}_{j}_profile')
                elif self.is_case2:
                    for i in range(self.nWT):
                        self.m.addConstr(self.P_wt[i,j] == (1 + self.P_wt_uncertainty[i,j] ) *                                       
                                         self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                        name=f'const_wt{i}_{j}_profile_uncertainty')
    
                    for i in range(self.nPV):
                        self.m.addConstr(self.P_pv[i,j] ==(1 + self.P_pv_uncertainty[i,j] ) *                                           
                                         self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                        name=f'const_pv{i}_{j}_profile_uncertainty')
                elif self.is_case4:
                    if self.is_igdt_risk_averse:
                        for i in range(self.nWT):
                            self.m.addConstr(self.P_wt[i,j] == (1 - self.alpha) *                                       
                                             self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                            name=f'const_wt{i}_{j}_profile_uncertainty')
        
                        for i in range(self.nPV):
                            self.m.addConstr(self.P_pv[i,j] ==(1) *                                           
                                             self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                            name=f'const_pv{i}_{j}_profile_uncertainty')
                        # for i in range(self.nPV):
                        #     self.m.addConstr(self.P_pv[i,j] ==(1 - self.alpha) *                                           
                        #                      self.pv_list[i].max_power * self.pv_list[i].profile[j],
                        #                     name=f'const_pv{i}_{j}_profile_uncertainty')
                            
                    if self.is_igdt_risk_seeking:
                        for i in range(self.nWT):
                            self.m.addConstr(self.P_wt[i,j] == (1 + self.alpha) *                                       
                                             self.wt_list[i].max_power * self.wt_list[i].profile[j],
                                            name=f'const_wt{i}_{j}_profile_uncertainty')
        
                        for i in range(self.nPV):
                            self.m.addConstr(self.P_pv[i,j] ==(1) *                                           
                                             self.pv_list[i].max_power * self.pv_list[i].profile[j],
                                            name=f'const_pv{i}_{j}_profile_uncertainty')  
                        # for i in range(self.nPV):
                        #     self.m.addConstr(self.P_pv[i,j] ==(1 + self.self.alpha) *                                           
                        #                      self.pv_list[i].max_power * self.pv_list[i].profile[j],
                        #                     name=f'const_pv{i}_{j}_profile_uncertainty') 
                     
                else:                   
                    raise Exception("No Considered Case at add_res_constraints")
                    
            print("add_res_constraints completed")
                
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
    
    def set_igdt_params(self, base_obj, delta):
        self.base_obj = base_obj
        self.delta = delta        
    
    def add_igdt_constraints(self):
        
        if self.is_igdt_risk_averse:
            self.revenue =gp.quicksum(self.da_smp[j]*self.Pbid[j] for j in range(self.nTimeslot))
            self.m.addConstr(self.revenue >= (1-self.delta)*self.base_obj,
                             name = 'const_igdt_risk_averse_cost')
            
        elif self.is_igdt_risk_seeking:
            self.revenue =gp.quicksum(self.da_smp[j]*self.Pbid[j] for j in range(self.nTimeslot))
            self.m.addConstr(self.revenue >= (1+self.delta)*self.base_obj,
                             name = 'const_igdt_risk_seeking_cost')
        
            print("add_igdt_risk_seeking_constraints sucessfully") 
        
        
        
    def set_Objectives(self):
        
        if self.is_case3 or self.is_case4:
            if self.is_igdt_risk_averse:
                self.obj = self.alpha
                self.set_obj = self.m.setObjective(self.obj, GRB.MAXIMIZE)
            elif self.is_igdt_risk_seeking:
                self.obj = self.alpha
                self.set_obj = self.m.setObjective(self.obj, GRB.MINIMIZE)               
            else:                
                raise Exception("No Considered Case at set_Objectives")  
            
        elif self.is_case1 or self.is_case2:
            self.obj = gp.quicksum(self.da_smp[j]*self.Pbid[j] for j in range(self.nTimeslot))
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

        return sol, self.m.objVal
        
    def get_sol(self):
        
        P_BidSol = np.zeros([self.nTimeslot])
        
        P_essDisSol = np.zeros([self.nESS, self.nTimeslot])
        P_essChgSol = np.zeros([self.nESS, self.nTimeslot])
        U_essDisSol = np.zeros([self.nESS, self.nTimeslot])
        U_essChgSol = np.zeros([self.nESS, self.nTimeslot])
        if self.is_res_var:
            P_wtSol = np.zeros([self.nWT, self.nTimeslot])
            P_pvSol = np.zeros([self.nPV, self.nTimeslot])
        
        if self.is_uncertainty:
            P_wt_uncertaintySol = np.zeros([self.nWT, self.nTimeslot])
            P_pv_uncertaintySol = np.zeros([self.nPV, self.nTimeslot])
            smp_uncertaintySol = np.zeros([self.nTimeslot])
            
        for j in range(self.nTimeslot):
            
            try:
                P_BidSol[j] = self.m.getVarByName(f"Pbid[{j}]").X
            except:
                P_BidSol[j] = self.Pbid[j]
            
            for i in range(self.nESS):
                P_essDisSol[i,j] = self.m.getVarByName(f"P_essDis[{i},{j}]").X
                P_essChgSol[i,j] = self.m.getVarByName(f"P_essChg[{i},{j}]").X
                U_essDisSol[i,j] = self.m.getVarByName(f"U_essDis[{i},{j}]").X
                U_essChgSol[i,j] = self.m.getVarByName(f"U_essChg[{i},{j}]").X

            if self.is_res_var:    
                for i in range(self.nWT):
                    P_wtSol[i,j] = self.m.getVarByName(f"P_wt[{i},{j}]").X
                for i in range(self.nPV):
                    P_pvSol[i,j] = self.m.getVarByName(f"P_pv[{i},{j}]").X 
            
            if self.is_uncertainty:
                try:
                    smp_uncertaintySol[j] = self.m.getVarByName(f"smp_uncertainty[{j}]").X
                except:
                    abc = 3
                for i in range(self.nWT):
                    P_wt_uncertaintySol[i,j] = self.m.getVarByName(f"P_wt_uncertainty[{i},{j}]").X
                for i in range(self.nPV):
                    P_pv_uncertaintySol[i,j] = self.m.getVarByName(f"P_pv_uncertainty[{i},{j}]").X
                    
                        
                
        P_dict = {'bid': P_BidSol}
        U_dict = {}
        
        if self.is_res_var:
            P_dict['essDis'] = P_essDisSol
            P_dict['essChg'] = P_essChgSol
            U_dict['essDis'] = U_essDisSol
            U_dict['essChg'] = U_essChgSol
            
            P_dict['wt'] = P_wtSol
            P_dict['pv'] = P_pvSol
            
        if self.is_uncertainty:
            P_dict['wt_uncertainty'] = P_wt_uncertaintySol
            P_dict['pv_uncertainty'] = P_pv_uncertaintySol
            U_dict['smp_uncertainty'] = smp_uncertaintySol
            
        self.P_dict = P_dict
        self.U_dict = U_dict
        
        return P_dict, U_dict
    
    def optimize(self):
        
        self.add_Variables()
        self.add_bid_constraints()
        self.add_smp_constraints()
        self.add_res_constraints()
        self.add_ess_contraints()
        self.add_igdt_constraints()
        obj_eq = self.set_Objectives()
        
        mip_gap = 0.0001
        feas_tol = 1e-4
        sol, obj = self.solve([mip_gap, feas_tol])
        
        P_dict, U_dict = self.get_sol()
        
        return sol, obj, P_dict, U_dict

    
    def check_res_var_sol(self,var):

        print("**********")
        print("Check_RES_Constaraint_Solution")
        print("Check_RES_Constaraint_Solution")
        print("Check_RES_Constaraint_Solution")
        print("**********")

        P_wt = self.P_dict['wt'] 
        P_pv = self.P_dict['pv']
        P_essDis = self.P_dict['essDis']
        P_essChg = self.P_dict['essChg']
        
        self.check_set = {}
        if self.wt_list:
        
            wt_lb = np.zeros([self.nWT, self.nTimeslot])
            wt_ub = np.zeros([self.nWT, self.nTimeslot])
            
            for i in range(self.nWT):
                print("***********")
                
                if self.is_uncertainty:
                    print("P_wt[i,step], lb, ub, max_power, profile, uncert")
                else:
                    print("P_wt[i,step], lb, ub")
                    
                for step in range(self.nTimeslot):
                    lb, ub = self.check_ub_lb('WT', var, i, step, P_wt[i,step])
                    wt_lb[i,step] = lb
                    wt_ub[i,step] = ub
                                    
                print("***********")
            self.wt_bound_list = [P_wt,wt_lb,wt_ub]
            self.check_set['wt'] = self.wt_bound_list
            
        if self.pv_list:
        
            pv_lb = np.zeros([self.nPV, self.nTimeslot])
            pv_ub = np.zeros([self.nPV, self.nTimeslot])
            
            for i in range(self.nPV):
                print("***********")
                
                if self.is_uncertainty:
                    print("P_pv[i,step], lb, ub, max_power, profile, uncert")
                else:
                    print("P_pv[i,step], lb, ub")
                    

                for step in range(self.nTimeslot):
                    lb, ub = self.check_ub_lb('PV', var, i, step, P_pv[i,step])
                    pv_lb[i,step] = lb
                    pv_ub[i,step] = ub
            
                print("***********")
            self.pv_bound_list = [P_pv,pv_lb,pv_ub]
            self.check_set['pv'] = self.pv_bound_list
                    
        return self.check_set
    
    
    def check_ub_lb(self, res_type, var, i, step, component):

        if res_type == 'WT':
        
            min_power = self.wt_list[i].min_power
            max_power = self.wt_list[i].max_power
            uncert = self.wt_uncertainty[step]
            profile = self.wt_list[i].profile[step]
            Puncert = self.P_dict['wt_uncertainty'][i,step]
            
        elif res_type == 'PV':
            
            min_power = self.pv_list[i].min_power
            max_power = self.pv_list[i].max_power
            uncert = self.pv_uncertainty[step]
            profile = self.pv_list[i].profile[step]
            Puncert = self.P_dict['pv_uncertainty'][i,step]

        lb = min_power * (1-uncert)        
        
        if var == False and self.is_case_risk_averse:
            ub = max_power * (1-uncert)
        else:
            ub = max_power * (1+uncert)
        
        
        is_violate_lb = component < lb
        if is_violate_lb:            
            print("%%%%%%%%%%")
            print(f"violate_{res_type}_lb[{i},{step}]:", component, lb)
            print("%%%%%%%%%%")   
            
        is_violate_ub = component > ub
        if is_violate_ub:            
            print("%%%%%%%%%%")
            print(f"violate_{res_type}_ub[{i},{step}]:", component, ub)
            print("%%%%%%%%%%")   
        
        
        if self.is_uncertainty:
            print(f"var_P{res_type}_lb_ub[{i},{step}]_maxpower_profile_Puncertainty:", 
                  component, lb, ub, 
                  max_power, profile, Puncert)
            
        else:
            print(f"var_P{res_type}_lb_ub[{i},{step}]:", 
                  component, lb, ub)
        
        return lb, ub
    
    def check_bid_const_sol(self):
        
        self.lhs_bid = []
        self.rhs_bid = []
        
        P_bid = self.P_dict['bid']
        P_wt = self.P_dict['wt'] 
        P_pv = self.P_dict['pv']
        P_essDis = self.P_dict['essDis']
        P_essChg = self.P_dict['essChg']
        
        if self.is_res_var:
            
            print("check_bid_const_sol")
            print("check_bid_const_sol")
            print("check_bid_const_sol")
            
        
            for j in range(self.nTimeslot):
                
                lh = P_bid[j]
                rh = 0
                for i in range(self.nWT):
                    
                    rh += P_wt[i,j]
  
                for i  in range(self.nPV):
                    
                    rh += P_pv[i,j]

                for i in range(self.nESS):
                    
                    rh += P_essDis[i,j]

                    rh -= P_essChg[i,j]

                if lh > rh + 0.0001 or lh < rh - 0.0001:
                    
                    print(f"violation occurs at time {j}")
                    print(lh, rh)
                    
                self.lhs_bid.append(lh)
                self.rhs_bid.append(rh)
                
        const_list_set = list(zip(self.lhs_bid, self.rhs_bid))
        for i in range(len(const_list_set)):
            print(const_list_set[i])
                
        return self.lhs_bid, self.rhs_bid
    
    def check_res_const_sol(self):
        
        self.lhs = []
        self.rhs = []
        
        
        P_wt = self.P_dict['wt'] 
        P_pv = self.P_dict['pv']
        P_wt_uncertainty = self.P_dict['wt_uncertainty']
        P_pv_uncertainty = self.P_dict['pv_uncertainty']
        
        
        if self.is_res_var:
            
            print("Risk_averse_Check_res_const_sol")
            print("Risk_averse_Check_res_const_sol")
            print("Risk_averse_Check_res_const_sol")
            for j in range(self.nTimeslot):
                
                for i in range(self.nWT):
                    lh = P_wt[i,j]
                    
                    if self.is_case2:
                        
                        if self.is_case_risk_averse:
                            
                            rh = ((1-P_wt_uncertainty[i,j]) 
                                  * self.wt_list[i].max_power 
                                  * self.wt_list[i].profile[j])
                        
                        else:
                            rh = ((1+P_wt_uncertainty[i,j]) 
                            * self.wt_list[i].max_power 
                            * self.wt_list[i].profile[j])
                    else:
                        rh = self.wt_list[i].max_power * self.wt_list[i].profile[j]
                                     
                    if lh > rh + 0.001:
                        print(f"violation occurs on WT{i} at time {j} ")
                        print(lh, rh)                        
                        
                    self.lhs.append(lh)
                    self.rhs.append(rh)
                
                for i in range(self.nPV):
                    
                    lh = P_pv[i,j]
                    
                    if self.is_case2:
                        if self.is_case_risk_averse:
                            
                            rh = ((1-P_pv_uncertainty[i,j]) 
                            * self.pv_list[i].max_power 
                            * self.pv_list[i].profile[j])
                            
                        else:
                            
                            rh = ((1+P_wt_uncertainty[i,j]) 
                            * self.wt_list[i].max_power 
                            * self.wt_list[i].profile[j])
                    
                    else:
                        
                        rh = ((1+P_pv_uncertainty[i,j]) 
                        * self.pv_list[i].max_power 
                        * self.pv_list[i].profile[j])    
                        
                        
                    if lh > rh + 0.001:
                        print(f"violation occurs on PV{i} at time {j} ")
                        print(lh, rh)
                        
                            
     
                    self.lhs.append(lh)
                    self.rhs.append(rh)

                    
            bid_type = self.case_dict['bid_type']
            print(f'{bid_type}-check_res_const')
            print("(LHS, RHS) - WT + PV ")
            const_list_set = list(zip(self.lhs, self.rhs))
            for i in range(len(const_list_set)):
                print(const_list_set[i])
            
            return self.lhs, self.rhs
                