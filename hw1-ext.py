from gurobipy import *
import matplotlib.pyplot as plt
import time
import copy
### Read data from file you choose: commont/uncomment to choose the different files

#import hw1data2 as data
import hw1data2 as data
class jeffrey_model:
    def __init__(self):
        self.Fset = data.Fset  # set of facilities (list of strings)
        self.Hset = data.Hset  # set of warehouses (list of strings)
        self.Cset = data.Cset  # set of customers (list of strings)
        self.Sset = data.Sset  # set of scenarios (list of strings)
        self.arcExpCost = data.arcExpCost  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
        self.facCap = data.facCap   # facility capacities (dictionary mapping F to floats) 
        self.curArcCap = data.curArcCap  # current arc capacities (dictionary mapping (i,j) to floats, where either 
                                         # i is facility, j is warehouse, or i is warehouse and j is customer
        self.unmetCost = data.unmetCost  # penalty for unment customer demands (dicationary mapping C to floats)
        self.demScens = data.demScens  # demand scenarios (dictionary mapping (i,k) tuples to floats, where i is customer, k is
                                        #scenario
        self.FHArcs = [(i,j) for i in self.Fset for j in self.Hset]  ## arcs from facilities to warehouses
        self.HCArcs = [(i,j) for i in self.Hset for j in self.Cset]   ## arcs from warehouses to customers
        self.AllArcs = self.FHArcs + self.HCArcs
        self.m=Model("Extensive")
        self.mean_m=Model("Mean-Val")
        self.avar_m=Model("Mean-Avar")
        self.increase={}
        self.x={}
        self.z={}
        self.extra={}
        self.gamma={}
        self.scen_cost={}
        self.stoch_sol=0
        self.scenario_cost_stoch=[0]*len(self.Sset)
        self.scenario_cost_mean=[0]*len(self.Sset)
        self.scenario_cost_avar=[0]*len(self.Sset)
        
    def solve_extensive(self):
        Fset=self.Fset
        Hset=self.Hset
        Cset=self.Cset
        Sset=self.Sset
        arcExpCost=self.arcExpCost
        facCap=self.facCap
        curArcCap=self.curArcCap
        unmetCost=self.unmetCost
        demScens=self.demScens
        AllArcs=self.AllArcs
        unmetCostscaled ={}
        for c in Cset:
            for k in Sset:
                unmetCostscaled[(c,k)]=unmetCost[c]/float(len(Sset))
        start_time=time.time()
        #first stage vars
        self.increase=self.m.addVars(AllArcs,obj=arcExpCost,name="increase") # increase on each arc
        #second stage
        self.x=self.m.addVars(AllArcs,Sset,name="flow") #flow on each arc
        self.z=self.m.addVars(Cset,Sset,obj=unmetCostscaled, name="unmetdemand") #unmet demand
        self.m.modelSense = GRB.MINIMIZE
        #add constraints
        # facility capacity constraints
        self.m.addConstrs((self.x.sum(i,'*',k)<=facCap[i] for i in Fset for k in Sset),name="capacity")
        #flow balance constraints at hub
        self.m.addConstrs((self.x.sum('*',j,k)==self.x.sum(j,'*',k) for j in Hset for k in Sset),name="balance")
        #fulfiling demandof each customer
        self.m.addConstrs((self.x.sum('*',j,k)+self.z[j,k]>=demScens[j,k] for j in Cset for k in Sset),name="demand")
        #arc capacity 
        self.m.addConstrs((self.x[(*a),k]<=curArcCap[a]+self.increase[a]  for a in AllArcs for k in  Sset),name="arc-cap")
        #bjective
        self.m.update()
        self.m.optimize()
        self.m.write("ext.lp")
        self.stoch_sol = self.m.objVal
        print('Stochastic Solution Cost :{}'.format(self.stoch_sol))
        print('\nExpansion in arcs')
        for a in AllArcs:
            if self.increase[a].x >0.0001:
                print('Arc {}: {}'.format(a,self.increase[a].x))
        print('\nCost of increasing capacity : {}'.format(sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)))
        print('\nDemand unmet (avg):') 
        for c in Cset:
            unmet_avg = sum(self.z[c,k].x for k in Sset)/(len(Sset))
            print('Customer {}: {}'.format(c, unmet_avg))
        print('\nExtensive model time:{}'.format(time.time() - start_time))
        
        #calculate cost in each scenario

        scenario_num= 0
        for k in Sset:
            self.scenario_cost_stoch[scenario_num] = sum(self.z[c,k].x*unmetCost[c] for c in Cset) + sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)
            scenario_num +=1
    
    def solve_mean_val(self):
        Fset=self.Fset
        Hset=self.Hset
        Cset=self.Cset
        Sset=self.Sset
        arcExpCost=self.arcExpCost
        facCap=self.facCap
        curArcCap=self.curArcCap
        unmetCost=self.unmetCost
        demScens=self.demScens
        AllArcs=self.AllArcs
        start_time=time.time()
        #create variables
        #first stage
        self.mv_increase=self.mean_m.addVars(AllArcs,obj=arcExpCost,name="increase") # increase on each arc
        #second stage
        self.mv_x=self.mean_m.addVars(AllArcs,name="flow") #flow on each arc
        self.mv_z=self.mean_m.addVars(Cset,obj=unmetCost, name="unmetdemand") #unmet demand
        self.mean_m.modelSense = GRB.MINIMIZE
        #add constraints
        # facility capacity constraints
        self.mean_m.addConstrs((self.mv_x.sum(i,'*')<=facCap[i] for i in Fset),name="capacity")
        #flow balance constraints at hub
        self.mean_m.addConstrs((self.mv_x.sum('*',j)==self.mv_x.sum(j,'*') for j in Hset),name="balance")
        #fulfiling demandof each customer
        self.mean_m.addConstrs((self.mv_x.sum('*',j)+self.mv_z[j]>=sum(demScens[j,k] for k in Sset)/len(Sset) for j in Cset),name="demand")
        #arc capacity 
        self.mean_m.addConstrs((self.mv_x[a]<=curArcCap[a]+self.mv_increase[a]  for a in AllArcs),name="arc-cap")
        #bjective
        self.mean_m.update()
        # =============================================================================
        # mean_m.write("file.lp")
        # =============================================================================
        #We now plus in mean value solution in stoch model
        
        self.mean_m.optimize()
        mean_val_sol = self.mean_m.objVal
        self.m.addConstrs(self.increase[a]==self.mv_increase[a].x for a in AllArcs)
        self.m.update()
        self.m.optimize()
        print('\nMean Value Solution\n')
        print('\nExpansion in arcs')
        for a in AllArcs:
            if self.increase[a].x >0.0001:
                print('Arc {}: {}'.format(a,self.increase[a].x))
        print('\nCost of increasing capacity in Mean Value Solution: {}'.format(sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)))
        print('\nDemand unmet (avg) in Mean Value Solution:') 
        for c in Cset:
            unmet_avg = sum(self.z[c,k].x for k in Sset)/(len(Sset))
            print('Customer {}: {}'.format(c, unmet_avg))
        print('\n Cost of Mean Value Solution:{}'.format(self.m.objVal))
        print('VSS:{}'.format(self.m.objVal-self.stoch_sol))
        
        #cost in each scenario for mean value solution
        scenario_num= 0
        for k in Sset:
            self.scenario_cost_mean[scenario_num] = sum(self.z[c,k].x*unmetCost[c] for c in Cset) + sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)
            scenario_num +=1

            
    def plot_ext_mean(self):
        plt.hist([self.scenario_cost_stoch,self.scenario_cost_mean], range=(0,18000), bins=18)
        plt.title('Stochastic vs Mean Value Solution')
        plt.xlabel('Cost')
        plt.ylabel('Number of Scenarios')
        plt.legend(['Stochastic Model solution','Mean value solution'])
        plt.savefig('hw1-hist.png')
        
    def solve_mean_avar(self,lamb,alpha):
        Fset=self.Fset
        Hset=self.Hset
        Cset=self.Cset
        Sset=self.Sset
        arcExpCost=self.arcExpCost
        facCap=self.facCap
        curArcCap=self.curArcCap
        unmetCost=self.unmetCost
        demScens=self.demScens
        AllArcs=self.AllArcs


        start_time=time.time()
        #first stage vars
        self.increase=self.avar_m.addVars(AllArcs,name="increase") # increase on each arc
        self.gamma=self.avar_m.addVar(lb=-GRB.INFINITY,name="gamma") # avar variable
        #second stage
        self.x=self.avar_m.addVars(AllArcs,Sset,name="flow") #flow on each arc
        self.z=self.avar_m.addVars(Cset,Sset, name="unmetdemand") #unmet demand
        self.extra=self.avar_m.addVars(Sset, name="extra") #avar variable
        self.scen_cost=self.avar_m.addVars(Sset, name="scen-cost") #storing mean objective cost
        self.avar_m.modelSense = GRB.MINIMIZE
        #add constraints
        # facility capacity constraints
        self.avar_m.addConstrs((self.x.sum(i,'*',k)<=facCap[i] for i in Fset for k in Sset),name="capacity")
        #flow balance constraints at hub
        self.avar_m.addConstrs((self.x.sum('*',j,k)==self.x.sum(j,'*',k) for j in Hset for k in Sset),name="balance")
        #fulfiling demandof each customer
        self.avar_m.addConstrs((self.x.sum('*',j,k)+self.z[j,k]>=demScens[j,k] for j in Cset for k in Sset),name="demand")
        #arc capacity 
        self.avar_m.addConstrs((self.x[(*a),k]<=curArcCap[a]+self.increase[a]  for a in AllArcs for k in  Sset),name="arc-cap")
        #AVAR auxillary constraints
        self.avar_m.addConstrs((
                self.scen_cost[k]==sum(self.increase[a]*arcExpCost[a]for a in AllArcs)+
                sum(unmetCost[j]*self.z[j,k] for j in Cset) for k in Sset),name="scen_cost-con" )
        self.avar_m.addConstrs(self.extra[k]>=self.scen_cost[k]-self.gamma for k in Sset)
        #objective
        self.avar_m.setObjective(float(1/len(Sset))*sum(self.scen_cost[k] for k in Sset)+
                                 lamb*(self.gamma+float(1/((1-alpha)*len(Sset)))*sum(self.extra[k] for k in Sset)))
        self.avar_m.update()
        self.avar_m.optimize()
        self.avar_m.write("avar.lp")
# =============================================================================
#         self.mean_avar_sol = self.avar_m.objVal
#         print('Mean AVAR Solution Cost :{}'.format(self.mean_avar_sol))
# =============================================================================
        print('\nExpected Cost: {}'.format(sum(self.scen_cost[k].x for k in Sset)/len(Sset)))
        print('\nAV@R of COST : {}'.format(self.gamma.x + sum(self.extra[k].x for k in Sset)/(len(Sset)*(1-alpha))))

# =============================================================================
#         print('\nExpansion in arcs in Mean AVAR')
#         for a in AllArcs:
#             if self.increase[a].x >0.0001:
#                 print('Arc {}: {}'.format(a,self.increase[a].x))
#         print('\nCost of increasing capacity in mean-AVAR: {}'.format(sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)))
#         print('\nDemand unmet (avg) in mean AVAR:') 
#         for c in Cset:
#             unmet_avg = sum(self.z[c,k].x for k in Sset)/(len(Sset))
#             print('Customer {}: {}'.format(c, unmet_avg))
# =============================================================================
        print('\nAvar model time:{}'.format(time.time() - start_time))
        
        #calculate cost in each scenario

        scenario_num= 0
        for k in Sset:
            self.scenario_cost_avar[scenario_num] = sum(self.z[c,k].x*unmetCost[c] for c in Cset) + sum(self.increase[a].x*arcExpCost[a] for a in AllArcs)
            scenario_num +=1
        return self.scenario_cost_avar
            
    def plot_avar_mean(self):
        alpha=0.9
        lambda_0=copy.deepcopy(self.solve_mean_avar(0,alpha))
        lambda_1=copy.deepcopy(self.solve_mean_avar(10,alpha))
        lambda_100=copy.deepcopy(self.solve_mean_avar(1000,alpha))
        plt.hist([lambda_0,lambda_1,lambda_100], range=(0,18000), bins=18)
        plt.title('AV@R')
        plt.xlabel('Cost')
        plt.ylabel('Number of Scenarios')
        plt.legend(['Lambda=0','Lambda=1','Lambda=100'])
        plt.savefig('hw1-hist.png')

#%%
