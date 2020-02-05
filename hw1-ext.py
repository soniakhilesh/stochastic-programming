from gurobipy import *
import matplotlib.pyplot as plt
import time
### Read data from file you choose: commont/uncomment to choose the different files

import hw1data1 as data
#import hw1data2 as data
#MODULAR APPROACH
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
        self.increase={}
        self.x={}
        self.z={}
        self.stoch_sol=0
        self.scenario_cost_stoch=[0]*len(Sset)
        self.scenario_cost_mean=[0]*len(Sset)
        
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
        self.stoch_sol = self.m.objVal
        print('Stochastic Solution Cost :{}'.format(stoch_sol))
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
        scenario_num = 0
        for k in Sset:
        	self.scenario_cost_mean[scenario_num] = sum(z[i,k].x*unmetCost[i] for i in Cset) + sum(increase[i,j].x*arcExpCost[i,j] for i,j in AllArcs)
        	scenario_num +=1
            
    def plot(self):
        plt.hist([self.scenario_cost_stoch,self.scenario_cost_mean], range=(0,18000), bins=18)
        plt.title('Stochastic vs Mean Value Solution')
        plt.xlabel('Cost')
        plt.ylabel('Number of Scenarios')
        plt.legend(['Stochastic Model solution','Mean value solution'])
        plt.savefig('hw1-hist.png')

#%%
#NON MODULAR APPROACG
Fset = data.Fset  # set of facilities (list of strings)
Hset = data.Hset  # set of warehouses (list of strings)
Cset = data.Cset  # set of customers (list of strings)
Sset = data.Sset  # set of scenarios (list of strings)
arcExpCost = data.arcExpCost  # arc expansion costs (dictionary mapping F,H and H,C pairs to floats)
facCap = data.facCap   # facility capacities (dictionary mapping F to floats) 
curArcCap = data.curArcCap  # current arc capacities (dictionary mapping (i,j) to floats, where either 
                                 # i is facility, j is warehouse, or i is warehouse and j is customer
unmetCost = data.unmetCost  # penalty for unment customer demands (dicationary mapping C to floats)
demScens = data.demScens  # demand scenarios (dictionary mapping (i,k) tuples to floats, where i is customer, k is
                                #scenario

### This is just a check of the data. You should probably comment out/delete these lines once you see the structure 

# =============================================================================
# print(Fset)
# print(Hset)
# print(Cset)
# print(Sset)
# print(arcExpCost)
# print(facCap)
# print(curArcCap)
# print(unmetCost)
# print(demScens)
# 
# =============================================================================

### Define sets of arcs (used as keys to dictionaries)
FHArcs = [(i,j) for i in Fset for j in Hset]  ## arcs from facilities to warehouses
HCArcs = [(i,j) for i in Hset for j in Cset]   ## arcs from warehouses to customers
AllArcs = FHArcs + HCArcs
# =============================================================================
# print(AllArcs)
# =============================================================================
unmetCostscaled ={}
for c in Cset:
    for k in Sset:
        unmetCostscaled[(c,k)]=unmetCost[c]/float(len(Sset))


##### Start building the Model #####
#create a model
start_time=time.time()
m=Model("JLWRC")

#create variables
#first stage
increase=m.addVars(AllArcs,obj=arcExpCost,name="increase") # increase on each arc
#second stage
x=m.addVars(AllArcs,Sset,name="flow") #flow on each arc
z=m.addVars(Cset,Sset,obj=unmetCostscaled, name="unmetdemand") #unmet demand
m.modelSense = GRB.MINIMIZE
#add constraints
# facility capacity constraints
m.addConstrs((x.sum(i,'*',k)<=facCap[i] for i in Fset for k in Sset),name="capacity")
#flow balance constraints at hub
m.addConstrs((x.sum('*',j,k)==x.sum(j,'*',k) for j in Hset for k in Sset),name="balance")
#fulfiling demandof each customer
m.addConstrs((x.sum('*',j,k)+z[j,k]>=demScens[j,k] for j in Cset for k in Sset),name="demand")
#arc capacity 
m.addConstrs((x[(*a),k]<=curArcCap[a]+increase[a]  for a in AllArcs for k in  Sset),name="arc-cap")
#bjective
m.update()
# =============================================================================
# m.write("file.lp")
# =============================================================================
m.optimize()
stoch_sol = m.objVal
print('Stochastic Solution Cost :{}'.format(stoch_sol))
print('\nExpansion in arcs')
for a in AllArcs:
    if increase[a].x >0.0001:
        print('Arc {}: {}'.format(a,increase[a].x))
print('\nCost of increasing capacity : {}'.format(sum(increase[a].x*arcExpCost[a] for a in AllArcs)))
print('\nDemand unmet (avg):') 
for c in Cset:
    unmet_avg = sum(z[c,k].x for k in Sset)/(len(Sset))
    print('Customer {}: {}'.format(c, unmet_avg))
print('\nExtensive model time:{}'.format(time.time() - start_time))

#calculate cost in each scenario
scenario_cost_stoch=[0]*len(Sset)
scenario_num= 0
for k in Sset:
	scenario_cost_stoch[scenario_num] = sum(z[c,k].x*unmetCost[c] for c in Cset) + sum(increase[a].x*arcExpCost[a] for a in AllArcs)
	scenario_num +=1

#now define mean value model


start_time_mean=time.time()
mean_m=Model("mean-val")

#create variables
#first stage
mv_increase=mean_m.addVars(AllArcs,obj=arcExpCost,name="increase") # increase on each arc
#second stage
mv_x=mean_m.addVars(AllArcs,name="flow") #flow on each arc
mv_z=mean_m.addVars(Cset,obj=unmetCost, name="unmetdemand") #unmet demand
mean_m.modelSense = GRB.MINIMIZE
#add constraints
# facility capacity constraints
mean_m.addConstrs((mv_x.sum(i,'*')<=facCap[i] for i in Fset),name="capacity")
#flow balance constraints at hub
mean_m.addConstrs((mv_x.sum('*',j)==mv_x.sum(j,'*') for j in Hset),name="balance")
#fulfiling demandof each customer
mean_m.addConstrs((mv_x.sum('*',j)+mv_z[j]>=sum(demScens[j,k] for k in Sset)/len(Sset) for j in Cset),name="demand")
#arc capacity 
mean_m.addConstrs((mv_x[a]<=curArcCap[a]+mv_increase[a]  for a in AllArcs),name="arc-cap")
#bjective
mean_m.update()
# =============================================================================
# mean_m.write("file.lp")
# =============================================================================
#We now plus in mean value solution in stoch model

mean_m.optimize()
mean_val_sol = mean_m.objVal
m.addConstrs(increase[a]==mv_increase[a].x for a in AllArcs)
m.update()
m.optimize()
print('\nMean Value Solution\n')
print('\nExpansion in arcs')
for a in AllArcs:
    if increase[a].x >0.0001:
        print('Arc {}: {}'.format(a,increase[a].x))
print('\nCost of increasing capacity in Mean Value Solution: {}'.format(sum(increase[a].x*arcExpCost[a] for a in AllArcs)))
print('\nDemand unmet (avg) in Mean Value Solution:') 
for c in Cset:
    unmet_avg = sum(z[c,k].x for k in Sset)/(len(Sset))
    print('Customer {}: {}'.format(c, unmet_avg))
print('\n Cost of Mean Value Solution:{}'.format(m.objVal))
print('VSS:{}'.format(m.objVal-stoch_sol))

#cost in each scenario for mean value solution
scenario_cost_mean = [0]*len(Sset)
scenario_num = 0
for k in Sset:
	scenario_cost_mean[scenario_num] = sum(z[i,k].x*unmetCost[i] for i in Cset) + sum(increase[i,j].x*arcExpCost[i,j] for i,j in AllArcs)
	scenario_num +=1

#creating histogram
plt.hist([scenario_cost_stoch,scenario_cost_mean], range=(0,18000), bins=18)
plt.title('Stochastic vs Mean Value Solution')
plt.xlabel('Cost')
plt.ylabel('Number of Scenarios')
plt.legend(['Stochastic Model solution','Mean value solution'])
plt.savefig('hw1-hist.png')
