from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt

### Read data from file you choose: commont/uncomment to choose the different files

import hw1data1 as data
#import hw1data2 as data

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


##### Start building the Model #####
#create a model
m=Model("JLWRC")

#create variables
i=m.addVars(AllArcs,obj=arcExpCost,name="increase",lb=0) # increase on each arc
x=m.addVars(AllArcs,Sset,name="flow",lb=0) #flow on each arc
z=m.addVars(Cset,Sset,obj=unmetCost*(1/len(Sset)),name="demand",lb=0) #unmet demand

#add constraints
# facility capacity constraints
m.addConstrs((x.sum(i,'*',k)<=facCap[i] for i in Fset for k in Sset),"capacity")
#flow balance constraints at hub
m.addConstrs((x.sum('*',j,k)==x.sum(j,'*',k) for j in Hset for k in Sset),"balance")
#fulfiling demandof each customer
m.addConstrs((x.sum('*',j,k)+z[j,k]==demScens(j,k) for j in Cset for k in Sset),"demand")
#arc capacity 
m.addConstrs((x[i,j,k]<=curArcCap[i,j]+i[i,j]  for i,j in AllArcs for k in  Sset),"arc-cap")
m.setObjective(GRB.MAXIMIZE)
m.update()
m.optimize()
