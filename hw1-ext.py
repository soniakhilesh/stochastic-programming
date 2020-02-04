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
m=model("JLWRC")

#create variables
i=m.addVars(AllArcs,obj=arcExpCost,name="increase") # increase on each arc
x=m.addVars(AllArcs,Sset,name="flow") #flow on each arc
z=m.addVars(Cset,Sset,name="demand") #flow on each arc

#add constraints