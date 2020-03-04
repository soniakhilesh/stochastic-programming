#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 21:07:37 2020

@author: soni6
"""

from gurobipy import *

prob=[0.25,0.5,0.25]
theta_list=[1,2,3]

m=Model('Master')
x=m.addVar(lb=0,ub=20,name="x")
theta=m.addVars(theta_list,lb=0,name="theta")

m.setObjective(1.5*x+sum(theta[i]*prob[i-1] for i in theta_list))
m.modelSense = GRB.MINIMIZE
m.addConstr(theta[1]>=4-4*x,name="cut1")
m.addConstr(theta[2]>=16-4*x,name="cut2")
m.addConstr(theta[3]>=32-4*x,name="cut3")
#iter 1 cuts
m.addConstr(theta[1]>=4*x-10,name="cut4")
#m.addConstr(theta[2]>=4*x-22,name="cut5")
#iter 2 cuts
#m.addConstr(theta[2]>=4*x-38,name="cut6")

m.optimize()
print('x',x.x)
print('theta',theta[1].x,theta[2].x,theta[3].x)
print('obj val',m.objVal)