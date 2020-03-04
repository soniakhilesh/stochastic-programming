from gurobipy import *
import time
import numpy as np

### Read data from file you choose: commont/uncomment to choose the different files

#import hw1data1 as data
import hw2data as data

class facility:
    def __init__(self):        
        self.F = data.F  # set of facilities (list of strings)
        self.D = data.D  # set of districts (list of strings)
        self.P = data.P  # set of item types (list of strings)
        self.S = data.S  # set of scenarios (list of strings)
        self.B = data.B
        self.fcap = data.fcap # capacity of each facility (dictionary mapping facility to number)
        self.a = data.a # space used per unit of item type (dictionary mapping item type to number)
        self.pcost = data.pcost # cost per unit or each item (dictionary mapping item type to number)
        self.closefac = data.closefac # set of facilities that are close enough for each district and item type (dictionary mapping (district,itemtype) to a list of facilities 
        self.facdistpairs = {} # For each product p, make a list of eligible district/facility pairs
        for p in self.P: 
            eligpairs = []
            for d in self.D:
                for f in self.closefac[d,p]:
                    eligpairs.append((f,d))
            self.facdistpairs[p] = tuplelist(eligpairs)
        ## define models
        self.m_ext=Model("Extensive")
        self.m_master=Model("Master")
        ##decision variables
        self.m_sub=Model("Sub")
        self.m_alt=Model("Alt")
        self.x={} #items p to place in fac f
        self.y={} #items p sent from f to d in scen k
        self.z={} #unmet demand of p in d in scen k
        self.gamma={} #penalty due to item p in dist d in scen k
        self.theta={} #aux variable to generate vars
        #alternative cut vars
        self.alpha2={}
        self.alpha_dem_1={}
        self.alpha_dem_2={}
        self.alpha_pen1_1={}
        self.alpha_pen1_2={}
        self.alpha_pen2_1={}
        self.alpha_pen2_2={}
        self.alpha_fopen_1={}
        self.alpha_fopen_2={}
        #storing second stage constraints
        self.demand={}
        self.penalty2={}
        self.fopen={}
        self.aux={}
        ##parameters
        self.demscens = data.demscens # requirements in each scenario (dictionary mappint (district,itemtype,scenario) to  number)
        self.availscens = data.availscens # availability of facilities in each scenario (dictionary mapping (facility,scenario) to 0,1 indicator of facility availability


    def model_extensive(self):
        F=self.F #faci
        D=self.D #dist
        P=self.P #item
        S=self.S #scen
        B=self.B #budget
        a=self.a #space
        fcap=self.fcap #fac cap
        demscens=self.demscens;
        closefac=self.closefac;
        availscens=self.availscens;
        pcost=self.pcost; #item cost
        #first stage
        self.x=self.m_ext.addVars(F,P,name='allocate');
        #second stage
        self.y=self.m_ext.addVars(F,P,D,S,name='flow');
        self.z=self.m_ext.addVars(D,P,S,name='unmet');
        self.gamma=self.m_ext.addVars(D,P,S,name='penalty');
        #define constraints
        #demand
        self.m_ext.addConstrs((self.z[d,p,k]+sum(self.y[f,p,d,k] for f in self.closefac[(d,p)])>=demscens[d,p,k] for d in D for p in P for k in S),name="demand")
        #budget
        self.m_ext.addConstr((sum(sum(self.x[f,p]*pcost[p] for f in F) for p in P)<=B),name="budget")
        #capacity
        self.m_ext.addConstrs((sum(a[p]*self.x[f,p] for p in P)<=fcap[f] for f in F),name="cap")    
        #penalty
        self.m_ext.addConstrs((self.gamma[d,p,k]>=self.z[d,p,k] for d in D for p in P for k in S),name="penalty1")        
        self.m_ext.addConstrs((self.gamma[d,p,k]>=0.5*demscens[d,p,k]+10*(self.z[d,p,k]-0.5*demscens[d,p,k]) for d in D for p in P for k in S),name="penalty2")        
        #open facility
        self.m_ext.addConstrs((self.y.sum(f,p,'*',k)<=self.x[f,p]*availscens[f,k] for p in P for f in F for k in S),name="open")
        #objective
        self.m_ext.setObjective((1/len(S))*self.gamma.sum('*','*','*'),GRB.MINIMIZE)
        self.m_ext.update()
        return self.m_ext
    
    def model_master(self):
        F=self.F #faci
        P=self.P #item
        S=self.S #scen
        B=self.B #budget
        a=self.a #space
        fcap=self.fcap #fac cap
        pcost=self.pcost; #item cost
        #decision vars
        self.x=self.m_master.addVars(F,P,name='allocate');
        self.theta=self.m_master.addVars(S,name='theta');
        #constraints
        #budget
        self.m_master.addConstr((sum(sum(self.x[f,p]*pcost[p] for f in F) for p in P)<=B),name="budget")
        #capacity
        self.m_master.addConstrs((sum(a[p]*self.x[f,p] for p in P)<=fcap[f] for f in F),name="cap")    
        #objective
        self.m_master.setObjective((1/len(S))*self.theta.sum('*'),GRB.MINIMIZE)
        self.m_master.setParam(GRB.Param.LogToConsole,0)
        self.m_master.update()
        #self.m_master.optimize() #solve master
        #return self.m_master

    def model_sub(self,k):
        #k is the scenario in set S
        F=self.F #faci
        D=self.D #dist
        P=self.P #item
        S=self.S #scen
        demscens=self.demscens;
        closefac=self.closefac;
        availscens=self.availscens;
        #decision vars without k index as defined for each scenario
        self.y=self.m_sub.addVars(F,P,D,name='flow');
        self.z=self.m_sub.addVars(D,P,name='unmet');
        self.gamma=self.m_sub.addVars(D,P,name='penalty');
        #constraints
        #demand
        self.demand=self.m_sub.addConstrs((self.z[d,p]+sum(self.y[f,p,d] for f in closefac[(d,p)])>=demscens[d,p,k] for d in D for p in P),name="demand")
        #penalty
        self.m_sub.addConstrs((self.gamma[d,p]>=self.z[d,p] for d in D for p in P),name="penalty1")        
        self.penalty2=self.m_sub.addConstrs((self.gamma[d,p]-10*self.z[d,p]>=-4.5*demscens[d,p,k] for d in D for p in P),name="penalty2")        
        #open facility
        self.fopen=self.m_sub.addConstrs((self.y.sum(f,p,'*')<=self.x[f,p].x*availscens[f,k] for f in F for p in P),name="open")
        #objective
        self.m_sub.setObjective(self.gamma.sum('*','*'),GRB.MINIMIZE)
        self.m_sub.setParam(GRB.Param.LogToConsole,0)

        self.m_sub.update()
        return self.demand,self.penalty2,self.fopen

    def benders_multi_cut(self):
        cutfound=True;
        F=self.F #faci
        D=self.D #dist
        P=self.P #item
        S=self.S #scen
        availscens=self.availscens;
        demscens=self.demscens;
        #solve master
        start = time.time()
        self.model_master() #build master model
        self.m_master.optimize() #solve master to build sub
        #build sub model for some k and then update necessary part of fmodel
        dem,pen,fo=self.model_sub(S[0])
        itera=0
        while cutfound==True:
            #update and re-solve master model
            self.m_master.update()
            print('Re solving Master', itera)
            self.m_master.optimize()
            print('lower bound = ', self.m_master.objval)
            itera+=1
            cutfound=False
            sscost = 0.0 
            for k in S:
                #change rhs of sub problem constraints
                for p in P:
                    for d in D:
                        self.demand[d,p].RHS=demscens[d,p,k]
                        self.penalty2[d,p].RHS=-4.5*demscens[d,p,k]
                    for f in F:
                        self.fopen[f,p].RHS=self.x[f,p].x*availscens[f,k]
                #resolve sub model
                self.m_sub.update()
                self.m_sub.optimize()
                obj=self.m_sub.objVal;
                dem=self.demand
                pen=self.penalty2
                fo=self.fopen
                rhs=0;
                if obj>self.theta[k].x+0.000001:
                    sscost += obj
                    cutfound=True #cut found
                    #calculate cut which is obj of dual
                    rhs+=sum(dem[d,p].pi*demscens[d,p,k] for d in D for p in P)
                    rhs+=sum(-4.5*pen[d,p].pi*demscens[d,p,k] for d in D for p in P)
                    xcoef={} # define xcoeff
                    for f in F:
                        for p in P:
                            xcoef[f,p]=availscens[f,k]*fo[f,p].pi;
                    rhs+=self.x.prod(xcoef)
                    #addcut
                    self.m_master.addConstr(self.theta[k] >= rhs)
            curub = (1/float(len(S)))*sscost #there is no first stage cost
            print('upper bound = ', curub)
        
        print('total time = ', time.time() - start)
        print('total iterations = ', itera)
        print('OPtimal Solution = ',self.m_master.objval)
            
    def alternative_sub(self,k):
        #k is the scenario in set S
        F=self.F #faci
        D=self.D #dist
        P=self.P #item
        demscens=self.demscens;
        closefac=self.closefac;
        availscens=self.availscens;
        
        #decision vars without k index as defined for each scenario
        self.y=self.m_alt.addVars(F,P,D,name='flow');
        self.z=self.m_alt.addVars(D,P,name='unmet');
        self.gamma=self.m_alt.addVars(D,P,name='penalty');
        
        #new vars
        self.alpha2=self.m_alt.addVar(lb=0,name='alpha-aux')
        #2 var for dem con
        self.alpha_dem_1=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-dem-1')
        self.alpha_dem_2=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-dem-2')
        #2 vars for pen1 con
        self.alpha_pen1_1=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-pen1-1')
        self.alpha_pen1_2=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-pen1-2')
        #2 var for pen2 con
        self.alpha_pen2_1=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-pen2-1')
        self.alpha_pen2_2=self.m_alt.addVars(D,P,lb=-GRB.INFINITY,ub=0,name='alpha-pen2-2')
        #2 vars for pen con
        self.alpha_fopen_1=self.m_alt.addVars(F,P,lb=-GRB.INFINITY,ub=0,name='alpha-open-1')
        self.alpha_fopen_2=self.m_alt.addVars(F,P,lb=-GRB.INFINITY,ub=0,name='alpha-open-2')
        
        #objective
        self.m_alt.setObjective(self.alpha2,GRB.MINIMIZE)
        self.m_alt.setParam(GRB.Param.LogToConsole,0)

       
        #constraints
        #aux constraint
        self.aux=self.m_alt.addConstr((self.alpha2-self.gamma.sum('*','*')>=-self.theta[k].x),name="aux")
        #demand
        self.demand=self.m_alt.addConstrs((self.z[d,p]+sum(self.y[f,p,d] for f in closefac[(d,p)])-self.alpha_dem_1[d,p]+self.alpha_dem_2[d,p]>=demscens[d,p,k] for d in D for p in P),name="demand")
        #penalty
        self.m_alt.addConstrs((self.gamma[d,p]-self.alpha_pen1_1[d,p]+self.alpha_pen1_2[d,p]>=self.z[d,p] for d in D for p in P),name="penalty1")        
        self.penalty2=self.m_alt.addConstrs((self.gamma[d,p]-10*self.z[d,p]-self.alpha_pen2_1[d,p]+self.alpha_pen2_2[d,p]>=-4.5*demscens[d,p,k] for d in D for p in P),name="penalty2")        
        #open facility
        self.fopen=self.m_alt.addConstrs((self.y.sum(f,p,'*')-self.alpha_fopen_1[f,p]+self.alpha_fopen_2[f,p]<=self.x[f,p].x*availscens[f,k] for f in F for p in P),name="open")
        #add new constraints
        self.m_alt.addConstrs((self.alpha2+self.alpha_dem_1[d,p]+self.alpha_dem_2[d,p]==0 for d in D for p in P),name="new-aux1")        
        self.m_alt.addConstrs((self.alpha2+self.alpha_pen1_1[d,p]+self.alpha_pen1_2[d,p]==0 for d in D for p in P),name="new-aux2")        
        self.m_alt.addConstrs((self.alpha2+self.alpha_pen2_1[d,p]+self.alpha_pen2_2[d,p]==0 for d in D for p in P),name="new-aux3")
        self.m_alt.addConstrs((self.alpha2+self.alpha_fopen_1[f,p]+self.alpha_fopen_2[f,p]==0 for f in F for p in P),name="new-aux4")
        #update the model
        self.m_alt.update()
        return self.demand,self.penalty2,self.fopen
        
    def benders_alt_cut(self):
        cutfound=True;
        F=self.F #faci
        D=self.D #dist
        P=self.P #item
        S=self.S #scen
        availscens=self.availscens;
        demscens=self.demscens;
        #solve master
        start = time.time()
        self.model_master() #build master model
        self.m_master.optimize() #solve master to build sub
        #build sub model for some k and then update necessary part of fmodel
        dem,pen,fo=self.alternative_sub(S[0])
        itera=0
        while cutfound==True:
            #update and re-solve master model
            self.m_master.update()
            print('Re solving Master', itera)
            self.m_master.optimize()
            print('lower bound = ', self.m_master.objval)
            itera+=1
            cutfound=False
            sscost = 0.0 
            for k in S:
                #change rhs of sub problem constraints
                self.aux.RHS=-self.theta[k].x;
                for p in P:
                    for d in D:
                        self.demand[d,p].RHS=demscens[d,p,k]
                        self.penalty2[d,p].RHS=-4.5*demscens[d,p,k]
                    for f in F:
                        self.fopen[f,p].RHS=self.x[f,p].x*availscens[f,k]
                #resolve sub model
                self.m_alt.update()
                self.m_alt.optimize()
                #if model is feasible
                if self.m_alt.status==2: #optimal solution exists
                    #add optimality cuts
                    obj=self.m_alt.objVal;
                    dem=self.demand
                    pen=self.penalty2
                    fo=self.fopen
                    rhs=0;
                    if obj>0.000001:
                        #sscost += obj+self.theta[k].x
                        cutfound=True #cut found
                        #calculate cut which is obj of dual
                        rhs+=sum(dem[d,p].pi*demscens[d,p,k]/self.aux.pi for d in D for p in P)
                        rhs+=sum(-4.5*pen[d,p].pi*demscens[d,p,k]/self.aux.pi for d in D for p in P)
                        xcoef={} # define xcoeff
                        for f in F:
                            for p in P:
                                xcoef[f,p]=availscens[f,k]*fo[f,p].pi/self.aux.pi;
                        rhs+=self.x.prod(xcoef)
                        #addcut
                        self.m_master.addConstr(self.theta[k] >= rhs)
                    
            #curub = (1/float(len(S)))*sscost #there is no first stage cost
            #print('upper bound = ', curub)
        
        print('total time = ', time.time() - start)
        print('total iterations = ', itera)
        print('OPtimal Solution = ',self.m_master.objval)
      
        

        
        

        
        
