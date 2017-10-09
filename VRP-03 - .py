# -*- coding: utf-8 -*-
"""
This program offers a heuristic solution to the Vehicle Routng Problem problem. The algorithm is based on two steps.
The first step uses an improved adaptive Simulated Annealing meta-heuristic to minimize the number of routes. 
The basic SA is modified to adaptively adjust its search span based on the quality of previous solutions obtained, 
i.e. if significant improvements do not occur repetitively, search span and depth are increased (and vice versa).
This adaptive feature significant improves the speed of algorithm.
The second part, inport the solution found via SA and minimizes the total transportation cost via a Large Neighborhood Search heuristic. 
The algorith is designed to handle very large distribution networks with multiple stores and dates. 
The algorithm also offers a simple solution to the shortcoming of many VRP solutions in optimizing fleet size. 
By adding a mechanism to assign routes to trucks the algorithm is capable of optimizing the fleet size as part of its objective function. 
"""


"""
VRP Algorithm: part 1 - Simulated Annealing with Adaptive Neighborhood Search
"""

import csv

import os, re
import pandas, numpy
from math import *
import datetime
import time
#from HTMLParser import HTMLParser

import gmplot
import geocoder
import random

import matplotlib.pyplot as plt
import networkx as nx
import copy, functools
# <codecell>
global customersplus, file_address, costs_achieved, costs_achievedII, customers, store_id, arcs,date, drive_time, dist_mtx,tg1removal,max_no_improvement
# <codecell>
tg1removal=0 #indicates whether tg1 is removed
temperature=100
beta=2   #neighbor selection performance bias
max_iterations=50
alpha=0.95
NSize=5 # Neighborhood Size
### Parameters
store_id=0
date=''
#date='2015-01-02'
service_time=0.25 #15minutes
max_stops=10
wage=5
driving_cost_per_km=0.25
max_route_duration=12
#
p=2
maxiter1=2
maxiter2=5
max_no_improvement=10

try: os.chdir("C:\\...\\VRP Analysis")
except: os.chdir("C:\\...\\VRP Analysis")

# <codecell>

def dist(lat1,lon1,lat2,lon2):
    if (lat1==lat2) and (lon2==lon1): return 0
    else:  return acos(cos(radians(90-lat1))*cos(radians(90-lat2))+sin(radians(90-lat1))*sin(radians(90-lat2)) *cos(radians(lon1-lon2)))*6371*1.3  #multiplied by 1.3 to account for road network features

def network(adjacency_matrix):
    mylabels=dict(zip(range(len(adjacency_matrix)),adjacency_matrix.columns)) #a dictionary of labels
    rows, cols = numpy.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.DiGraph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500, labels=mylabels, with_labels=True)
    plt.show()



# <codecell>
def successor(arcs_mtx,i):
    succ=arcs_mtx.columns[(arcs_mtx.loc[i,:]==1)]
    if len(succ)==0:
        return ''
    else: return succ[0]
    
def predecessor(arcs_mtx,i):
    pred=arcs_mtx.columns[(arcs_mtx.loc[:,i]==1)]
    if len(pred)==0: 
        print 'no predecessor!'
        return ''
    else: return pred [0]
 
# <codecell>
#representation of a solution
#zeros_matrix = numpy.array([numpy.zeros(row_count+1)]*(row_count+1))

def solution(arcs_mtx,printroute=False):
    global temperature, beta, max_iterations, alpha, NSize, service_time, max_stops, wage,driving_cost_per_km,max_route_duration, customersplus, customers,store_id, drive_time, dist_mtx
    arcs_mtx1=arcs_mtx.copy()
    #### A. Initialize solution matrix
    sol=pandas.DataFrame(index=customersplus,columns=['OrderID','Start_Time','End_Time','route','stop_number','arrival_time','drive_time','drive_distance','delivery_start','infeasible_flag','route_start_time','route_end_time','return_distance', 'route_duration'])
    #fill some basic information: order id and TW
    for ii in customers:
        for col in ['OrderID','Start_Time','End_Time']:
            sol.loc[ii,col]= orders[orders['Customer_ID2']==ii][col].iloc[0]

    sol['infeasible_flag']=sol['route']=0
    sol.loc[store_id,'Start_Time']=sol.loc[store_id,'arrival_time']= 7
    sol.loc[store_id,'End_Time']= 23
    
    ### B. follow routes, create the solution, and check feasibility
    for route_number in range(1,row_count+1):
        # find route head
        for i in customers:
            if (arcs_mtx1.loc[store_id,i]==1) and (sol.loc[i,'route']==0):       #skip the ones that are already visited or are not linked to depo
                break
            elif i==customers[row_count-1]:
                if printroute: print '\nComplete!'
                return sol      
        
        if printroute: print '\nroute', route_number, ': ',
        
        prev=store_id
        stop_number=1
    
        #loop throught the route & find successors    
        for stop_number in range(1,row_count+1):
            if printroute: print i,
            sol.loc[i,'route']=route_number
            
            #calculate variables
            sol.loc[i,'drive_time']= drive_time.loc[prev,i] #just for debugging; delete later
            sol.loc[i,'drive_distance']= dist_mtx.loc[prev,i] #just for debugging; delete later
            
            if stop_number==1:
                sol.loc[i,'route_start_time']= sol.loc[i,'End_Time']-drive_time.loc[prev,i]-service_time  #depart as late as possible
                sol.loc[i,'arrival_time']= sol.loc[i,'route_start_time']+ drive_time.loc[prev,i]
                route_start_time=sol.loc[i,'route_start_time'] #save start time
            else:
                sol.loc[i,'arrival_time']= sol.loc[prev,'arrival_time']+ drive_time.loc[prev,i] + service_time 
            
            #time constraints
            delivery_start=max(sol.loc[i,'arrival_time'],sol.loc[i,'Start_Time'])
            sol.loc[i,'delivery_start']=delivery_start
            if delivery_start>sol.loc[i,'End_Time']:  sol.loc[i,'infeasible_flag']=1 #check feasibility according to time windows
        
            #load constraints
            sol.loc[i,'stop_number']=stop_number
            if stop_number>max_stops: sol.loc[i,'infeasible_flag']=1 #check feasibility of maximum stops on a route
            
            
            prev=i
            i=successor(arcs_mtx1,i)
            #find route end
            if i==store_id: 
                sol.loc[prev,'route_end_time']=sol.loc[prev,'arrival_time']+ drive_time.loc[prev,i] + service_time
                sol.loc[prev,'return_distance']=dist_mtx.loc[prev,i]
                sol.loc[prev,'route_duration']=sol.loc[prev,'route_end_time']-route_start_time
                if sol.loc[prev,'route_duration']>max_route_duration:  sol.loc[prev,'infeasible_flag']=1 #check feasibility of route duration
                break
    return sol


    
# <codecell>                    
# Cost Function for the SA part. Emphasizes route reduction.
def cost1(solution_mtx):
    #solution_mtx=solution_mtx.copy()
    NOrders=len(solution_mtx)-1
    NRoutes=solution_mtx.route.max()
    if NOrders/NRoutes<4: NTrucks=100
    else: NTrucks=TruckCounter(solution_mtx)
     #for speed improvement do not deal with Truck number in the beginning
    ### calculate the cost factor for distribution of stops
    
    stops=list()
    for r in range(1,NRoutes+1):
        stops_on_route=solution_mtx[solution_mtx.route==r].stop_number.max()
        stops.append(stops_on_route)
    stops_cost_factor = sum([i**2 for i in stops])
        
    total_distance=solution_mtx['drive_distance'].sum(skipna=True) + solution_mtx['return_distance'].sum(skipna=True)
    total_time=solution_mtx['route_duration'].sum(skipna=True)
    shipping_cost=total_distance*driving_cost_per_km + total_time*wage + NTrucks*100
    service_cost=NOrders*wage*service_time
    
    total_cost=NTrucks*10**10+NRoutes*10**8+(10**7-stops_cost_factor*10**4)+shipping_cost +service_cost  #give a large weight to Nroutes so that route number is minimized first
    #print '\nTotal Distance:',total_distance,'\nTotal Time:',total_time,'\nTotal Cost:',total_cost
    return total_cost

def cost2(solution_mtx):
    global driving_cost_per_km, wage
    #solution_mtx=solution_mtx.copy()
    NRoutes=solution_mtx.route.max()
    NTrucks=TruckCounter(solution_mtx)
    NOrders=len(solution_mtx)-1
    
    total_distance=solution_mtx['drive_distance'].sum(skipna=True) + solution_mtx['return_distance'].sum(skipna=True)
    total_time=solution_mtx['route_duration'].sum(skipna=True)
    shipping_cost=total_distance*driving_cost_per_km + total_time*wage + NTrucks*100
    service_cost=NOrders*wage*service_time
 
    total_cost=NTrucks*10**8+NRoutes*10**6+shipping_cost + service_cost  #give a large weight to Nroutes so that route number is minimized first
    #print '\nTotal Distance:',total_distance,'\nTotal Time:',total_time,'\nTotal Cost:',total_cost
    return int(total_cost)

def cost2_from_arcs(arcs_mtx):
    return cost2(solution(arcs_mtx))
           
# <codecell>                    
def feasible(solution_mtx):
    f=solution_mtx['infeasible_flag'].sum(skipna=True)
    if f==0: return 1
    elif f>0: return 0
    else: return 'Error'

# <codecell>                    
# Function: Relocate i after j
def relocate(arcs_mtx,i,j):
    arcs_new=arcs_mtx.copy()
    # Take i from where it was
    pred_i=predecessor(arcs_mtx,i)
    succ_i=successor(arcs_mtx,i)
    if pred_i!=succ_i: arcs_new.loc[pred_i,succ_i]=1             # create an arc pred(i) -> succ(i) ; only if both are not the depo
    arcs_new.loc[i,succ_i]=0                              # remove the arcs: pred(i) -> i and i -> succ(i)
    arcs_new.loc[pred_i,i]=0
    
    # insert i after j
    succ_j=successor(arcs_new,j)
    arcs_new.loc[j,succ_j]=0
    arcs_new.loc[j,i]=1
    arcs_new.loc[i,succ_j]=1
    
    numpy.fill_diagonal(arcs_new.values, 0) #set diagonal to zero
    return arcs_new.astype(int)

# <codecell>
def Initialize():
    global zips, allorders, orders, customers, customersplus, row_count,arcs, store_id, row_count,drive_time, dist_mtx, tg1removal
    if tg1removal==1: 
        orders=allorders[allorders.tg1==0].copy()
    else:
        orders=allorders.copy()
    orders=orders[(orders['Store_ID']==store_id) & (orders['Delivery Date']==date)] #select store and date
    orders['Customer_ID2']=1000+orders['Customer_ID'].astype('category').cat.codes #categorical coding of customer ID
    orders['Customer_ID2'].value_counts() #Frequency Table
    
    #merge with lat-long information from Sites.xlsx
    merged=orders.merge(zips[zips['Type']=='Store'], left_on='PA_Zip',right_on='Zip', how='left')
    #merged.Lat[merged.Lat.isnull()] #check
    merged=merged.rename(columns={"Lat":"PA_lat","Long":"PA_lon"})

    merged=merged.merge(zips[zips['Type']=='Customer'], left_on='RA',right_on='Zip', how='left')
    #merged.Lat[merged.Lat.isnull()] #check
    merged=merged.rename(columns={"Lat":"RA_lat","Long":"RA_lon",})    
    orders=merged
    
    customers=orders.Customer_ID2.values.tolist()
    customers=list(set(customers))
    row_count=len(customers)
    print 'Orders Count =',row_count
    customersplus=[store_id]+customers
        
    sites=orders[['Customer_ID2','RA_lat','RA_lon']]
    sites.drop_duplicates(inplace=True)
    sites['type']='C'
    store_info=[store_id,orders[orders.Store_ID==store_id].PA_lat.iloc[0],orders[orders.Store_ID==store_id].PA_lon.iloc[0],'S']
    sites.loc[-1]=store_info
    sites.columns=['site_id','lat','lon','type']
    sites.set_index('site_id', inplace=True)
    sites

    # Create Distance Matrix
    zeros_matrix = numpy.array([numpy.zeros(row_count+1)]*(row_count+1)).astype(int)
    dist_mtx=pandas.DataFrame(zeros_matrix,index=customersplus,columns=customersplus)
    customersplus[customersplus==1026]
    #print dist_mtx
    for i in customersplus:
        for j in customersplus:
            #print i
            #print j
            lat1=sites.loc[i,'lat']
            lon1=sites.loc[i,'lon']
            lat2=sites.loc[j,'lat']
            lon2=sites.loc[j,'lon']
            dist_mtx.loc[i,j]=dist(lat1,lon1,lat2,lon2)
        
    # Create Drive Time Matrix
    drive_time=dist_mtx.copy()
    #print drive_time
    for i in customersplus:
        for j in customersplus:
            d=dist_mtx.loc[i,j]
            #print 'i=',i,'j=',j, 'd=',d
            if d<50: 
                drive_time.loc[i,j]= d/36.0
            else: 
                drive_time.loc[i,j]= 50/36.0+(d-50)/70.0
    
    #<codecell>
    # create an initial solution
    zeros_matrix = numpy.array([numpy.zeros(row_count+1)]*(row_count+1)).astype(int)
    arcs=pandas.DataFrame(zeros_matrix,index=[store_id]+customers,columns=[store_id]+customers)
    #arcs.loc[1001,1002]
    arcs.iloc[0,1:]=1
    arcs.iloc[1:,0]=1


def neighborhood(arcs_mtx,size):
    global customers
    neighbors=list()
    ns=min(size,len(customers)-1) # corrected neighborhood size
    arcs_neighbor=arcs_mtx.copy() #initialize
    sample=random.sample(customers, 1+ns) #random select
    #print 'Candidate Customer:',sample[0]
    for j in sample[1:ns+1]:
        i=sample[0]
        
        relocated=relocate(arcs_neighbor,i,j).copy()
        #network(relocated)
        arcs_neighbor=relocated
        neighbors.append(arcs_neighbor)           
        
        print '.',

    return neighbors

# <codecell>

def SA():
    global temperature, beta, max_iterations, alpha, NSize, service_time, max_stops, wage,driving_cost_per_km,max_route_duration, customers, customersplus, best_arcs,best_cost
    global dist_mtx,drive_time,orders,row_count,sites,arcs, costs_achieved
    # Simulated Annealing
    print '\n- Begin SA'
    Initialize()
    #print 'row_count=', row_count
    originalNSize=NSize
    best_arcs=arcs.copy()
    best_cost=cost1(solution(best_arcs))
    current_arcs=best_arcs.copy()
    if len(best_arcs)<=2 : return best_arcs
    #network(best_arcs)
    costs_achieved=[best_cost]
    no_improvement=0
    empty_counter=0
        
    for iteration in range(max_iterations):
        no_improvement+=1
        if no_improvement>max_no_improvement: 
            print 'stopped due to no improvement'
            return best_arcs
        
        try:
            print 'Iteration',iteration, #' Best cost',best_cost,
            if no_improvement>2 and iteration<max_iterations/2: NSize=min(NSize+5,15)#expand search
            elif empty_counter/iteration>.3: NSize=min(NSize+5,15)
            else: NSize=originalNSize #narrow search
            
            l=neighborhood(current_arcs, NSize)   #create neighborhood
            if len(l)>0 and feasible(solution(l[0]))==1: l=l+neighborhood(l[0], NSize)
            
#            neig_size=len(l)
#            # neighborhood booster
#            if no_improvement>2: l=l+neighborhood(current_arcs, NSize)
                
            l_solutions=map(solution,l)
            l_feasible=map(feasible,l_solutions)
            
            #Delete infeasible neighbors 
            j=0
            for items in range(len(l)):
                if len(l)==0: break
                if l_feasible[j]==0: 
                    del l_solutions[j]
                    del l[j]
                else:
                    j+=1
            # IF NO FEASIBLE NEIGHBOR
            if len(l)==0: 
                print '(neighborhood empty)' 
                empty_counter+=1
                continue
            
            #network(current_arcs)
            rank=map(cost1,l_solutions)     #compute costs
            rank=pandas.DataFrame(rank,columns=['cost'])      #convert to DF
            rank.sort_values('cost', ascending=True, inplace=True)       #rank neighbors
            #print 'rank:', rank
            
            current_cost=rank.iloc[0]['cost']
            delta=best_cost-current_cost
            if delta>0:
                print 'Best Cost --->', current_cost
                best_arcs=l[rank.index[0]].copy()
                best_cost=current_cost
                current_arcs=l[rank.index[0]].copy()
                
                if delta>10**8: no_improvement=0 #reset counting of no improvements in NRoutes
                #network(best_arcs)
                #print '*'*50
            else:
                r=random.random()**beta 
                r=int(r*len(rank))    #select random integer between 0 and Nsize => randomly select a neighbor
                neighbor_cost=rank.iloc[r,0]
                
                delta=best_cost-neighbor_cost
                if delta>=0: 
                    current_arcs=l[rank.index[r]].copy()   #find the arcs matrix of the selected neighbor
                    current_cost=neighbor_cost
                elif random.random()<=exp(delta/temperature): 
                    current_arcs=l[rank.index[r]].copy()
                    current_cost=neighbor_cost
                    
            temperature=temperature*alpha
            #reporting
                        #print 'Current Cost:',current_cost
            costs_achieved.append(current_cost)
        except Exception as e: 
            print e
            print '#'
            pass
    #print '*'*50
    #l[rank.index[r]]
    #print best_arcs
    #print best_cost
    #plt.plot(costs_achieved)
    #plt.ylabel('cost')
    #plt.show()
    #network(best_arcs)
    return best_arcs
    

# <codecell>

# -*- coding: utf-8 -*-
"""
VRP algorithm, Part 2- Large neighborhood Search (LNS)
"""

# <codecell>

def TruckCounter(solution_mtx):
    # create a summary schedule of routes with start and end times
    solution_mtx=solution_mtx.sort_values(['route'])
    sss=solution_mtx.groupby(by='route')['route_start_time'].max()
    eee=solution_mtx.groupby(by='route')['route_end_time'].max()
    routes_schedule=pandas.concat([sss, eee], axis=1)
    routes_schedule=routes_schedule.sort_values(['route_start_time'])
    routes_schedule['truck']=0
    
    #assign routes to trucks
    for truck_count in range(1,len(routes_schedule)+1):
        end=0
        for r in routes_schedule.index: #find routes on the same truck
            candidates=routes_schedule[(routes_schedule['route_start_time']>end) & (routes_schedule['truck']==0)]
            if len(candidates)==0: break
            current_route=candidates.index[0]  
            routes_schedule.loc[current_route,'truck']=truck_count
            end=routes_schedule.loc[current_route,'route_end_time']
            #print 'current_route =',current_route

    return routes_schedule.truck.max()
# <codecell>

def relatedness(sol_df1, i,j):
    global wage, dist_mtx
    v=int(sol_df1.loc[i,'route']==sol_df1.loc[j,'route'])
    c=dist_mtx.loc[i,j]*driving_cost_per_km+drive_time.loc[i,j]*wage
    return 1/(c+v+.001)
    
def RemoveSet(base_list,removing_list):
    A=copy.copy(base_list)
    for i in removing_list:
        A.remove(i)
    return A

def SelectCustomers(solution_df, n):
    S1=[random.choice(customers)]  #initialize 
    for i in range(2,n+1):
        C=random.choice(S1)
        nonS1=RemoveSet(customers,S1) #remove S from customers
        relatedness_ranking=sorted(nonS1, key=functools.partial(relatedness, solution_df,C),reverse=True) #sort nonS1 customers based on relatedness to C
        rank_number=int(random.random()**beta*len(nonS1))  # Select a biased-random related node
        next_node=relatedness_ranking[rank_number]
        S1.append(next_node)
    return S1

def neighborhood2(arcs_mtx,S,size):
    #print S
    for retry in range(1):      # if no feasible neighborhood is found retry
        neighbors=list()
        nonS=RemoveSet(customers,S)
        ns=min(size,len(nonS)-1) # corrected neighborhood size
        arcs_neighbor=arcs_mtx.copy() #initialize
        
        for iter in range(ns):
            relocated=arcs_neighbor.copy()
            #print '... relocating',
            for fromm in S:
                to=random.choice(nonS)
                relocated=relocate(relocated,fromm,to).copy()
                
            sol1=solution(relocated)
        
            if feasible(sol1):
                #print '... relocating',fromm,'->',to,
                #print 'cost=', cost2(sol1)
                print '.',
                arcs_neighbor=relocated.copy()
                neighbors.append(arcs_neighbor)
                #network(relocated)
        if len(neighbors)>0: break
    #print 'Neighborhood Size=',len(neighbors)
    return neighbors


# <codecell>

def LNS(arcs_mtx):
    global costs_achievedII,max_no_improvement
    best_arcsII=arcs_mtx.copy()
    best_costII=cost2_from_arcs(best_arcsII)
    costs_achievedII=[best_costII]
    if len(best_arcsII)<=2 : return best_arcsII #if there is only 1 order there is no need to solve
    print '\n- Begin LNS','\nstarting_cost=',best_costII
    #network(best_arcsII)
    
    no_improvement=0 # counter of number of iterations with less than 1% imrovement
    
    for iter1 in range(maxiter1):
        for N in range(1,p+1):
            for iter2 in range(maxiter2):
                if no_improvement>max_no_improvement: 
                    print 'stopped due to no improvement!'
                    return best_arcsII # stop if there is less than 1% imrovement improvement after 10 iterations
            
                print 'IterI=','N=',N,iter1,'IterII=',iter2,
                S=SelectCustomers(solution(best_arcsII), N)
                sizeLNS=min(4,len(customers)*len(S))
                
                NR=neighborhood2(arcs_mtx,S,sizeLNS)
                
                if len(NR)==0: 
                    no_improvement+=1
                    continue
                NR_ranking=sorted(NR, key=cost2_from_arcs) #sort by cost 
                NR_best_cost=cost2_from_arcs(NR_ranking[0]) # select the best
                delta=best_costII-NR_best_cost
                
                if delta>0:
                    print 'Best Cost Improved:',best_costII,
                    best_arcsII=NR_ranking[0]
                    best_costII=NR_best_cost
                    print '->',best_costII
                    #network(best_arcsII)
                    costs_achievedII.append(best_costII)
                
                
                percent_improv=delta/(best_costII%(10**6)+0.000001)
                #print 'delta=',delta,'denom=',(best_costII%(10**6)+0.000001),percent_improv,'% improved'
                if percent_improv<0.01 : no_improvement+=1 #count if there is less than 1% imrovement 
                else: no_improvement=0
                print 'no_improvement=',no_improvement
    print 'best_cost=',best_costII
    print costs_achievedII
    #plt.plot(costs_achievedII)
    #plt.ylabel('cost')
    #plt.show()
    return best_arcsII

# <codecell>
##############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################
# <codecell>
# Break the original large file by store.
.
file_address='DID-Reduced Size Super Small.xlsx'
allorders=pandas.read_excel(file_address)
#stores=allorders.Store_ID.drop_duplicates().tolist()
#stores=sorted(stores)
stores=[2, 3, 4, 6, 7, 8, 9, 11, 13, 203, 212, 460, 462, 463, 464]


for store_id in stores:
    subset=allorders[allorders.Store_ID==store_id].copy()
    subset['Delivery Date']=subset['Delivery Date'].astype(str)
    name= 'DID_'+str(store_id)+'.csv'
    subset.to_csv(name)

# <codecell>
# create summary dataframes to store results

for store_id in stores:
    filename='DID_'+str(store_id)
    allorders=pandas.read_csv(filename+'.csv')
    orders=allorders.copy()
    
    #only select dates in which there are tg1 orders
    alldates=allorders['Delivery Date'].drop_duplicates()
    selected_dates=allorders.groupby(by=['Delivery Date'])['tg1'].max()
    selected_dates=selected_dates[selected_dates==1].index
    allorders=allorders[allorders['Delivery Date'].isin(selected_dates)]
    
    #create teh summary data frame
    iterables = [selected_dates,[0,1]]
    ind=pandas.MultiIndex.from_product(iterables, names=['date', 'tg1removal']) # create a multilevel index
    
    summary=pandas.DataFrame(index=ind,columns=['store_id','raw_cost','actual_cost','N','status','Nroutes','NTrucks','CPO', 'ts'])
    summary['status']=' '
    summary['store_id']=store_id
    summary.sort_index()
    summary.to_csv(filename+'_summary.csv')
    summary.head()


# <codecell>
# Main Loop

zips=pandas.read_csv('Sites.csv')
stores=[2, 3, 4, 6, 7, 8, 9, 11, 13, 203, 212, 460, 462, 463, 464]


for store_id in stores:
    Save=0
    progress=0
    CPO=-1    
    #load output file
    filename='DID_'+str(store_id)
    summary=pandas.read_csv(filename+'_summary.csv', index_col=[0,1])
    summary.head()
    selected_dates=summary.index.get_level_values('date').unique().tolist() #get unique dates
    total_iterations=len(selected_dates)*2
    allorders=pandas.read_csv(filename+'.csv')
    orders=allorders.copy()

    for date in selected_dates:
        NSize=5
        progress+=2
        percent_progress=round(100*progress/(total_iterations+0.0),2)
        
        #skip if there is no tg1
        try:
            y=allorders[(allorders['Store_ID']==store_id) & (allorders['Delivery Date']==date)]
            #print 'y=',y
            if len(y)==0:
                node1=(date,0)
                node2=(date,1)
                summary.loc[node1,'N']= summary.loc[node2,'N']=0
                summary.loc[node1,'status']=summary.loc[node2,'status']='skipped_no orders'
                print 'skipped_no orders on', date, 'for store',store_id
                continue
            if y['tg1'].max()==0:
                node1=(date,0)
                node2=(date,1)
                summary.loc[node1,'status']=summary.loc[node2,'status']='skipped_no tg1 orders'
                print 'skipped_no tg1 orders on', date
                continue
                           
            
            for tg1removal in [0,1]:
                print '*'*100,'\n','*'*5,'store_id=',store_id,'date=',date, '(', len(selected_dates), 'days',percent_progress,'%)','\n','*'*100
                node=(date,tg1removal)
                #summary.loc[node,'status']
                Save=1 #Save indicator
                if summary.loc[node,'status']!=' ': 
                    Save=0
                    continue #skip if it is already done
                ts1=time.time() #time stamp
                b=LNS(SA())
                ts2=time.time() #time stamp
          
                ss=solution(b)
                c=cost2(ss)
                actual_cost=c%10**6
                Nroutes=ss.route.max()
                NTrucks=TruckCounter(ss)
                NOrders=len(b)-1
                CPO=actual_cost/NOrders
                
                summary.loc[node,'raw_cost']=c
                summary.loc[node,'actual_cost']=actual_cost
                summary.loc[node,'N']=NOrders
                summary.loc[node,'status']='Complete'
                summary.loc[node,'Nroutes']=Nroutes
                summary.loc[node,'NTrucks']=NTrucks
                summary.loc[node,'CPO']=CPO
                summary.loc[node,'ts']=ts2-ts1
                network(b)
                print '\nStore:',store_id,'\tdate:',date, 'tg1removal:',tg1removal, 'N:',NOrders,'Nroutes:',Nroutes,'NTrucks:',NTrucks,'actual_cost:',actual_cost,'CPO:',CPO
                #print summary[node:]
                
        except Exception as e:
            #raise ValueError('A very specific bad thing happened')
            node1=(date,0)
            node2=(date,1)
            print 'Error on', date, e
            summary.loc[node1,'status']=summary.loc[node2,'status']=e
            
            pass
        
        
        if Save==1: #Save the file every time a new date is computed
            try: 
                summary.to_csv(filename+'_summary.csv')
            except: 
                summary.to_csv('C:\\Users\\sgolara.ASURITE.000\\Desktop\\VRP\\'+filename+'_summary.csv')
                print 'Dropbox locked. saved on desktop.'
            print '(saved)','(',progress,'/',total_iterations,')'
        
        
    print summary.status.head()
    summary.to_csv(filename+'_summary.csv')
    summary.to_csv('backup'+filename+'_summary.csv')

