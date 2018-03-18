#load required libraries
from igraph import *
from gurobipy import *
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import random
from sklearn.metrics import roc_curve, auc

###################################################################################################################
#default arguments
args = {}
FOLDS = 5
KO_FILE = './reimand.txt'
PPI_FILE = './yeast_ppi.txt'
PDI_FILE = './yeast_pdi.txt'
KPI_FILE = './yeast_kpi.txt'
KEGG_FILE = './yeast_kegg.txt'
Q = 0.001
LFC = 2
ROUNDS = 10
VARIANT = 'ASP'

#defining empty model
y = {}
r = {}
c = {}
x = {}
f = {}
a = {}
o = {}
m = Model()
####################################################################################################################

#helper function to read the yeast interactions
def read_net(filepathPPI, filepathPDI, filepathKPI, filepathKEGG):
    net = None
    pdi = None
    kpi = None
    kegg = None
    if filepathPPI is not None:
        net = pd.read_table(filepath_or_buffer=filepathPPI,delimiter='\t',header=None,skiprows=1)
        net.iloc[:,0] = net.iloc[:,0].astype(str)
        net.iloc[:,1] = net.iloc[:,1].astype(str)
        
    if filepathPDI is not None:
        pdi = pd.read_table(filepath_or_buffer=filepathPDI,delimiter='\t',header=None,skiprows=1)
        if filepathPDI == './S_cerevisiae-pdi2.net':
            pdi.iloc[:,0] = pdi.iloc[:,0].astype(int)
            pdi.iloc[:,0] = pdi.iloc[:,0].astype(str)
            pdi.iloc[:,1] = pdi.iloc[:,1].astype(int)
            pdi.iloc[:,1] = pdi.iloc[:,1].astype(str)
    if filepathKPI is not None:
        kpi = pd.read_table(filepath_or_buffer=filepathKPI,delimiter='\t',header=None,skiprows=1)
        if filepathPDI == './S_cerevisiae-kpi2.net':
            kpi.iloc[:,0] = kpi.iloc[:,0].astype(int)
            kpi.iloc[:,0] = kpi.iloc[:,0].astype(str)
            kpi.iloc[:,1] = kpi.iloc[:,1].astype(int)
            kpi.iloc[:,1] = kpi.iloc[:,1].astype(str)
    if filepathKEGG is not None:
        kegg = pd.read_table(filepathKEGG,delimiter='\t',header=None,skiprows=1)
    
    
    if net is not None:
        total = net
        total = pd.concat([total,pdi])
        total = pd.concat([total,kpi])
        total = pd.concat([total,kegg])
    elif pdi is not None:
        total = pdi
        total = pd.concat([total,net])
        total = pd.concat([total,kpi])
        total = pd.concat([total,kegg])
    elif kpi is not None:
        total = kpi
        total = pd.concat([total,net])
        total = pd.concat([total,pdi])
        total = pd.concat([total,kegg])
    elif kegg is not None:
        total = kegg
        total = pd.concat([total,net])
        total = pd.concat([total,pdi])
        total = pd.concat([total,kpi])
        
    vset = map(str,pd.unique(total.iloc[:,0:2].values.ravel()))
    G = Graph()
    G.add_vertices(len(vset))
    G.vs['name'] = vset
    Edge_sign = {}
    Edge_orient = {}
    Edge_type = {}
    Edge_so = {}
    if net is not None:
        for i in tqdm(range(net.shape[0])):
            u = None
            v = None
            try:
                u = G.vs.find(net.iloc[i,0]).index 
                v = G.vs.find(net.iloc[i,1]).index
            except IndexError:
                next
            if G.vs[u]['name'] != G.vs[v]['name']:
                ss = net.iloc[i,2] if str(net.iloc[i,2]) != 'nan' else None
                if ss is not None:
                    Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = (G.vs[u]['name'],G.vs[v]['name'])
                else:
                    Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = ss
                Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'ppi'

    if pdi is not None:
        for i in tqdm(range(pdi.shape[0])):
            u = None
            v = None
            try:
                u = G.vs.find(pdi.iloc[i,0]).index 
                v = G.vs.find(pdi.iloc[i,1]).index
            except IndexError:
                next
            ss = pdi.iloc[i,2] if str(pdi.iloc[i,2]) != 'nan' else None
            if G.vs[u]['name'] != G.vs[v]['name']:
                if frozenset([G.vs[u]['name'],G.vs[v]['name']]) in Edge_type.keys():
                    if Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] == (G.vs[v]['name'],G.vs[u]['name']):
                        Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    
                    if Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] != ss:
                        Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                        Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    
                    if Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] != 'pdi':
                        Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                else:
                    Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = ss
                    if ss is not None:
                        Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = G.vs[u]['name']
                    else:
                        Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = (G.vs[u]['name'],G.vs[v]['name'])
                    Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'pdi'
    
    if kpi is not None:
        for i in tqdm(range(kpi.shape[0])):
            u = None
            v = None
            try:
                u = G.vs.find(kpi.iloc[i,0]).index 
                v = G.vs.find(kpi.iloc[i,1]).index
            except IndexError:
                next
            ss = kpi.iloc[i,2] if str(kpi.iloc[i,2]) != 'nan' else None
            if G.vs[u]['name'] != G.vs[v]['name']:
                if frozenset([G.vs[u]['name'],G.vs[v]['name']]) in Edge_type.keys():
                    if Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] == (G.vs[v]['name'],G.vs[u]['name']):
                        Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                  
                    if Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] != ss:
                        Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                        Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    
                    Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'kpi'
                else:
                    Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = ss
                    if ss is not None:
                        Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = (G.vs[u]['name'],G.vs[v]['name'])
                    else:
                        Edge_orient[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'kpi'

    if kegg is not None:
        for i in tqdm(range(kegg.shape[0])):
            u = None
            v = None
            try:
                u = G.vs.find(kegg.iloc[i,0]).index
                v = G.vs.find(kegg.iloc[i,1]).index
            except IndexError:
                next
            ss = None 
            if str(kegg.iloc[i,3]) != 'indirect effect': 
                if str(kegg.iloc[i,2]) == 'activation':
                    ss = 0
                elif str(kegg.iloc[i,2]) == 'inhibition':
                    ss = 1
            else:
                next
            
            if G.vs[u]['name'] != G.vs[v]['name'] and str(kegg.iloc[i,3]) != 'indirect effect':
                if frozenset([G.vs[u]['name'],G.vs[v]['name']]) in Edge_type.keys():
                    if ss is not None:
                        Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = ss
                        Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                        Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'kegg'
                else:
                    Edge_sign[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = ss
                    Edge_so[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = None
                    Edge_type[frozenset([G.vs[u]['name'],G.vs[v]['name']])] = 'kegg'
    
    G.add_edges(map(tuple,Edge_type.keys()))
    G.es['name'] = Edge_type.keys()
    G.es['sign'] = Edge_sign.values()
    G.es['orient'] = [None]*len(Edge_sign.values())
    G.es['type'] = Edge_type.values()
    G.es['parent'] = Edge_so.values()
                
    return(G)



###################################################################################################################
#Compute shortest path subgraphs Gs using BFS
def sp_graph(G,s,Ts):
    ds = dict([(v,G.vcount()) for v in G.vs['name']])
    dTs = {}
    for v in G.bfsiter(vid = G.vs.find(s).index,advanced=True):
        ds[v[0]['name']] = v[1]
    for t in Ts:
        dt = dict([(v,G.vcount()) for v in G.vs['name']])
        for v in G.bfsiter(vid = G.vs.find(t).index,advanced=True):
            dt[v[0]['name']] = v[1]
        dTs[t] = dt
    Gs = Graph()
    Gs.vs['name'] = []
    Gs.es['name'] = []
    for e in G.es:
        u = list(e['name'])[0]
        v = list(e['name'])[1]
        sat = [(ds[u] + 1 + dTs[t][v]) == ds[t] or (ds[v] + 1 + dTs[t][u]) == ds[t] for t in Ts]
        if any(sat):
            if u not in Gs.vs['name']:
                Gs.add_vertex(u, layer = ds[u])
            if v not in Gs.vs['name']:
                Gs.add_vertex(v, layer = ds[v])
            Gs.add_edge(u,v,name = frozenset([u,v]),orient = e['orient'],sign = e['sign'])
    return Gs,ds,dTs

###################################################################################################################
#Helper function that calls the shortest path subgraph function to generate shortest path subgraphs per source 
#given a graph, set of sources and associated targets as input
def preprocess(G,S,P):
    H = {}
    R = {}
    TNP = 0
    for s in tqdm(np.unique(S)):
        res = sp_graph(G,s,P[s])
        if len(res[0].vs['name']) > 1:
            H[s] = res
            R[s] = np.intersect1d(P[s],H[s][0].vs['name'])
            TNP += len(R[s])
    print '%d pairs in play'%(TNP)
    return G,H,R

###################################################################################################################
#helper to define ILP variables in Gurobi 7.0.2 (Free Academic Lisence)
def setVars(H,R, variant):
    if variant == 'AP' or variant == 'ASP' or variant == 'AllSP':
        SS = H.keys()
        for s in tqdm(SS):
            GG = H[s][0]
            for e in GG.es:
                u = list(e['name'])[0]
                v = list(e['name'])[1]
                assert np.abs(H[s][1][u] - H[s][1][v]) == 1
                if e['name'] not in x.keys():
                    x[frozenset([u,v])] = m.addVar(vtype = GRB.BINARY, name = 'x_%s%s'%(u,v))
            for v in GG.vs['name']:
                c[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'c%s_%s'%(s,v))
                r[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'r%s_%s'%(s,v))
                if v in R[s]:
                    y[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'y%s_%s'%(s,v))
    
    elif variant == 'AdirSP':
        SS = H.keys()
        for s in tqdm(SS):
            GG = H[s][0]
            for e in GG.es:
                u = list(e['name'])[0]
                v = list(e['name'])[1]
                assert np.abs(H[s][1][u] - H[s][1][v]) == 1
                f[(s,v,u)] = m.addVar(vtype = GRB.BINARY, name = 'f%s_%s_%s'%(s,v,u))
                f[(s,u,v)] = m.addVar(vtype = GRB.BINARY, name = 'f%s_%s_%s'%(s,u,v))
                a[(s,u,v)] = m.addVar(vtype = GRB.BINARY, name = 'a%s_%s_%s'%(s,u,v))
                a[(s,v,u)] = m.addVar(vtype = GRB.BINARY, name = 'a%s_%s_%s'%(s,v,u))
                if e['name'] not in x.keys():
                    o[(u,v)] = m.addVar(vtype = GRB.BINARY, name = 'o%s_%s'%(u,v))
                    o[(v,u)] = m.addVar(vtype = GRB.BINARY, name = 'o%s_%s'%(v,u))
                    x[frozenset([u,v])] = m.addVar(vtype = GRB.BINARY, name = 'x_%s%s'%(u,v))
            for v in GG.vs['name']:
                c[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'c%s_%s'%(s,v))
                r[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'r%s_%s'%(s,v))
                if v in R[s]:
                    y[(s,v)] = m.addVar(vtype = GRB.BINARY, name = 'y%s_%s'%(s,v))

###################################################################################################################
#helper to define ILP constraints in Gurobi 7.0.2 (Free Academic Lisence)        
def setConstr(H,R, variant, sign):
    SS = H.keys()
    if variant == 'AP' or variant == 'ASP':
        for s in tqdm(SS):
            GG = H[s][0]
            ds = H[s][1]
            m.addConstr(r[(s,s)] == 0, name = 'r%s_%s = 0'%(s,s))
            m.addConstr(c[(s,s)] == 0, name = 'c%s_%s = 0'%(s,s))
            for e in GG.es:
                l = list(e['name'])
                u = l[0]
                v = l[1]

                if ds[u] == ds[v] + 1:    
                    m.addConstr(r[(s,u)] - c[(s,u)] <= 2 - r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= 2 - r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,v))
                    m.addConstr(r[(s,u)] - c[(s,u)]  <= r[(s,v)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= x_%s%s + r_%s%s'%(s,u,s,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] >= r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s >= x_%s%s - r_%s%s'%(s,u,s,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] >= x[frozenset([u,v])] - r[(s,v)],
                                name = 'r_%s%s + c_%s%s >= r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,v))

                else:   
                    m.addConstr(r[(s,v)] - c[(s,v)] <= 2 - r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= 2 - r_%s%s - x_%s%s'%(s,v,s,v,s,u,u,v))
                    m.addConstr(r[(s,v)] - c[(s,v)]  <= r[(s,u)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= x_%s%s + r_%s%s'%(s,v,s,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] >= r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s >= x_%s%s - r_%s%s'%(s,v,s,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] >= x[frozenset([u,v])] - r[(s,u)],
                                name = 'r_%s%s + c_%s%s >= r_%s%s - x_%s%s'%(s,v,s,v,s,u,u,v))

            for t in R[s]:
                m.addConstr(c[(s,t)] + y[(s,t)] <= 1, name = 'c_%s%s + y_%s%s <= 1'%(s,t,s,t))
                m.addConstr(r[(s,t)] == sign[(s,t)], name = 'r%s_%s = %d'%(s,t,sign[(s,t)]))

            for v in GG.vs['name']:
                if ds[v] >= 1:
                    Np = [u for u in GG.vs[GG.neighbors(GG.vs.find(v).index)]['name'] if ds[u]+1==ds[v]]
                    m.addConstr(quicksum(c[(s,u)] - 1 for u in Np)+1 <= c[(s,v)],
                                name = 'c_%s%s atleast1'%(s,v))

    elif variant == 'AllSP':
        for s in tqdm(SS):
            GG = H[s][0]
            ds = H[s][1]
            m.addConstr(r[(s,s)] == 0, name = 'r%s_%s = 0'%(s,s))
            m.addConstr(c[(s,s)] == 0, name = 'c%s_%s = 0'%(s,s))
            for e in GG.es:
                l = list(e['name'])
                u = l[0]
                v = l[1]
                if ds[u] == ds[v] + 1:
                    m.addConstr(c[(s,v)] <= c[(s,u)],
                                name = 'c_%s%s <= c_%s%s'%(s,v,s,u))
                    m.addConstr(r[(s,u)] - c[(s,u)] <= 2 - r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= 2 - r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,v))
                    m.addConstr(r[(s,u)] - c[(s,u)] <= r[(s,v)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= x_%s%s + r_%s%s'%(s,u,s,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] >= r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s >= x_%s%s - r_%s%s'%(s,u,s,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] >= x[frozenset([u,v])] - r[(s,v)],
                                name = 'r_%s%s + c_%s%s >= r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,v))


                else:
                    m.addConstr(c[(s,u)] <= c[(s,v)],
                               name = 'c_%s%s <= c_%s%s'%(s,u,s,v))
                    m.addConstr(r[(s,v)] - c[(s,v)] <= 2 - r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= 2 - r_%s%s - x_%s%s'%(s,v,s,v,s,u,u,v))
                    m.addConstr(r[(s,v)] - c[(s,v)] <= r[(s,u)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s <= x_%s%s + r_%s%s'%(s,v,s,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] >= r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s >= x_%s%s - r_%s%s'%(s,v,s,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] >= x[frozenset([u,v])] - r[(s,u)],
                                name = 'r_%s%s + c_%s%s >= r_%s%s - x_%s%s'%(s,v,s,v,s,u,u,v))


            for t in R[s]:
                m.addConstr(c[(s,t)] + y[(s,t)] <= 1, name = 'c_%s%s + y_%s%s <= 1'%(s,t,s,t))
                m.addConstr(r[(s,t)] == sign[(s,t)], name = 'r%s_%s = %d'%(s,t,sign[(s,t)]))

    elif variant == 'AdirSP':
        for e in tqdm(G.es):
            l = list(e['name'])
            u = l[0]
            v = l[1]
            try:
                #constr_o[e['name']] = m.addConstr(o[e['orient']] == 1, name = 'orient %s%s'%(tuple(e['name'])[0],tuple(e['name'])[1]))
                m.addConstr(o[(u,v)]+o[(v,u)] == 1,name='o_%s%s + o_%s%s = 1'%(u,v,v,u))
            except KeyError:
                next
        for s in tqdm(SS):
            GG = H[s][0]
            ds = H[s][1]
            #dTs = H[s][2]
            m.addConstr(r[(s,s)] == 0, name = 'r%s_%s = 0'%(s,s))
            m.addConstr(c[(s,s)] == 0, name = 'c%s_%s = 0'%(s,s))
            for e in GG.es:
                l = list(e['name'])
                u = l[0]
                v = l[1]
                m.addConstr(f[(s,v,u)] <= o[(v,u)], name = 'f_%s%s%s <= o_%s%s'%(s,u,v,u,v))
                m.addConstr(f[(s,u,v)] <= o[(u,v)], name = 'f_%s%s%s <= o_%s%s'%(s,u,v,u,v))
                if ds[u] == ds[v] + 1:
                    if ds[u] >= 2:
                        Np = [w for w in GG.vs[GG.neighbors(GG.vs.find(v).index)]['name'] if ds[w]+1==ds[v]]
                        m.addConstr(f[(s,v,u)] <= quicksum(f[(s,w,v)] for w in Np), name = 'flow cons at %s %s %s'%(s,v,u))

                    m.addConstr(a[(s,v,u)] >= c[(s,v)],
                                name = 'a_%s%s%s <= c_%s%s'%(s,v,u,s,v))
                    m.addConstr(a[(s,v,u)] >=1-f[(s,v,u)],
                                name = 'a_%s%s%s <= 1-f_%s%s%s'%(s,u,v,s,u,v))
                    m.addConstr(a[(s,v,u)] <= 1- f[(s,v,u)]+c[(s,v)],
                                name = 'a_%s%s%s>= 1-f_%s%s%s + c_%s%s'%(s,v,u,s,v,u,s,v))    
                    m.addConstr(r[(s,u)] - c[(s,u)] - 1 + f[(s,v,u)] <= 2 - r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s - 1 + f_%s%s%s <= 2 - r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,s,v,u,v))
                    m.addConstr(r[(s,u)] - c[(s,u)] - 1 + f[(s,v,u)]  <= r[(s,v)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s - 1 + f_%s%s%s<= x_%s%s + r_%s%s'%(s,u,s,u,s,v,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] + 1 - f[(s,v,u)] >= r[(s,v)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s + 1 - f_%s%s%s>= x_%s%s - r_%s%s'%(s,u,s,u,s,v,u,u,v,s,v))
                    m.addConstr(r[(s,u)] + c[(s,u)] + 1 - f[(s,v,u)] >= x[frozenset([u,v])] - r[(s,v)],
                                name = 'r_%s%s + c_%s%s + 1 - f_%s%s%s>= r_%s%s - x_%s%s'%(s,u,s,u,s,v,u,s,v,u,v))

                else:
                    if ds[v] >= 2:
                        Np = [w for w in GG.vs[GG.neighbors(GG.vs.find(u).index)]['name'] if ds[w]+1==ds[u]]
                        m.addConstr(f[(s,u,v)] <= quicksum(f[(s,w,u)] for w in Np), name = 'flow cons at %s %s %s'%(s,u,v))
                    m.addConstr(a[(s,u,v)] >= c[(s,u)],
                                name = 'a_%s%s%s <= c_%s%s'%(s,u,v,s,u))
                    m.addConstr(a[(s,u,v)] >=1-f[(s,u,v)],
                                name = 'a_%s%s%s <= 1-f_%s%s%s'%(s,u,v,s,u,v))
                    m.addConstr(a[(s,u,v)] <= 1- f[(s,u,v)]+c[(s,u)],
                                name = 'a_%s%s%s>= 1-f_%s%s%s + c_%s%s'%(s,u,v,s,u,v,s,u))    
                    m.addConstr(r[(s,v)] - c[(s,v)] - 1 + f[(s,u,v)] <= 2 - r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s - 1 + f_%s%s%s <= 2 - r_%s%s - x_%s%s'%(s,v,s,v,s,u,v,s,u,u,v))
                    m.addConstr(r[(s,v)] - c[(s,v)] - 1 + f[(s,u,v)]  <= r[(s,u)] + x[frozenset([u,v])],
                                name = 'r_%s%s - c_%s%s - 1 + f_%s%s%s<= x_%s%s + r_%s%s'%(s,v,s,v,s,u,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] + 1 - f[(s,u,v)] >= r[(s,u)] - x[frozenset([u,v])],
                                name = 'r_%s%s + c_%s%s + 1 - f_%s%s%s>= x_%s%s - r_%s%s'%(s,v,s,v,s,u,v,u,v,s,u))
                    m.addConstr(r[(s,v)] + c[(s,v)] + 1 - f[(s,u,v)] >= x[frozenset([u,v])] - r[(s,u)],
                                name = 'r_%s%s + c_%s%s + 1 - f_%s%s%s>= r_%s%s - x_%s%s'%(s,v,s,v,s,u,v,s,u,u,v))

            for t in R[s]:
                m.addConstr(c[(s,t)] + y[(s,t)] <= 1, name = 'c_%s%s + y_%s%s <= 1'%(s,t,s,t))
                m.addConstr(r[(s,t)] == sign[(s,t)], name = 'r%s_%s = %d'%(s,t,sign[(s,t)]))

            for v in GG.vs['name']:
                if ds[v] >= 1:
                    Np = [u for u in GG.vs[GG.neighbors(GG.vs.find(v).index)]['name'] if ds[u]+1==ds[v]]
                    m.addConstr(quicksum(a[(s,u,v)] - 1 for u in Np)+1 <= c[(s,v)],
                                name = 'c_%s%s sat'%(s,v))
        
    
    m.update()
 
###################################################################################################################   
#helper function to contract an input graph into an acyclic graph (A-path, Houri et. al. 2012)
def contract(G):
    res = G.biconnected_components()
    aa = G.cut_vertices()
    aa = G.vs[aa]['name']
    F1 = []
    
    for gg in res.subgraphs():
        if gg.ecount() >= 3:
            assert all([d >= 2 for d in gg.degree()])
            if all([e['sign'] == 0 for e in gg.es]):
                print gg
            for e in gg.es:
                F1.append(e['name'])
    
    
    F = [e for e in G.es if e['name'] not in F1]
    
    G1 = Graph()
    G1.vs['name'] = []
    block = 0
    for gg in res.subgraphs():
        if gg.ecount() >= 3:
            G1.add_vertex('block%d'%(block))
            for u in gg.vs['name']:
                if u in aa:
                    if u not in G1.vs['name']:
                        G1.add_vertex(u)
                        G1.add_edge(u,'block%d'%(block),name = frozenset(['block%d'%(block),u]),sign = None,orient=None,type='block')
                    
             
            block += 1
            
    for e in tqdm(F):
        l = list(e['name'])
        u = l[0]
        v = l[1]
        if u not in G1.vs['name'] and v not in G1.vs['name']:
            G1.add_vertex(u)
            G1.add_vertex(v)
            G1.add_edge(u,v,sign = e['sign'],orient = e['orient'], name = e['name'],type = e['type'])
        elif u in G1.vs['name'] and v not in G1.vs['name']:
            G1.add_vertex(v)
            G1.add_edge(u,v,sign = e['sign'],orient = e['orient'], name = e['name'],type = e['type'])
        elif u not in G1.vs['name'] and v in G1.vs['name']:
            G1.add_vertex(u)
            G1.add_edge(u,v,sign = e['sign'],orient = e['orient'], name = e['name'],type = e['type'])
        else:
            G1.add_edge(u,v,sign = e['sign'],orient = e['orient'], name = e['name'],type = e['type'])
    
    to_del = [v['name'] for v in G1.vs if len(G1.neighbors(v.index)) == 0]
    G1.delete_vertices(to_del)
    return F,G1

###################################################################################################################
def main():
    G = read_net(PPI_FILE,PDI_FILE,KPI_FILE,KEGG_FILE)
    #Uncomment for A-path
    #F,G = contract(G)
    pairs = pd.read_table(KO_FILE)
    pairs = pairs[pairs.iloc[:,3]<Q]
    pairs = pairs[np.array(map(np.abs,pairs.iloc[:,2]))>LFC]
    S = map(np.str,pairs.iloc[:,0])
    T = map(np.str,pairs.iloc[:,1])
    filter = [idx for idx in tqdm(range(len(T))) if S[idx] in G.vs['name'] and T[idx] in G.vs['name'] and S[idx] != T[idx]]
    S = [S[idx] for idx in filter]
    T = [T[idx] for idx in filter]
    M = [pairs.iloc[idx,2] for idx in filter]
    P = {}
    for s in np.unique(S):
        P[s] = [T[idx] for idx in np.where(np.array(S) == s)[0]]
    sign = {}
    weight = {}
    for i in tqdm(range(len(T))):
        sign[(S[i],T[i])] = 0 if M[i] < 0 else 1
        
    G, H, R = preprocess(G,S,P)
    setVars(H,R,VARIANT)
    setConstr(H,R,VARIANT,sign)


    #prediction in F folds
    EE = [e for e in G.es if e['name'] in x.keys() and e['sign'] is not None]
    fx = {}
    for e in EE:
        fx[e['name']] = 0
    accuracy = []
    random.shuffle(EE)
    N = len(EE)
    step = (N/FOLDS)
    N = step*FOLDS
    test = []
    for i in range(step,(N+step),step):
        C = {}
        test = range((i-step),i)
        print test
        train = np.setdiff1d(range(N),test)
        m.reset()
        for j in train:
            C[EE[j]['name']] = m.addConstr(x[EE[j]['name']] == EE[j]['sign'])
        m.update()
        for round in range(ROUNDS):
            m.reset()
            m.setObjective(quicksum((1 + 0.1*np.random.randn())*y[(s,t)] for s,t in y.keys()),GRB.MAXIMIZE)
            m.setParam(GRB.param.TimeLimit, 14400)
            m.setParam(GRB.param.MIPGap, 0.1)
            m.optimize()
            for j in test:
                fx[EE[j]['name']] += x[EE[j]['name']].X
        m.reset()
        for j in train:
            m.remove(C[EE[j]['name']])
        m.update()


    print 'Validating predictions...'
    predd = []
    label = []
    pi_type = []
    for j in tqdm(range(len(EE))):
        e = EE[j]
        if e['name'] in fx.keys() and e['sign'] is not None:
            predd.append(fx[e['name']])
            label.append(e['sign'])
            if e['type'] == 'pdi':
                pi_type.append('pdi')
            elif e['type'] == 'kpi':
                pi_type.append('kpi')
            elif e['type'] == 'kegg':
                pi_type.append('kegg')

    label = np.array(label)
    predd = np.array(predd)/ROUNDS
    pi_type = np.array(pi_type)

    fpr, tpr, thresholds = roc_curve(label, predd, pos_label=1)
    print 'Performance (AUC) All (%d +, %d -): %.2f'%(len(label)-sum(label), sum(label),  auc(fpr, tpr))
    fpr, tpr, thresholds = roc_curve(label[pi_type == 'pdi'], predd[pi_type == 'pdi'], pos_label=1)
    print 'Performance (AUC) PDI (%d +, %d -): %.2f'%(len(label[pi_type == 'pdi'])-sum(label[pi_type == 'pdi']), sum(label[pi_type == 'pdi']),  auc(fpr, tpr))
    fpr, tpr, thresholds = roc_curve(label[pi_type == 'kpi'], predd[pi_type == 'kpi'], pos_label=1)
    print 'Performance (AUC) KPI (%d +, %d -): %.2f'%(len(label[pi_type == 'kpi'])-sum(label[pi_type == 'kpi']), sum(label[pi_type == 'kpi']),  auc(fpr, tpr))
    fpr, tpr, thresholds = roc_curve(label[pi_type == 'kegg'], predd[pi_type == 'kegg'], pos_label=1)
    print 'Performance (AUC) KEGG (%d +, %d -): %.2f'%(len(label[pi_type == 'kegg'])-sum(label[pi_type == 'kegg']), sum(label[pi_type == 'kegg']),  auc(fpr, tpr))
    dat = pd.DataFrame({'pred':predd,'label':label,'type':pi_type})
    print dat.iloc[:10,]
    dat.to_csv('./output.txt',sep = '\t', index = False)


def getopts(argv):
    opts = {}  
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
        argv = argv[1:]
    return opts


if __name__ == '__main__':
    from sys import argv
    args = getopts(argv)
    if '-i' in args:
        KO_FILE = args['-i']
    if '-v' in args:
        VARIANT = args['-v']
    if '-f' in args:
        FOLDS = int(args['-f'])
    if '-r' in args:
        ROUNDS = int(args['-r'])
    main()
