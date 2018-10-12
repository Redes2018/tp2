import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import os
#-----------------------------------------------------
#               FUNCIONES PARA EL TP2:
#-----------------------------------------------------
def overlap(G,F):

#Le pasamos dos grafos G y F y nos devuelve
#cantidad de enlaces en comun nenlaces_comunes tipo int

    nodos_G=set(G.nodes)
    nodos_F=set(F.nodes)

    nodos_FG=nodos_G.intersection(nodos_F)
    nodos_FG=list(nodos_FG)

    G=np.array(nx.to_numpy_matrix(G,nodelist=nodos_FG))
    F=np.array(nx.to_numpy_matrix(F,nodelist=nodos_FG))

    suma=G+F

    #Buscamos los lugares donde suma sea 2
    nenlaces_comunes=len(np.where(suma==2)[0])
    
    return nenlaces_comunes


#-------------------------------------------------------
def esenciales(G,ess):
    
# A partir de un grafo G y una lista de proteinas esenciales ess
# devuelve el grafo con el atributo essential (True o False) y dos listas:
# de las proteinas esenciales y las no esenciales de G
    nodos_G = set(G.nodes()) # set de nodos de G
    nodos_ess_G = nodos_G.intersection(set(ess)) # nodos esenciales de G (como interseccion entre nodos de G y esenciales)
    nodos_no_ess_G = nodos_G.difference(set(ess)) # nodos no esenciales de G (como diferencia entre nodos de G y esenciales)
    
    # Agrego el atributo correspondiente a cada nodo
    G.add_nodes_from(nodos_ess_G, essential=True)
    G.add_nodes_from(nodos_no_ess_G, essential=False)

    #Pasamos a listas
    nodos_ess_G=list(nodos_ess_G)
    nodos_no_ess_G=list(nodos_no_ess_G)
	
    return (G,nodos_ess_G,nodos_no_ess_G)

#--------------------------------------------------------------------------------------
def frac_ess(G):
    
# Toma un grafo G con atributo essential y devuelve dos vectores, k y nodos_frac
    grados_dict = dict(G.degree())
    ess_dict = nx.get_node_attributes(G,'essential')

    k_lista = list(grados_dict.values()) # lista de grados de nodos en orden
    k = np.unique(k_lista) # vector de grado de nodos sin repetir

    L = len(k)
    nodos_ess = np.zeros(L)
    nodos_total = np.zeros(L)
    nodos_frac = np.zeros(L)

    for i,grado in enumerate(k):
        nodos_total[i] = k_lista.count(grado)
		# cuenta cuantas veces aparece cada grado en k_lista
		
    for proteina in ess_dict:
        if ess_dict[proteina] == True:
            i = np.where(k == grados_dict[proteina])
            nodos_ess[i]=nodos_ess[i]+1

    nodos_frac = nodos_ess / nodos_total
    nodos_frac=list(nodos_frac)
	
	# Devuelve dos listas: el vector k con grados (eje x) y el vector nodos_frac con la fraccion de nodos esenciales (eje y)
    k=list(k)
    nodos_frac=list(nodos_frac)
	# Podria agregarse que tambien devuelva nodos totales y nodos esenciales
	
    return(k,nodos_frac)
#-------------------------------------------------------------------------------
def log_Pe(G):

# Toma un grafo G con atributo essential y devuelve dos vectores, k y ln(1-Pe)
    grados_dict = dict(G.degree())
    ess_dict = nx.get_node_attributes(G,'essential')

    k_lista = list(grados_dict.values()) # lista de grados de nodos en orden
    k = np.unique(k_lista) # vector de grado de nodos sin repetir

    L = len(k)
    nodos_ess = np.zeros(L)
    nodos_total = np.zeros(L)
    nodos_frac = np.zeros(L)

    for i,grado in enumerate(k):
        nodos_total[i] = k_lista.count(grado)
        # cuenta cuantas veces aparece cada grado en k_lista
    for proteina in ess_dict:
        if ess_dict[proteina] == True:
            i = np.where(k == grados_dict[proteina])
            nodos_ess[i] += 1

    nodos_ess=np.array(nodos_ess, dtype=float) #pasamos a float 
    nodos_frac = nodos_ess / nodos_total
    
    menos_Pe = np.log(1 - nodos_frac)
    
	# Devuelve dos listas : el vector k con grados (eje x) y el vector log(1-Pe) (eje y)
    menos_Pe=list(menos_Pe)
    k=list(k)
        
    return (k, menos_Pe)

#-------------------------------------------------------------------------------
