import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import os
from funciones import *

carpeta = (os.getcwd()+'/tc02Data/')

#-------------------------------------------------------------------------------------------------------------------
#                        TP2: PROGRAMA PRINCIPAL
#-------------------------------------------------------------------------------------------------------------------
#1)Proteinas Escenciales:
#-------------------------------------------------------------------------------------------------------------------
# Creamos un dataframe con las proteinas esenciales, pero el primer elemento es parte del encabezado
# podria eliminarlo con header pero si agregamos lineas a header no podemos seleccionar columnas con usecols
data = pd.read_csv(carpeta+'Essential_ORFs_paperHe.txt', sep='\t', header=0,skipfooter=4,usecols=[1],engine='python')

# Como algunos nombres de proteinas terminan con "  " hay que borrar los espacios vacios
data['ORF_name'] = data['ORF_name'].map(lambda x: x.strip()) 


# Nos podemos quedar con la serie de pandas o puedo hacer una lista, eliminando el primer elemento que sobraba
ess = data["ORF_name"].tolist()
del ess[0]
# ess es una lista con los nombres de las proteinas esenciales


#--------------------------------------------------------------------------------------------------------------------
#2)Grafos:
#--------------------------------------------------------------------------------------------------------------------
#Creamos una variable grafos que contenga los cuatro grafos que vamos a analizar.
grafos = []
archivos = ['Y2H','AP-MS','LIT','LIT_Reguly']
for j,arc in enumerate(archivos):
    data = pd.read_csv(carpeta+'yeast_'+arc+'.txt', sep='\t', header=None)
    grafo = nx.Graph()
    for i in range(len(data)):
        grafo.add_edges_from([(data[0][i],data[1][i])])
    grafos.append(grafo)

#--------------------------------------------------------------------------------------------------------------------
#3)Usamos la funcion esenciales para asignar la propiedad esencial a los grafos:
#--------------------------------------------------------------------------------------------------------------------
for j,grafo in enumerate(grafos):
    grafos[j]=esenciales(grafos[j],ess)[0]
    

