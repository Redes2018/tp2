#Funciones para el tp2:

#-----------------------------------------------------
def overlap(G,F):

#Le pasamos dos grafos G y F y nos devuelve
#cantidad de enlaces en común nenlaces_comunes tipo int

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
#def nombre funcion(variablesin):
#Breve resumen

