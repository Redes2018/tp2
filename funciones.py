import sys
sys.path.append('/usr/local/lib/python2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sc
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

    #Uso lo que calculamos para hacer la figura 1A.
    frac_ess_nodes=np.zeros(L)
    hub_cutoff=np.zeros(L)
    for i, grado in enumerate(k):
        frac_ess_nodes[i]=np.sum(nodos_ess[i:])/float(np.sum(nodos_total[i:]))
        hub_cutoff[i]=np.sum(nodos_total[i:])/float(G.number_of_nodes())

    #Coeficienets correlacion
    kendall_tau=sc.kendalltau(k,frac_ess_nodes)[0]
    k_pvalue=sc.kendalltau(k,frac_ess_nodes)[1]
    kendall=[kendall_tau,k_pvalue]

    spearman_ro=sc.spearmanr(k,frac_ess_nodes)[0]
    s_pvalue=sc.spearmanr(k,frac_ess_nodes)[1]
    spearman=[spearman_ro,s_pvalue]
        
    # Devuelve seis listas:
    k=list(k)
    nodos_frac=list(nodos_frac)
    hub_cutoff=list(hub_cutoff)
    frac_ess_nodes=list(frac_ess_nodes)
        
    # Nota: nodos_frac[i] es la fraccion de nodos esenciales que hay en el grado k[i]
    #       frac_ess_nodes[i] es la fraccion de nodos esenciales mayor o iguales al grado k[i] que es la de la figura 1A.
	
    return(k,nodos_frac,hub_cutoff,frac_ess_nodes,kendall,spearman)
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
def ec(G,ess):
    #Toma un grafo G y va eliminando nodos de acuerdo a la medida de centralidad
    #eigenvalues centrality(ec). Devuelve dos listas, la fraccion de nodos elimi
    #nadas y el tamano de la componente conexa mas grande
    
    #Usamos la funcion nx.eigenvector_centrality que ya nos calcula el autovector
    #de autovalor mas grande.
    nodes_total=G.number_of_nodes()
    centrality = nx.eigenvector_centrality(G) #diccionario de centralidades
    centrality_sorted=sorted(centrality.items(), key=lambda x: x[1],reverse=True) #ordenados de mayor a menor

    nodes_to_remove=[centrality_sorted[i][0] for i in range(0,len(centrality_sorted))]

    #Eliminamos de mayor ec a menor
    remove_node_fraction=[]
    lcc=[]

    remove_node_fraction.append(0)
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    G0 = Gcc[0]
    lcc.append(G0.number_of_nodes()/float(nodes_total))
    
    for i in range(0,nodes_total-1):
        G.remove_node(nodes_to_remove[i])
        remove_node_fraction.append((i+1)/float(nodes_total))
        Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
        G0 = Gcc[0]
        lcc.append(G0.number_of_nodes()/float(nodes_total))
        
    #Outputs
    forn=list(remove_node_fraction) #fraction of removed nodes
    lcc_ec=list(lcc)                #largestconnectedcomponent

    return (forn,lcc_ec)

#-------------------------------------------------------------------------------
def spbc(G,ess):
    #Toma un grafo G y va eliminando nodos de acuerdo a la medida de centralidad
    #shortest path betweenness centrality(spbc). Devuelve dos listas, la fraccion de nodos elimi
    #nadas y el tamano de la componente conexa mas grande
    
    #Usamos la funcion nx.eigenvector_centrality que ya nos calcula el autovector
    #de autovalor mas grande.
    nodes_total=G.number_of_nodes()
    centrality = nx.betweenness_centrality(G,normalized=True) #diccionario de centralidades
    centrality_sorted=sorted(centrality.items(), key=lambda x: x[1],reverse=True) #ordenados de mayor a menor

    nodes_to_remove=[centrality_sorted[i][0] for i in range(0,len(centrality_sorted))]

    #Eliminamos de mayor spbc a menor
    remove_node_fraction=[]
    lcc=[]

    remove_node_fraction.append(0)
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
    G0 = Gcc[0]
    lcc.append(G0.number_of_nodes()/float(nodes_total))
    
    for i in range(0,nodes_total-1):
        G.remove_node(nodes_to_remove[i])
        remove_node_fraction.append((i+1)/float(nodes_total))
        Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
        G0 = Gcc[0]
        lcc.append(G0.number_of_nodes()/float(nodes_total))
        
    #Outputs
    forn=list(remove_node_fraction) #fraction of removed nodes
    lcc_spbc=list(lcc)                #largestconnectedcomponent

    return (forn,lcc_spbc)
#-----------------------------------------------------------------------------------
def G_simulation(G,alfa,beta):
    #Toma un grafo G y dos probabilidades alfa y beta. Asignamos la propiedad
    #esencial a los enlaces con probabilidad alfa. Asigna la probabilidad
    #esencial a los nodos con probabilidad beta. Devuelve un grafo G_model con
    #la propiedad esencial distribiuda en enlaces y nodos segun ese proceso.
    nodos=list(G.nodes())
    enlaces=list(G.edges())

    random_ess_enlace=np.random.choice(['True','False'],p=[alfa, 1-alfa], size=(1,len(enlaces)))[0]
    random_ess_nodo=np.random.choice(['True','False'],p=[beta, 1-beta], size=(1,len(nodos)))[0]

    #Creamos un diccionario de enlaces con la prop essential
    ess_enlace_dic={}
    for e, enlace in enumerate(enlaces):
        ess_enlace_dic[enlace]=random_ess_enlace[e]

    #Creamos un diccionario de nodos con la prop essential
    ess_nodo_dic={}
    for n, nodo in enumerate(nodos):
        ess_nodo_dic[nodo]=random_ess_nodo[n]
    
    #Asignamos la propiedad esencial a enlaces
    nx.set_edge_attributes(G,ess_enlace_dic,'essential')

    #Asignamos la propiedad esencial a nodos
    nx.set_node_attributes(G,ess_nodo_dic,'essential')
    
    return G
#-----------------------------------------------------------------------------------
def rewiring(G):
    #Funcion de grafo G que toma un grafo y realiza un recableado
    #manteniendo el grado de cada nodo y devuelve ibeps la cantidad
    #de enlaces entre nodos esenciales luego de recablear.
    
    nodos=list(G.nodes)
    enlaces=list(G.edges)
    grados_dict = dict(G.degree())
    k_nodo= list(grados_dict.values()) # lista de grados de cada nodo en nodos(ordenados)

    #Eliminamos todos los enlaces:
    G.remove_edges_from(enlaces)
    
    #Inicializo k_control:
    k_control=np.array(k_nodo) #cuando creo un enlace entre nodoi y nodoj se le resta un 1 a los lugares i y j de k_control
    nodos_new=nodos
    while(len(nodos_new)>1):
        pair=list(np.random.choice(nodos_new, size=(1,2), replace=False)[0])#elegimos dos nodos al azar con reposicion
        if pair[0]!=pair[1]: #evitamos crear autoloops
            #Creamos el enlace
            G.add_edge(pair[0], pair[1])
            #Actualizamos variables de control: k_control y nodos_new
            k_control[nodos.index(pair[0])]=k_control[nodos.index(pair[0])]-1 #actualizamos k_control en la pos i
            k_control[nodos.index(pair[1])]=k_control[nodos.index(pair[1])]-1 #actualizamos k_control en la pos j
            index_nonzerok=[i for i in range(0,len(k_control)) if k_control[i]!=0] #buscamos los lugares de k_control donde no hayan quedado zeros
            nodos_new=[nodos[index_nonzerok[i]] for i in range(0,len(index_nonzerok)-1)] #actualizamos la lista de nodos asi no volvemos a tomar nodos que ya recableamos por completo
            
    #Contamos ibeps
    ibeps=0
    ess_dict = nx.get_node_attributes(G,'essential')
    enlaces=list(G.edges())
    
    for enlace in enlaces:
        if ess_dict[enlace[0]] == True & ess_dict[enlace[1]]==True:
            ibeps=ibeps+1
            
    return(G,ibeps)
#-----------------------------------------------------------------------------------
def alfa_He(G,N,name):
    #Recibe un grafo G y un entero N.
    #Calculamos los ibeps_real del grafo G original.
    #Realizamos N rewirings del grafo original y en cada uno
    #contamos al final la cantidad de ibeps_simul(enlaces entre
    #nodos esenciales).
    #Hacemos histograma
    #Calculamos el valor de alfa. Devolvemos alfa y error_alfa
    #y m que es la media del histograma de ibeps_simul.
    
    numero_de_rewirings=N
    ibeps_simul=[]

    #Contamos los ibeps en el grafo original
    ibeps_real=0
    ess_dict = nx.get_node_attributes(G,'essential')
    enlaces=list(G.edges())
    
    for enlace in enlaces:
        if ess_dict[enlace[0]] == True & ess_dict[enlace[1]]==True:
            ibeps_real=ibeps_real+1
            
    for i in range(0,numero_de_rewirings):
        print('rewiring nro {}'.format(i))
        D=G.copy()
        ibeps_simul.append(rewiring(D)[1])

    #Histograma
    plt.figure(1)
    plt.hist(ibeps_simul,color = 'lightseagreen',edgecolor='darkcyan', linewidth=1.2 ,label='Red Simulada',normed=1)
    plt.axvline(x=ibeps_real,color='r',label='Red Real')
    plt.xlabel('$Numero$ $de$ $IBEPS$')
    plt.ylabel('$Frecuencia$')
    plt.legend(loc='upper center')
    plt.title(name)
    plt.savefig(name+'_ibeps_grafo.png')
    plt.show()
    plt.close(1)

    #Estimacion de alfa:
    m=np.mean(ibeps_simul)
    alfa=(ibeps_real-m)/len(enlaces)
    error_alfa=np.std(ibeps_simul)/len(enlaces)

    #Guardamos los datos:
    output={}
    output['ibeps']=ibeps_simul
    output['alfa']=alfa
    output['error_alfa']=error_alfa
    output['m']=m
    
    df= pd.DataFrame()
    df['Date'] = output.keys()
    df['DateValue'] = output.values()
    df.to_csv(name+'_ibeps_data.txt', sep='\t')

    return(alfa,error_alfa,m)                               
#-----------------------------------------------------------------------------------
def beta_He(G,m,N,name):
    
    #1)Contamos nodos esenciales
    nodos=list(G.nodes())
    ess_dict = nx.get_node_attributes(G,'essential')
    nodos_ess_real=[nodos[i] for i in range(0,len(nodos)) if ess_dict[nodos[i]]==True]
    numero_nodos_ess_real=len(nodos_ess_real)
    
    #Contamos los ibeps en el grafo original
    ibeps_real=0
    enlaces=list(G.edges())
    
    for enlace in enlaces:
        if ess_dict[enlace[0]] == True & ess_dict[enlace[1]]==True:
            ibeps_real=ibeps_real+1
    #print(ibeps_real)
    #print(m)

    numero_nodos_beta=[]
    numero_de_iteraciones=N
    for i in range(0,numero_de_iteraciones):
        print('iteracion {}'.format(i))
        #2)Hacemos una copia del grafo y borramos la esencialidad de los nodos:
        D=G.copy()
        ess_nodo_dic={}
        for n, nodo in enumerate(nodos):
            ess_nodo_dic[nodo]=''
        nx.set_node_attributes(D,ess_nodo_dic,'essential')

        #3)Asignamos esencialidad a (enlaces_PPI=ibeps-m) enlaces
        enlaces=list(D.edges())
        enlaces_PPI=ibeps_real-m
        #print(enlaces_PPI)
        enlaces_PPI_control=enlaces_PPI
        enlaces_new=enlaces
        while (enlaces_PPI_control>0):
            idx = np.random.choice(len(enlaces_new),1)[0]
            enlace_elegido=enlaces_new[idx]
            #Asignamos la propiedad esencial a ese enlace
            D.add_edge(enlace_elegido[0],enlace_elegido[1],essential=True)
            #Asignamos la propiedad esencial a cada nodo por estar en un enlace esencial:
            D.add_node(enlace_elegido[0],essential=True)
            D.add_node(enlace_elegido[1],essential=True)
            #Actualizamos variable de control
            enlaces_new.remove(enlace_elegido)
            enlaces_PPI_control=enlaces_PPI_control-1
        
        #4)Ahora marcamos nodos de forma random hasta llegar a tener en la red
        #numero_nodos_ess_real

        #Contamos cuantos esenciales ya tengo y cuantos faltan completar
        #nodos=list(D.nodes())
        ess_dict = nx.get_node_attributes(D,'essential')
        nodos_ess_tengo=[nodos[k] for k in range(0,len(nodos)) if ess_dict[nodos[k]]==True]
        numero_nodos_ess_tengo=len(nodos_ess_tengo)
        #Nota: nodos_ess_simul va cambiando en cada corrida, lo cual creo esta bien
        #ya que al asignar enlaces esenciales, la esencialidad de los nodos va cambiando
        #Si asigno un dos enlaces a una misma proteina entonces la esencialidad de esa proteina sigue siendo un
        #lo cual es distinto si asigno dos enlaces esenciales no a la misma proteina, la cantidad de esenciales
        #sera mayor.
        
        numero_nodos_ess_quiero=numero_nodos_ess_real
        #print(numero_nodos_ess_quiero)
        numero_nodos_ess_acompletar=numero_nodos_ess_quiero-numero_nodos_ess_tengo
        #print(numero_nodos_ess_tengo)
        #print(numero_nodos_ess_acompletar)
        
        #Nota: puede pasar que cuando asigne enlacesPPI eso me dio una cantidad de
        #enlaces esenciales mas alta que la cantidad de enlaces de mi red real.
        #en este caso creo que no se puede aplicar el algoritmo para encontrar beta.
        #ya que el siguiente paso no tiene sentido. Lo que me parece es que si pasa
        #esto necesito de mas iteraciones para caer en una situacion en donde
        #el numero de enlaces esenciales por PPI sea menor que los esenciales en la red real.
        
        #Informamos del error si pasa esto:
        if numero_nodos_ess_quiero < numero_nodos_ess_tengo:
            print('Error: Asignamos enlaces PPI y el numero de nodos esenciales resulto mayor que el de la red real.')


        #Completamos esa cantidad:
        numero_nodos_ess_acompletar_control=numero_nodos_ess_acompletar
        nodos_new=list(D.nodes())
        numero_nodos_ess_tenia=numero_nodos_ess_tengo
        
        while (numero_nodos_ess_acompletar_control > 0):
            nodo_elegido=np.random.choice(nodos_new)#elegimos un nodo al azar
            if ess_dict[nodo_elegido]=='': #Nos fijamos si el elegido ya era esencial.''(vacío) es que no era esencial
                D.add_node(nodo_elegido,essential=True) #lo ponemos en escencial
                #Actualizamos las variables de control
                numero_nodos_ess_acompletar_control=numero_nodos_ess_acompletar_control-1 #restamos uno ya que si pase de false a true sumamos un escencial mas a mi red 
            #Actualizamos las variables de control
            nodos_new.remove(nodo_elegido)#no lo vuelvemos a elegir
            ess_dict = nx.get_node_attributes(D,'essential')
            nodos_ess_actual=[nodos[j] for j in range(0,len(nodos)) if ess_dict[nodos[j]]==True] #check
            #print(len(nodos_ess_actual)) #check

        #5)Calculamos el numero de nodos esenciales marcados por otros factores:
        nodos=list(D.nodes())
        nodos_ess_finales=[nodos[j] for j in range(0,len(nodos)) if ess_dict[nodos[j]]==True]
        numero_nodos_ess_finales=len(nodos_ess_finales)
        #print(numero_nodos_ess_finales) #check
        #print(numero_nodos_ess_tenia)   #check
        numero_nodos_ess_otros_factores=numero_nodos_ess_finales-numero_nodos_ess_tenia
        numero_nodos_beta.append(numero_nodos_ess_otros_factores)

    #Estimacion de beta:
    beta=np.mean(numero_nodos_beta)/numero_nodos_ess_real
    beta_error=np.std(numero_nodos_beta)/numero_nodos_ess_real

    #Guardamos los datos:
    output={}
    output['numero_nodos_ess_real']=numero_nodos_ess_real
    output['betas']=numero_nodos_beta
    output['beta_mean']=beta
    output['beta_error']=beta_error
  
    
    df= pd.DataFrame()
    df['Date'] = output.keys()
    df['DateValue'] = output.values()
    df.to_csv(name+'_ibeps_data_beta.txt', sep='\t')

    
    return(beta,beta_error)

#-------------------------------------------------------------------------------
def pairs(G):
	# Toma un grafo G (con atributo Esencial) y devuelve la cantidad de pares de nodos no adyacentes
	# y la cantidad de pares de nodos no adyacentes con al menos 3 vecinos en común
	ess_dict = nx.get_node_attributes(G,'essential')
	nodos_lista = list(G.nodes())

	A = nx.to_numpy_matrix(G) # matriz de adyacencia de G
	T = len(A) # tamaño de la matriz, debe ser igual que la cantidad de nodos, me aseguro:

	if G.number_of_nodes()!=T:
		print('El grafo y la matriz no son del mismo tamaño')

	for i in range(T):
		A[i,i]=0 # Como hay auto loops, pongo ceros en la diagonal

	A2 = A**2 # Creo la matriz de adyacencia al cuadrado
	# El lugar i,j de esta matriz me dice cuantos caminos de longitud 2 hay entre el nodo i y el j
	
	# Como quiero quedarme con pares de nodos que tengan al menos 3 vecinos en común,
	# busco que haya al menos 3 caminos de longitud 2 entre ellos
	I,J = np.where(A2 >= 3)
	# obtengo los indices de los lugares
	
	# Cuento cantidad de pares de nodos con 3 o más vecinos en común
	pares = 0
	pares_iguales = 0

	for i in range(len(I)):
		if I[i]!=J[i] and A[I[i],J[i]] == 0: # que no estén en la diagonal y que no sean 1os vecinos
			pares +=1
			nodo1 = nodos_lista[I[i]]
			nodo2 = nodos_lista[J[i]]
			if ess_dict[nodo1] == ess_dict[nodo2]:
				pares_iguales +=1
	pares = int(pares/2)
	pares_iguales = int(pares_iguales/2)
	
	return (pares, pares_iguales)