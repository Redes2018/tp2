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
    #manteniendo el grado de cada nodo.
    #Estrategia que vimos en clase de Cherno de redes random:
    #1) Realizamos una copia del grafo G sin enlaces.
    #2) Vamos tomando pares de nodos al azar y creamos enlaces manteniendo el grado de cada nodo hasta agotar el cupo.
    #Este proceso nos va a dejar con tres tipos de enlaces entre nodos: e. simples, e. autoloops y e. multienlace. Estos dos últimos son enlaces problemáticos
    #y hay que eliminarlos.
    #3)Eliminamos los autoloops.
    #4)Eliminamos los multienlaces
    #5)Corroboramos que no queden ni autoloops ni multienlaces.
    #6)Ultimo chequeo para corroborar que no haya cambiado el grado de cada nodo.
    
    nodos=list(G.nodes)
    enlaces=list(G.edges)
    grados_dict = dict(G.degree())
    k_nodo= list(grados_dict.values()) # lista de grados de cada nodo en nodos(ordenados)
    print('grado G')
    print('enlaces originales en G sin autoloops: {}'.format(len(enlaces))) #hemos removido autoloops en el programa principal

    #Ahora nos quedamos con nodos de k distinto de 0 de esa forma mantengo el k=0 para los nodos aislados:
    index_nonzerok=[i for i in range(0,len(k_nodo)) if k_nodo[i]!=0] #buscamos los lugares de k_control donde no hayan quedado zeros
    k_nodo=[k_nodo[i] for i in range(0,len(k_nodo)) if k_nodo[i]!=0]
    k_nodo_antes=k_nodo
    nodos=[nodos[index_nonzerok[i]] for i in range(0,len(index_nonzerok))]
    enlaces=list(G.edges)
    
    #1) Creo un multigraph D que acepte multiedges
    D = nx.MultiGraph()

    #Agrego nodos:
    D.add_nodes_from(nodos) #Nota D solo va a tener los nodos de G que se hallan conectdos, no posee los nodos de G que se encuentran aislados
    
    #Inicializo k_control y nodos_new:
    k_control=np.array(k_nodo) #cuando creo un enlace entre nodoi y nodoj se le restara un 1 a los lugares i y j de k_control
    nodos_new=nodos

    #2)Agregamos enlaces de forma aleatoria al grafo D, manteniendo controlado que el cupo de cada nodo no puede exceder su grado.
    while(len(nodos_new)>0): 
        #Elijo uno random de pairs
        pair= np.random.choice(nodos_new,(1,2),replace=True)[0] #al poner replace True permito sacar dos numeros iguales.(eso va a crear un autoloop)
                   
        #Actualizamos variable de control: k_control
        if pair[0] == pair[1]:
            if k_control[nodos.index(pair[0])]>1:
                k_control[nodos.index(pair[0])]=k_control[nodos.index(pair[0])]-2 #solo actualizamos ese y le restamos un 2 ya que se creo un autoloop
                #creamos el autoloop
                D.add_edge(pair[0], pair[1])
                #no es posible crear el autoloop si habia un 1 en k_nodos, como mínimo necesito tener un 2 o mayor en el vector de grado de ese nodo.
        else:
            #creamos el enlace el cual no es un autoloop
            D.add_edge(pair[0], pair[1])
            k_control[nodos.index(pair[0])]=k_control[nodos.index(pair[0])]-1 #actualizamos k_control en la pos i
            k_control[nodos.index(pair[1])]=k_control[nodos.index(pair[1])]-1 #actualizamos k_control en la pos j

    
        #Actualizamos variable de control: nodos_new
        if k_control[nodos.index(pair[0])]==0 or k_control[nodos.index(pair[1])]==0: #solo actualizo k_control cuando alguno de los valores llega a cero
            index_nonzerok=[i for i in range(0,len(k_control)) if k_control[i]>0] #buscamos los lugares de k_control donde hayan elementos dinstintos a cero
            index_equalzero=[i for i in range(0,len(k_control)) if k_control[i]==0]#buscamos los lugares de k_control donde hayan elementos igual a cero
            nodos_new=[nodos[index_nonzerok[i]] for i in range(0,len(index_nonzerok))] #actualizamos la lista de nodos asi no volvemos a tomar nodos que ya recableamos por completo o sea aquellos que alcanzaron k_control[i]=0   
      

    print('grafico D')
    enlaces=list(D.edges())
    #Enlaces problemáticos:
    #Selfloops:
    '''
    print('autoloops inicial: {}'.format(len(list(D.nodes_with_selfloops()))))
    autoloops=list(D.nodes_with_selfloops())
    '''
    print('autoloops inicial: {}'.format(len(list(D.selfloop_edges()))))
    autoloops=list(D.selfloop_edges())
    #print(autoloops)
    
    #Multiples y Simples:
    enlaces_multiples=[]
    enlaces_simples=[]
    
    for i in range(0,len(nodos)):
        for j in range(i+1,len(nodos)):
            if(D.number_of_edges(nodos[i],nodos[j]))>1:
                for k in range(0,D.number_of_edges(nodos[i],nodos[j])-1):#en este for agregamos al vector enlaces_multiples tantos enlaces como multiplicidd tenga el mismo menos 1(porque si hay 3 enlaces entre dos nodos solo hay que sacar 2 de ellos )                    
                    enlaces_multiples.append([nodos[i],nodos[j]]) #agrego multiplicidad -1 de enlaces muliples
                enlaces_simples.append([nodos[i],nodos[j]]) #agrego uno simple
            elif (D.number_of_edges(nodos[i],nodos[j]))==1:
                enlaces_simples.append([nodos[i],nodos[j]])
    
    print('multiples inicial: {}'.format(len(enlaces_multiples)))
    print('simples inicial: {}'.format(len(enlaces_simples)))
    
    #Comparamos grados en esta etapa intermedia si queremos:
    grados_dict = dict(D.degree())
    k_nodo_despues= list(grados_dict.values())

    #Hasta acá el programa lo que hizo fue reconectar las puntas conservando el constraint de los grado de los nodos.
    #El problema de esto es que aparecieron autoloops en un mismo nodo y multienlaces entre nodos distintos.
    #Estos enlaces los vamos a llamar enlaces problemáticos.

    #Por ultimo hay que eliminar estos enlaces que son problemáticos:
    numero_autoloops=len(autoloops)
    numero_enlaces_multiples=len(enlaces_multiples)

    #3) Eliminemos autoloops primero:
    print('Recableando autoloops...')
    while(numero_autoloops >0):
        for al in autoloops:
            idx = np.random.choice(len(enlaces_simples),1)[0] #elijo un enlace dentro de los simples o sea no problematicos
            enlace_elegido=enlaces_simples[idx]
            if (enlace_elegido[0]!=al[0]) & (enlace_elegido[1]!=al[0]): #acepto ese enlace al azar si ninguno es el nodo donde esta el autoloop. esto evita crear un nuevo autoloop en el ismo nodo.
                #Hago el swap:
                #Creo dos nuevos
                D.add_edge(al[0],enlace_elegido[0])
                D.add_edge(al[0],enlace_elegido[1])
                #Elimino dos
                D.remove_edge(al[0],al[0])
                D.remove_edge(enlace_elegido[0],enlace_elegido[1])
                #Recalculamos autoloops y numero_autoloops en cada paso:
                autoloops=list(D.selfloop_edges())
                numero_autoloops=len(autoloops)
                #Tengo que actualizar enlaces simples:
                enlaces_simples.remove([enlace_elegido[0],enlace_elegido[1]])
                enlaces_simples.append([al[0],enlace_elegido[0]])
                enlaces_simples.append([al[0],enlace_elegido[1]])
    print('autoloops intermedio: {}'.format(len(list(D.nodes_with_selfloops()))))

    
    #Actualizamos los enlaces multiples:
    #Multiples y Simples:
    enlaces_multiples=[]
    enlaces_simples=[]
  
    for i in range(0,len(nodos)):
        for j in range(i+1,len(nodos)):
            if(D.number_of_edges(nodos[i],nodos[j]))>1:
                for k in range(0,D.number_of_edges(nodos[i],nodos[j])-1):#en este for agregamos al vector enlaces_multiples tantos enlaces como multiplicidd tenga el mismo menos 1(porque si hay 3 enlaces entre dos nodos solo hay que sacar 2 de ellos ).
                    enlaces_multiples.append([nodos[i],nodos[j]])#agrego multiplicidad -1 de enlaces muliples
                enlaces_simples.append([nodos[i],nodos[j]]) #agrego uno simple
            elif (D.number_of_edges(nodos[i],nodos[j]))==1:
                enlaces_simples.append([nodos[i],nodos[j]])
    print('multiples intermedio: {}'.format(len(enlaces_multiples)))
    print('simples intermedio: {}'.format(len(enlaces_simples)))

    
    #4) Eliminamos los enlaces multiples:
    numero_enlaces_multiples=len(enlaces_multiples)
    print('Recableando multiples...')
    while(numero_enlaces_multiples >0):
        for em in enlaces_multiples:
            idx = np.random.choice(len(enlaces_simples),1)[0] #elijo un enlace dentro de los simples o sea no problematicos
            enlace_elegido=enlaces_simples[idx]
            loscuatronodos=[em[0],em[1],enlace_elegido[0],enlace_elegido[1]]
            A = nx.to_pandas_adjacency(D)
            a1=A[em[0]][enlace_elegido[0]]
            a2=A[em[0]][enlace_elegido[1]]
            a3=A[em[1]][enlace_elegido[0]]
            a4=A[em[1]][enlace_elegido[1]]
            adjacencynumber=a1+a2+a3+a4
            #A continuación solo recableamos si los 4 nodos son distintos sino no, porque puedo vovler a crear un autoloop y  ademas...
            #solo recableamos si son enlaces adyacentes sino no, esto evita que se vuelvan a formar mutienlaces.
            controlnumber=adjacencynumber + len(np.unique(loscuatronodos))
            if (controlnumber==4):
                #Hago el swap:
                #Creo dos nuevos
                D.add_edge(em[0],enlace_elegido[0])
                D.add_edge(em[1],enlace_elegido[1])
                #Elimino dos
                D.remove_edge(em[0],em[1])
                D.remove_edge(enlace_elegido[0],enlace_elegido[1])
                #Tengo que actualizar enlaces simples:
                enlaces_simples.remove([enlace_elegido[0],enlace_elegido[1]])
                enlaces_simples.append([em[0],enlace_elegido[0]])
                enlaces_simples.append([em[1],enlace_elegido[1]])
                #Tengo que actualizar enlaces_multiples
                enlaces_multiples.remove([em[0],em[1]])
                numero_enlaces_multiples=len(enlaces_multiples)

    #5)Nos fijamos que no hallan quedado autoloops:
    autoloops=list(D.nodes_with_selfloops())
    print('autoloops final: {}'.format(len(list(D.nodes_with_selfloops()))))
    
    #Por ultimo me fijo los multiples al final:(deberia ser cero)
    enlaces_multiples=[]
    enlaces_simples=[]
    for i in range(0,len(nodos)):
        for j in range(i+1,len(nodos)):
            if(D.number_of_edges(nodos[i],nodos[j]))>1:
                for k in range(0,D.number_of_edges(nodos[i],nodos[j])-1):
                    enlaces_multiples.append([nodos[i],nodos[j]])#en este for agregamos al vector enlaces_multiples tantos enlaces como multiplicidd tenga el mismo.
                enlaces_simples.append([nodos[i],nodos[j]]) #agrego uno simple
            elif (D.number_of_edges(nodos[i],nodos[j]))==1:
                enlaces_simples.append([nodos[i],nodos[j]])
    print('multiples final: {}'.format(len(enlaces_multiples)))
    print('simples final: {}'.format(len(enlaces_simples)))   

    #6) Chequeo final para ver que se mantuvo el grado k de los nodos:
    grados_dict = dict(D.degree())
    k_nodo_despues= list(grados_dict.values())
    diferencia=np.array(k_nodo_despues)-np.array(k_nodo_antes)
    if (len(np.where(diferencia!=0)[0])==0):
        print('Rewiring exitoso')
    
    print('Listo')
    
    return(D)
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

#-----------------------------------------------------------------------------------

def fig3(G,ess,name):      
    
    #Creo copias superficiales de G porque solo me interesan sus nodos
    G_copy0 = G.copy()
    G_copy1 = G.copy()
    G_copy2 = G.copy()
    G_copy3 = G.copy()
    G_copy4 = G.copy()
    
    
    L = len(G) # tamaño del grafo (cantidad de nodos)
    nodos_elim = np.arange(L) # vector que va del 0 al L-1, cantidad de nodos que se va a ir eliminando (eje x)
    forn= nodos_elim/L #fracción de nodos eliminados
        
    lcc_dc= np.ones(L)#componente gigante obtenida mediante remoción de nodos por degree centrality (eje y)
    lcc_dc[0] = L # la primera componente es el tamaño original del grafo

    lcc_random= np.ones(L) #componente gigante obtenida mediante remocion de nodos aleatorios (eje y)
    lcc_random[0] = L # la primera componente es el tamaño original del grafo
    
    lcc_ess= np.ones(L) #componente gigante obtenida mediante remoción de nodos esenciales de mayor a menor "esencialidad" (eje y)
    lcc_ess[0] = L # la primera componente es el tamaño original del grafo
    
    
    #Le voy sacando nodos de forma aleatoria. Defino una lista de los nodos de G y la mezlco
    ls_r = list(G.nodes())
    random.shuffle(ls_r) #ls_r ya queda mezclada
    for i in range(L-1):
        G_copy0.remove_node(ls_r[i])
        lcc_random[i+1] = len(max(nx.connected_component_subgraphs(G_copy0),key=len))
        
    lcc_random = [ x/L for x in lcc_random]    

        
    #Le voy sacando nodos siguiendo de mayor a menor degree centrality y los ordeno
    ls_dc=nx.degree_centrality(G_copy1)
    ls_dc=sorted(ls_dc, key=ls_dc.__getitem__, reverse=True)  #Ordeno de menor a mayor, por eso revierto
    for i in range(L-1):
        G_copy1.remove_node(ls_dc[i])
        lcc_dc[i+1] = len(max(nx.connected_component_subgraphs(G_copy1),key=len))
    
    lcc_dc = [ x/L for x in lcc_dc]    
    
    #Le saco todas las esenciales de una
    from funciones import esenciales
    (G,ls_ess,lista_no_es) = esenciales(G_copy2,ess) #uso esta funcion para obtener las lista de las ess en G
    G_copy2.remove_nodes_from(ls_ess) #Le saco todos los esenciales
    lcc_ess = len(max(nx.connected_component_subgraphs(G_copy2),key=len))/L
    x=(len(ls_ess))/L  #defino mi x como la cant de nodos ess sacados/total
    
    from funciones import ec
    (forn_ec,lcc_ec) = ec(G_copy3, ess)
    
    from funciones import spbc
    (forn_spbc, lcc_spbc) = spbc(G_copy4, ess)
    


    
    #Guardamos los datos:
    output={}
    output[name+'fron_random']=forn
    output[name+'lcc_random']=lcc_random
    
    output[name+'forn_dc']=forn
    output[name+'lcc_dc']=lcc_dc
    
    output[name+'lcc_ess']=lcc_ess
    output[name+'forn_ess']= x
    
    output[name+'fron_ec']=forn_ec
    output[name+'lcc_ec']=lcc_ec
    
    output[name+'fron_spbc']=forn_spbc
    output[name+'lcc_spbc']=lcc_spbc
  
    
    df= pd.DataFrame()
    df['Date'] = output.keys()
    df['DateValue'] = output.values()
    df.to_csv(name+'_datos.txt', sep='\t')
    
    plt.figure()

    plt.plot(forn,lcc_random,'--r',label='Random')
    plt.plot(forn,lcc_dc,'--m',label='DC')
    plt.plot(x,lcc_ess,'b.',markersize="10",label='Ess')
    plt.plot(forn_ec,lcc_ec,'--c',label='EC')
    plt.plot(forn_spbc,lcc_spbc,'--k',label='SPBC')

    plt.xlabel('Cantidad de nodos eliminados')
    plt.ylabel('Tamano de la componente gigante')
    plt.title('Eliminacion de nodos segun distintas estrategias')
    plt.legend()
    plt.title(name)
    plt.savefig(name+'_grafico.png')
    plt.show()
    
    return(lcc_dc, lcc_random, lcc_ess, lcc_spbc, lcc_ec, forn, x)
           
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

#-------------------------------------------------------------------------------
def tabla3(G,name):   
    
    from funciones import esenciales
    from funciones import frac_ess

    E = G.copy()
    NE = G.copy()
    G1 = G.copy()
    L = len(G)

    (G1,ls_ess1, ls_no_ess1)=esenciales(G1,ess)

    #Le saco todas las esenciales de una
    NE.remove_nodes_from(ls_ess1) #Le saco todos los esenciales
    LCC_NE=max(nx.connected_component_subgraphs(NE),key=len)
    frac_sacando_E = LCC_NE.number_of_nodes()/L 
    degrees_no_es = [val for (node, val) in NE.degree()]
    k_no_es=np.unique(degrees_no_es)

    #Ahora le voy a sacar a G
    
    #Veamos los grados de los nodos esenciales que saque en G1
    E.remove_nodes_from(ls_no_ess1) #Le saco todos los no esenciales
    degrees_es = [val for (node, val) in E.degree()]
    k_es=np.unique(degrees_es)

    #Quiero que armar una lista de grados no esenciales que sean parecidos a los grados esenciales
    #y después encontrar sus nodos en G
    #Parecidos es que el grado difiera en al menos 1 unidad
    nueva_ls=[] 
    for i in range(len(k_no_es)-1):
        for j in range(len(k_es)-1):
            if abs(k_no_es[i] -k_es[j])<=1:
                nueva_ls.append(k_no_es[i])

    nueva_ls=np.unique(nueva_ls)

    #Veamos cuales son los nodos no-esenciales que cumplen esa cantidad de grados
    nodos_NE = [node for (node, val) in NE.degree()]
    grados_NE = [val for (node, val) in NE.degree()]
    nodos_NE_para_sacar=[]
    for i in range(len(nueva_ls)-1):
        for j in range(len(nodos_NE)-1):
            if nueva_ls[i] == grados_NE[j]:
                nodos_NE_para_sacar.append(nodos_NE[j])

    #Habiendo seleccionado los nodos, mezclo la lista y se la saco a G. Y hago estadística
    valores = []
    for l in range(0,1000): #mil repeticiones
        G0 = G.copy()
        random.shuffle(nodos_NE_para_sacar)
        G0.remove_nodes_from(nodos_NE_para_sacar)
        LCC_E = max(nx.connected_component_subgraphs(G0),key=len)
        frac_sacando_NE = LCC_E.number_of_nodes()/L 
        valores.append(frac_sacando_NE)


    valor_medio = sum(valores)/len(valores)
    frac_sacando_NE_mean = valor_medio
    error = np.std(valor_medio)
      
    #Guardamos los datos:
    output={}
    output[name+'frac_sacando_E']=frac_sacando_E
    output[name+'frac_sacando_NE_mean']=frac_sacando_NE_mean
    output[name+'error']=error
  
    
    df= pd.DataFrame()
    df['Date'] = output.keys()
    df['DateValue'] = output.values()
    df.to_csv(name+'_Tabla3.txt', sep='\t')
    
    return(frac_sacando_E, frac_sacando_NE_mean, error)
