def esenciales(G,ess):
# A partir de un grafo G y una lista de proteinas esenciales ess
# devuelve el grafo con el atributo essential (True o False) y una lista de las proteinas esenciales de G
    nodos_G = set(G.nodes()) # set de nodos de G
    nodos_ess_G = nodos_G.intersection(set(ess)) # nodos esenciales de G (como interseccion entre nodos de G y esenciales)
    nodos_no_ess_G = nodos_G.difference(set(ess)) # nodos no esenciales de G (como diferencia entre nodos de G y esenciales)
    
    # Agrego el atributo correspondiente a cada nodo
    G.add_nodes_from(nodos_ess_G, essential=True)
    G.add_nodes_from(nodos_no_ess_G, essential=False)
	
	return (G,nodos_ess_G)