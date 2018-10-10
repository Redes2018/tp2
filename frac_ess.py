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
			nodos_ess[i] += 1
		

	nodos_frac = nodos_ess / nodos_total
    
	# Devuelve el vector k con grados (eje x) y el vector nodos_frac con la fraccion de nodos esenciales (eje y)
	# Podria agregarse que tambien devuelva nodos totales y nodos esenciales
	
    return (k, nodos_frac)
