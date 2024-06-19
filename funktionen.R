## Functions that describe the node and weight degree of empirical bipartite networks. 
# Author: WjC
# Last update: 21.03.2024

library(bipartite)

## =========================================================== ##
###### interaction matrices ######
## =========================================================== ##

datab = read_csv("bee_interactions.csv")

li = data.frame( higher = datab$genus.species, lower = datab$FLOWER_SPECIES, webID = datab$DATE )
liin= frame2webs( li ) # Interaction matrices

## =========================================================== ##
###### Empirical node degree dynamics ######
## =========================================================== ##

ndeg = function( liin, tl = 1, mn = 1 ){ # Compute the degree of the nodes that change over time given a list of matrices.
  # liin = Imput List of matrices with interactions 
  # tl = Trophic level: (1) plants, (2) pollinators. 
  # mn = Mutualistic number. Define the number of interactions to establish a mutualistic relationship (link) 
  
  nt = length( liin ) # Total number of time steps. 
  
  lima = liin #  Avoid change the original matrices
  if ( tl == 2 ){ # 
    for( i in 1: nt) {
      lima[[ i ]] = t( lima[[ i ]] )
    }
  }
  
  nsp = dim( lima[[ 1 ]] )[ 1 ] # Number of species. 
  nspi = dim( lima[[ 1 ]] )[ 2 ] # Number of species to interact with .
  
  
  ## Copy list of imput matrices ##
  nli = list() # New list of matrices with 1 and 0 as elements. 
  for( i in 1:nt ){
    rname = rownames( lima[[ i ]] )
    cname = colnames( lima[[ i ]] )
    nmat = matrix( data = 0, nrow = dim( lima[[ i ]] )[ 1 ], ncol = dim( lima[[ i ]] )[ 2 ] ) # Auxiliar matrix
    rownames( nmat ) = rname
    colnames( nmat ) = cname
    for( j in 1:dim( lima[[ i ]] )[ 1 ] ){
      nmat[ j, which( lima[[ i ]][ j, ] >= 1 ) ] = 1
    }
    nli[[ i ]] = nmat    
  }
  # return( nli ) # Pr?fen
  
  ## time evolution ##
  
  sev = vector( mode = "list", length = nsp )  # Species evolution
  names( sev ) = rownames( lima[[ 1 ]] ) # name of nodes 
  for( i in 1:nsp ) {
    sev[[ i ]] =  matrix( data = 0, nrow = nt, ncol = nspi )
    for( j in 1:nt ){
      # sev[[ i ]][ j, ] = lima[[ j ]][ i, ] # To check (to use with weights).
      sev[[ i ]][ j, ] = nli[[ j ]][ i, ]
    }
  }
  
  ## Just links ## 
  sed = sev   # New list fo evolution matrix that represent mutualistic interactions (Links).
  for( j in 1:nsp ) { # Delete interaction after (and before) the number of interactions, nm, required to have a link. 
    for( k in 1:nspi ){
      x = which( sev[[ j ]][ , k ] >= 1 )
      sed[[ j ]][ x , k ] = 0
      sed[[ j ]][ x[ mn ], k ] = sev[[ j ]][ x[ mn ], k  ]
    }
  }
  
  # return( list( sev, sed ) ) # Pr?fen
  
  ## evolution of links ##
  evd = vector( mode = "list", length = nsp )  # List with the evolution of node's degree.
  nodeid = rownames( nli[[ 1 ]] ) # name of nodes
  names( evd ) = nodeid
  
  ds = vector( mode = "double", length = nt )    # Degree of species 
  for( i in 1: nsp ) {
    evd[[ i ]] =  matrix( data = 0, nrow = nt, ncol = 3 ) # , dimnames = list( seq(1,nt), c("dt" , "#link", "Tlink" ) ) )  # Table with information with node's dynamics
    # print( evd[[ i ]] )    # tucheq
    for( j in 1: nt ){
      ds[ j ] = sum( sed[[i]][ j, ] )
      # print( ds[ j ] )  # tucheq
    }
    it = which( ds != 0 ) # interaction times 
    aux = 0
    for( k in 1:length( it ) ){
      # print(it)
      evd[[ i ]][ k, 1 ] = it[ k ]
      evd[[ i ]][ k, 2 ] = ds[ it[ k ] ]
      aux = aux + evd[[ i ]][ k, 2 ]
      evd[[ i ]][ k, 3 ] = aux
    }
    
    evd[[ i ]] = matrix( evd[[ i ]][ 1: length( it ), ], nrow = length( it ), ncol = 3,
                         dimnames = list( seq( 1, length( it ) ), c("dt" , "#link", "Tlink" ) ) ) # change the size of information table (it must be a matrix)
    
  }
  
  return( evd )
  
} 

## =========================================================== ##
###### Empirical node strength dynamics ######
## =========================================================== ##  

nweight = function( liin, tl = 1){  # Compute the dynamical evolution of the strength (or weight) of nodes given a list of matrices.
  # Produce a list, evd, where each element represent the dynamical evolution of a node (plant/animal species). evd[i][,1] time when the interaction happehed, evd[[i]][,2] number of weight, evd[[i]][,3] cumulative weight. 
  # liin = Imput List of matrices with interactions 
  # tl = Trophic level: (1) plants, (2) pollinators. 
  
  nt = length( liin ) # Total number of time steps. 
  
  lima = liin #  Avoid change the original matrices
  if ( tl == 2 ){ # 
    for( i in 1: nt) {
      lima[[ i ]] = t( lima[[ i ]] )
    }
  }
  
  nsp = dim( lima[[ 1 ]] )[ 1 ] # Number of species. 
  nspi = dim( lima[[ 1 ]] )[ 2 ] # Number of species to interact with .
  
  
  ## time evolution ##
  sev = vector( mode = "list", length = nsp )  # Species evolution
  names( sev ) = rownames( lima[[ 1 ]] ) # name of nodes 
  for( i in 1:nsp ) {
    sev[[ i ]] =  matrix( data = 0, nrow = nt, ncol = nspi )
    for( j in 1:nt ){
      sev[[ i ]][ j, ] = lima[[ j ]][ i, ] # To check (to use with weights).
    }
  }
  
  ## evolution of links ##
  evd = vector( mode = "list", length = nsp )  # List with the evolution of node's degree.
  nodeid = rownames( lima[[ 1 ]] ) # name of nodes
  names( evd ) = nodeid
  
  ds = vector( mode = "double", length = nt )    # weight of node i 
  for( i in 1: nsp ) {
    evd[[ i ]] =  matrix( data = 0, nrow = nt, ncol = 3 )  # Table with infotmation with node's dynamics
    for( j in 1: nt ){
      ds[ j ] = sum( sev[[i]][ j, ] )
    }
    it = which( ds != 0 ) # interaction times 
    aux = 0
    for( k in 1:length( it ) ){
      evd[[ i ]][ k, 1 ] = it[ k ]
      evd[[ i ]][ k, 2 ] = ds[ it[ k ] ]
      aux = aux + evd[[ i ]][ k, 2 ]
      evd[[ i ]][ k, 3 ] = aux
    }
    # evd[[ i ]] = evd[[ i ]][ 1: length( it ), ] # change the size of information table
    evd[[ i ]] = matrix( evd[[ i ]][ 1: length( it ), ], nrow = length( it ), ncol = 3,
                         dimnames = list( seq( 1, length( it ) ), c("dt" , "#link", "Tlink" ) ) )
  }
  
  return( evd )
  
} 
