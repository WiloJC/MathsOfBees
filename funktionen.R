## Functions that describe the node and weight degree of empirical bipartite networks. 
# Last update: 21.03.2024

library(dplyr) 
library(stringr)


## =========================================================== ##
###### List of interaction matrices that evolve over time ######
## =========================================================== ##

addinf = function(dt, p, b, dOrdnung, datum, lili) { 
  # Take two list (plant and bee species) an the time when those species interact and produce a list of matrices where each matrix represent the time where the interactions where recorded
  # dt = Total number of time steps (total number of matrices)
  # p = Plant species names list
  # b = Bee species names list
  # dOrdung = Order of the matrices 
  # datum = Time when the plant and bee interact
  # lili = List of matrices that as input is empty but as output it has the interactions
  k = 1     # Add information to each matrix
  for( i in 1:dt ) {      
    for(j in k:length(b) ){
      if( dOrdnung[[1]][i] == datum[j] ) {
        lili[[i]][ p[j],b[j] ] = lili[[i]][ p[j],b[j] ] + 1
        # if( lili[[i]][ p[j],b[j] ] == dOrdnung[[2]][i]){
        #   k = j - dOrdnung + 1
        #   break
        # }
      }
    }
  }
  return(lili)
}

## =========================================================== ##
###### Nodes degree ######
## =========================================================== ##
gradod = function(matrice){
  # Compute the node degree in a bipartite network (Adjunted matrix). 
  # it produce two matrices, one  for the plants (rows) and one for the animals (columns). Each contains (1) the place of the node in the matrix  (identification); (2) the number of interactions; and(3) the number of links. 
    # matrice = Matrix that contains the interactions (it could contain the name of the species but it is not necessary). 
  d = dim(matrice)
  n = d[1] # Pflanzenartenzahl
  m = d[2] # Tierartenzahl
  
  # die Pflanzen (Rows)
  # Jede Art in der Matrix
  dp = vector( mode = "double", length = n )    
  dpa = matrix( data = 0, nrow = n, ncol = 2)
  nome = vector ( mode = "character", length = n )
  pn = rownames( matrice )
  kp = 0
  for( i in 1:n ) {
    dp[i] = sum( matrice[i,] )
    if(dp[i] != 0){
      kp = kp + 1
      dpa[kp,1] = i
      dpa[kp,2] = dp[i]
      nome[kp] = pn[i]
    }
  }
  # Nur Arten mit Wechselwirkungen.
  pf = matrix(data = 0, nrow = kp, ncol = 2)  
  rownames( pf ) = nome[1:kp]
  for( i in 1:kp ) {
    pf[i,1] = dpa[i,1]
    pf[i,2] = dpa[i,2]
  }

  # Tiere (das Tier) (Columns)
  db = vector( mode = "double", length = m )
  dba = matrix(data = 0, nrow = m, ncol = 2)
  nome1 = vector ( mode = "character", length = m )
  bn = colnames( matrice ) 
  kb = 0
  for( i in 1:m ) {
    db[i] = sum( matrice[,i] )
    if ( db[i] != 0){
      kb = kb + 1
      dba[kb,1] = i
      dba[kb,2] = db[i]
      nome1[kb] = bn[i]
    }
  }
  af = matrix(data = 0, nrow = kb, ncol = 2)
  rownames( af ) = nome1[1:kb]
  for( i in 1:kb ) {
    af[i,1] = dba[i,1]
    af[i,2] = dba[i,2]
  }
  
  colnames( af ) = c( "id", "#int" )
  colnames( pf ) = c( "id", "#int" )
  pdegree = pf
  adegree = af
  #return( list(pdegree = pf, adegree = af) )
  return( list( plant = pdegree, animal = adegree ) )
}

## =========================================================== ##
###### Most connected nodes ######
## =========================================================== ##

popular = function( lm, tl = 1 ){ # Compute the number of times that a node appears along the time series of matrices. Do not calculate the number of interactions that a node has. 
  # lm = List of matrices produced with the gradod function 
  # tl = taxonomy level: 1 = plants (rows), 2 = animals (columns)
  ni = c()  # List of node's identifications for all the matrices in the list of matrices produces by gradodd function. 
  n = 0  # total number of nodes in the lm list.
  for( i in 1:length(lm) ){
    for( j in 1: length( lm[[i]][[tl]][,1] ) ){
      n = n + 1
      ni[ n ] = lm[[i]][[tl]][j,1]
    }
  }
  wie = count( data.frame(ni), ni )
  wie = wie[ order( wie[[2]], decreasing = TRUE ), ]
  
  return(wie)
}

## =========================================================== ##
###### Node's weight dynamics ######
## =========================================================== ##

evol = function( node, lm, tl = 1 ){  
  # compute de node's degree dynamic
  # It produce a matrix (ndyn) that contains when the node has new links, ndyn[,1], the number new interaction in the node, ndyn[,2], and the total amount of links in the node at each time step ndyn[,3].
  # lm = List of matrices that contain the information provided by function gradod().
  # tl = Trophic level, 1 for plants (default) and 2 for animals. 
  # node = Number that identify the node along the list of matrices (could be a name but it was too slow).
  
  nt = 0    # nt = Number of times that the node has new links.
  dove = c()  # The position of the node in the gradod() matrix (It is for informational purposes only).
  quando = c()    # In which matrix of the lm list is the node.
  t = length(lm) # Total number of matrices in lm list (total number of time steps). 
  
  for( i in 1:t ){
    if( node %in% lm[[i]][[tl]][,1] ){
      nt = nt + 1
      dove[ nt ] =  match( node, lm[[i]][[tl]][,1] )
      quando[ nt ] = i
    }
  }
 
  ndyn = matrix( data = 0, nrow = nt, ncol = 5 ) # node dynamics
  ndyn[ ,1 ] = quando
  aux = 0
  aux1 = 0
  for( i in 1:nt ){
    ndyn[ i,2 ] = lm[[ quando[i] ]][[ tl ]][dove[i], 2 ] 
    aux = aux + ndyn[i,2]
    ndyn[ i,3 ] = aux
    aux1 = aux1 + 1
    ndyn[ i,4 ] = aux1
    ndyn[ i,5 ] = dove[i]
  }
  
  colnames( ndyn ) = c( "dt", "#int", "Tint", "Tlink", "spId" )
  
  return( ndyn )
}

## =========================================================== ##
###### Node's degree dynamics ######
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
## Functions used un the development of ndeg() function
# ndegtest = function( lima, tl = 1 ){
#   nt = length( lima ) # Total number of time steps. 
#   # nmat = matrix( data = 0, nrow = dim( lima[[1]] )[ 1 ], ncol = dim( lima[[1]] )[ 2 ] )
#   # print( "fi" )
#   
#   ## Copy list of imput matrices ##
#   nli = list() # New list of matrices with 1 and 0 as elements. 
#   for( i in 1:nt ){
#     rname = rownames( lima[[ i ]] )
#     cname = colnames( lima[[ i ]] )
#     nmat = matrix( data = 0, nrow = dim( lima[[ i ]] )[ 1 ], ncol = dim( lima[[ i ]] )[ 2 ] )
#     rownames( nmat ) = rname
#     colnames( nmat ) = cname
#     # print( i ) 
#     for( j in 1:dim( lima[[ i ]] )[ 1 ] ){
#       nmat[ j, which( lima[[ i ]][ j, ] >= 1 ) ] = 1
#     }
#     nli[[ i ]] = nmat    
#   }
# 
#   # return(evd)
#   return( nli )
# }
# 
# ndeg2 = function( lima, tl = 1 ){
#   nt = length( lima ) # Total number of time steps. 
#   nsp = dim( lima[[ 1 ]] )[ 1 ] # Number of species. 
#   nspi = dim( lima[[ 1 ]] )[ 2 ] # Number of species to interact with 
#   # nmat = matrix( data = 0, nrow = dim( lima[[1]] )[ 1 ], ncol = dim( lima[[1]] )[ 2 ] )
#   # print( "fi" )
#   
#   ## Copy list of imput matrices ##
#   nli = list() # New list of matrices with 1 and 0 as elements. 
#   for( i in 1:nt ){
#     rname = rownames( lima[[ i ]] )
#     cname = colnames( lima[[ i ]] )
#     nmat = matrix( data = 0, nrow = dim( lima[[ i ]] )[ 1 ], ncol = dim( lima[[ i ]] )[ 2 ] )
#     rownames( nmat ) = rname
#     colnames( nmat ) = cname
#     # print( i ) 
#     for( j in 1:dim( lima[[ i ]] )[ 1 ] ){
#       nmat[ j, which( lima[[ i ]][ j, ] >= 1 ) ] = 1
#     }
#     nli[[ i ]] = nmat    
#   }
#   
#   ## time evolution ##
#   
#   sev = vector( mode = "list", length = nsp )  # Species evolution
#   names( sev ) = rownames( lima[[ 1 ]] ) # name of nodes 
#   for( i in 1:nsp ) {
#     sev[[ i ]] =  matrix( data = 0, nrow = nt, ncol = nspi )
#     for( j in 1:nt ){
#       # sev[[ i ]][ j, ] = lima[[ j ]][ i, ] # To check
#       sev[[ i ]][ j, ] = nli[[ j ]][ i, ]
#     }
#   }
#   
#   ## Just links ## 
#   sed = sev
#   for( j in 1:nsp ) {
#     for( k in 1:nspi ){
#       x = which( sev[[ j ]][ , k ] >= 1 )
#       sed[[ j ]][ x , k ] = 0
#       sed[[ j ]][ x[ 1 ], k ] = sev[[ j ]][ x[ 1 ], k  ]
#     }
#   }
#   
#   return( list( sev, sed ) )
#   
# } 
  
## =========================================================== ##
###### Node's weight dynamics (method 2) ######
## =========================================================== ##  

nweight = function( liin, tl = 1, mn = 1){  # Compute the dynamical evolution of the strength (or weight) of nodes given a list of matrices.
  # Produce a list, evd, where each element represent the dynamical evolution of a node (plant/animal species). evd[i][,1] time when the interaction happehed, evd[[i]][,2] number of weight, evd[[i]][,3] cumulative weight. 
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
  
  
  ## time evolution ##
  sev = vector( mode = "list", length = nsp )  # Species evolution
  names( sev ) = rownames( lima[[ 1 ]] ) # name of nodes 
  for( i in 1:nsp ) {
    sev[[ i ]] =  matrix( data = 0, nrow = nt, ncol = nspi )
    for( j in 1:nt ){
      sev[[ i ]][ j, ] = lima[[ j ]][ i, ] # To check (to use with weights).
    }
  }
  
  # ### Test (Number of interactions after the establishment of links) ##
  # ## Number of interactions after the establishment of links ## 
  # sed = sev   # New list fo evolution matrix that represent weight after the start of the mutualistic interaction.
  # for( j in 1:nsp ) { # Delete interaction after (and before) the number of interactions, nm, required to have a link. 
  #   for( k in 1:nspi ){
  #     x = which( sev[[ j ]][ , k ] > 0 )
  #     sed[[ j ]][ x[ 1: ( length( x ) + 1 - mn  ) ] , k ] = 0
  #     # print("##########")
  #     # print( x )
  #     # print("xx xx")
  #     # print(  x[ 1: ( length( x ) + 1 - mn  ) ] )
  #   }
  # } Not working so far. 
  
  ## evolution of links ##
  evd = vector( mode = "list", length = nsp )  # List with the evolution of node's degree.
  nodeid = rownames( lima[[ 1 ]] ) # name of nodes
  names( evd ) = nodeid
  
  ds = vector( mode = "double", length = nt )    # weight of node i 
  for( i in 1: nsp ) {
    evd[[ i ]] =  matrix( data = 0, nrow = nt, ncol = 3 )  # Table with infotmation with node's dynamics
    for( j in 1: nt ){
      # ds[ j ] = sum( sed[[i]][ j, ] ) # Test (Number of interactions after the establishment of links) Not working so far
      ds[ j ] = sum( sev[[i]][ j, ] )
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
    # evd[[ i ]] = evd[[ i ]][ 1: length( it ), ] # change the size of information table
    evd[[ i ]] = matrix( evd[[ i ]][ 1: length( it ), ], nrow = length( it ), ncol = 3,
                         dimnames = list( seq( 1, length( it ) ), c("dt" , "#link", "Tlink" ) ) )
  }
  
  return( evd )
  
} 

## =========================================================== ##
###### Mean value of node's degree/weight i.e. <k> and <w> over time ######
## =========================================================== ##

meanv = function( listdw ){
  # It uses the list of matrices (listdw) generated by the ndeg() and nweight() functions. 
  # It produces a data frame with two columns (1) time step and (2) the mean value of the cumulative degree/weight of the network
  nspecies = length( listdw ) # Total number of species 
  
  ## Compute the time range. From time step 1 to the highest time step ##
  maxi = vector( mode = "double", length = nspecies ) # Vector of maximum time steps for each species
  for( i in 1:nspecies ){ maxi[i] = max( listdw[[ i ]][ , 1 ] ) } #
  mdt = max( maxi )
  ## C t t r ##
  
  mv = vector( mode = "double", length = mdt ) # Mean value of k_i per time step
  mv1 = vector( mode = "double", length = mdt )  # Cumulative mean value of k_i
  
  for( j in 1:mdt ){
    cont = 0 
    for( k in 1:nspecies ){
      if( j %in% listdw[[k]][,1] ){ 
        cont = cont + 1
        mv[ j ] = mv[ j ] + listdw[[k]][ which( listdw[[ k ]][,1] %in% j ) , 2 ]
        # print( paste( "dt", j, "sp", k, "place", which( listdw[[ k ]][,1] %in% j), "#int", mv[ j ], "#sp", cont, sep = "--"  ) )     # Test
      }
    }
    mv1[ j ] = mv[ j ] / cont
    # print( paste( j, mv[ j ], cont, mv1[ j ], sep = " " )  )    # Test
  }
  
  ausgang = data.frame( "dt" = seq( 1: mdt ), "E(k)" = mv1 )
  ausgang = na.omit( ausgang )
  for( i in 2: length( ausgang[ , 2] ) ){ ausgang[ i, 2 ] = ausgang[ i, 2 ] + ausgang[ i - 1, 2 ] }
  
  return( ausgang )
}
