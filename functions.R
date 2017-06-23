# wrap C's printf-like function
printf <- function(...)cat(sprintf(...))

# This fuction is responsible for generation of many different versions of an incomplete tree. 
# Args: -integer(nrep): number of replications to be expanded
#       -character(newickf): path to a newick file that holds the phylogenetic tree
#       -character(putf): path to a text file that hold a list of PUT MDCC pair per line
#       -character(outf): path to the output file
#       -integer(seed): a numeric value required to make the simulation reproductible
# Output: a file containing all the generated tree in newick format
expansion <- function(nrep=1, newickf="newick.tree",putf="put.list", outf="versions.tree", seed=sample(.Machine$integer.max,1)){
  
  if(file.access(newickf,4)){
    stop("The file \"",newickf,"\" doesn't exist or cannot be read.")
  }
  if(file.access(putf,4)){
    stop("The file \"",putf,"\" doesn't exist or cannot be read.")
  }
  if(!is.loaded("expansion")){
    dyn.load("lib/librsunplin.so")
  }
  ret <- .C("expansion",
            as.integer(nrep),
            newickf,
            putf,
            outf,
            as.integer(0),
            as.integer(0),
            as.integer(0),
            as.integer(0),
            as.integer(seed)
            )
  nspc <- ret[[5]]
  nput <- ret[[6]]
  tspc <- ret[[7]]
  tnod <- ret[[8]]
  printf("Number of species of original tree: %d\n",nspc)
  printf("Number of PUTs to be inserted: %d\n",nput)
  printf("Number of species after insertion: %d\n",tspc)
  printf("Number of nodes (including internals): %d\n",tnod)
  printf("Number of replications: %d\n",nrep)
  printf("Output file: \"%s\"\n",outf)
}

# This fuction is responsible for generation of many different versions of an incomplete tree and output 
# a file containing the patristic matrix for each one.
# Args: -integer(nrep): number of replications to be expanded
#       -character(newickf): path to a newick file that holds the phylogenetic tree
#       -character(putf): path to a text file that hold a list of PUT MDCC pair per line
#       -character(outf): path to the output file
#       -integer(seed): a numeric value required to make the simulation reproductible
# Output: a file containing all the generated matrices in CSV format
patrices <- function(nrep=1, newickf="newick.tree",putf="put.list", outf="matrices.csv", seed=sample(.Machine$integer.max,1)){
  
  if(file.access(newickf,4)){
    stop("The file \"",newickf,"\" doesn't exist or cannot be read.")
  }
  if(file.access(putf,4)){
    stop("The file \"",putf,"\" doesn't exist or cannot be read.")
  }
  if(!is.loaded("chained_patrices")){
    dyn.load("lib/librsunplin.so")
  }
  ret <- .C("chained_patrices",
            as.integer(nrep),
            newickf,
            putf,
            outf,
            as.integer(0),
            as.integer(0),
            as.integer(0),
            as.integer(0),
            as.integer(seed)
  )
  nspc <- ret[[5]]
  nput <- ret[[6]]
  tspc <- ret[[7]]
  tnod <- ret[[8]]
  printf("Number of species of original tree: %d\n",nspc)
  printf("Number of PUTs to be inserted: %d\n",nput)
  printf("Number of species after insertion: %d\n",tspc)
  printf("Number of nodes (including internals): %d\n",tnod)
  printf("Number of replications: %d\n",nrep)
  printf("Output file: \"%s\"\n",outf)
}