# This file contains helper functions for helping the foreach package run
# smoothly by including all the functions loaded in the environment into
# each worker.

ListFunctions <- function() {
    # Return the names of all the functions in the global environment.
    #
    # Returns:
    #   A character vector with the names of all the functions in the global environment.
    inlist <- ls(.GlobalEnv)
    type <- 'closure'
    typelist <- sapply(sapply(inlist,get),typeof)
    
    return(names(typelist[typelist==type]))
}