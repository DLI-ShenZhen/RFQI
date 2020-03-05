#' registerPackage
#' 
#' register packages into cluster
registerPackage = function(cl){
  
  pkg.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))
  
  lapply(pkg.names, function(x){
    clusterEvalQ(cl, require(x))
  })
  
}

#' registerParentVars
#' 
#' register varlist of parent environment and current environment into cluster
registerParentVars = function(cl){
  
  this.env <- environment()
  
  if ( identical(this.env,globalenv())==FALSE ){
    # register current environment variables
    clusterExport(cl,
                  ls(all.names = TRUE, env=this.env),
                  envir = this.env)
    # search parent environment
    this.env = parent.env(environment())
    # register parent environment variables
    clusterExport(cl,
                  ls(all.names = TRUE, env=this.env),
                  envir = this.env)
  }
  
  # ensure register globaleEnv varlist
  if (identical(this.env,globalenv())==FALSE){
    clusterExport(cl,
                  ls(all.names=TRUE, env=globalenv()),
                  envir=globalenv())
  }
  
}




#' Define the hack
#' 
mclapply.hack <- function(X, FUN, mc.cores, ...) {
  
  cores = mc.cores
  cl <- makeCluster( min(mc.cores, availableCores()-1) )
  
  tryCatch({
    
    registerPackage(cl)
    registerParentVars(cl)
    
    ## Run the lapply in parallel 
    return( parLapply(cl, X=X, fun=FUN, ...) )
    
  }, finally = {        
    ## Stop the cluster
    stopCluster(cl)
  })
}

## Warn the user if they are using Windows
# if( Sys.info()[['sysname']] == 'Windows' ){
#   message(paste(
#     "\n", 
#     "   *** Microsoft Windows detected ***\n",
#     "   \n",
#     "   For technical reasons, the MS Windows version of mclapply()\n",
#     "   is implemented as a serial function instead of a parallel\n",
#     "   function.",
#     "   \n\n",
#     "   As a quick hack, we replace this serial version of mclapply()\n",
#     "   with a wrapper to parLapply() for this R session. Please see\n\n",
#     "     https://www.r-bloggers.com/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/
#  \n\n",
#     "   for details.\n\n"))
# }

#' mclapply
#' 
#' If the OS is Windows, set mclapply to the
#' the hackish version. Otherwise, leave the
#' definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {mclapply.hack}, 
                    Linux   = {mclapply},
                    Darwin  = {mclapply})
