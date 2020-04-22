# Package environment to store a couple of variables that must be global
metaEnv <- new.env(parent=emptyenv())
assign("VERBOSE",NULL,envir=metaEnv)
assign("LOGGER",NULL,envir=metaEnv)
