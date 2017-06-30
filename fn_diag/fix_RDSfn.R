# Functions for reading and writing RDS (specifically within fn_diag)
# -----------------------------------------------------------------------------
savebroadRDS <- function(tname) {
  setwd(git.dir)
  setwd(reponame)
  saveRDS(fin.res[[1]]$result,
    paste0("fn_diag/", fin.res[[1]]$data, "-broad-", tname, ".rds")
  )
  saveRDS(fin.res[[2]]$result,
    paste0("fn_diag/", fin.res[[2]]$data, "-broad-", tname, ".rds")
  )
  saveRDS(fin.res[[3]]$result,
    paste0("fn_diag/", fin.res[[3]]$data, "-broad-", tname, ".rds")
  )
  saveRDS(fin.res[[4]]$result,
    paste0("fn_diag/", fin.res[[4]]$data, "-broad-", tname, ".rds")
  )
  saveRDS(fin.res[[5]]$result,
    paste0("fn_diag/", fin.res[[5]]$data, "-broad-", tname, ".rds")
  )
}

readbroadRDS <- function(tname) {
  setwd(git.dir)
  setwd(reponame)
  d1a <<- saveRDS(paste0("fn_diag/d1a-broad-", tname, ".rds"))
  d2a <<- saveRDS(paste0("fn_diag/d2a-broad-", tname, ".rds"))
  d3a <<- saveRDS(paste0("fn_diag/d3a-broad-", tname, ".rds"))
  d2b <<- saveRDS(paste0("fn_diag/d2b-broad-", tname, ".rds"))
  d3b <<- saveRDS(paste0("fn_diag/d3b-broad-", tname, ".rds"))
}
