update_beta <- function(model, response, preds, n, np, nj, betas, vari, psis, phis, var_beta) {
  if (model == "rw2") {
    tmp <- .C("update_beta", response=as.double(response), preds=as.double(preds),
              n=as.integer(n), np=as.integer(np), nj=as.integer(nj), betas=as.double(betas),
              vari=as.double(vari), psis=as.double(psis), phis=as.double(phis),
              var_beta=as.double(var_beta))
  } else {
    tmp <- .C("update_beta_iid", response=as.double(response), preds=as.double(preds),
              n=as.integer(n), np=as.integer(np), nj=as.integer(nj), betas=as.double(betas),
              vari=as.double(vari), psis=as.double(psis), phis=as.double(phis),
              var_beta=as.double(var_beta))
  }
  return(list(betas=tmp$betas, vari=tmp$vari))
}