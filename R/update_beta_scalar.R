update_beta_scalar <- function(response, preds, n, np, nj, alpha, betas, vari,
                               psis, phis, var_beta, var_alpha)
{
  tmp <- .C("update_beta_scalar", response=as.double(response), preds=as.double(preds),
            n=as.integer(n), np=as.integer(np), nj=as.integer(nj), alpha=as.double(alpha),
            betas=as.double(betas), vari=as.double(vari), psis=as.double(psis),
            phis=as.double(phis), var_beta=as.double(var_beta), var_alpha=as.double(var_alpha))
  return(list(alpha=tmp$alpha, betas=tmp$betas, vari=tmp$vari))
}