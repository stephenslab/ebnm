# mle_point_normal <- function(x, s, g, control,
#                              fix_pi0, fix_a, fix_mu,
#                              optmethod, use_grad, use_hess) {
#   fix_par <- c(fix_pi0, fix_a, fix_mu)
#
#   retlist <- mle_parametric(x, s, g, fix_par,
#                             pn_startpar, pn_precomp, pn_nllik, pn_gfromopt,
#                             optmethod, control, use_grad, use_hess)
#
#   # Check the solution pi0 = 1.
#   if (!fix_pi0 && fix_mu) {
#     pi0_val <- sum(-0.5 * log(2 * pi * s^2) - 0.5 * (x - g$mu)^2 / s^2)
#     if (pi0_val > retlist$val) {
#       retlist$pi0 <- 1
#       retlist$a <- 1
#       retlist$val <- pi0_val
#     }
#   }
#
#   return(retlist)
# }
