#' Eight Schools Data
#'
#' The eight schools dataset analyzed by Rubin (1981) and Gelman et al. (2014).
#' In each of eight US high schools, different randomized experiments were
#' conducted in order to estimate the effect of coaching programs on
#' standardized test scores. The data includes estimates of treatment effects
#' as well as standard errors.
#'
#' @format
#' A data frame with 8 rows and 2 columns:
#' \describe{
#'   \item{x}{The treatment effect estimate.}
#'   \item{s}{The standard error.}
#' }
#'
#' @source Rubin DB (1981). "Estimation in Parallel Randomized Experiments."
#'   \emph{Journal of Educational Statistics}, \strong{6}(4), 377-401.
#'
#' @seealso Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB (2014).
#' \emph{Bayesian data analysis}. Texts in Statistical Science Series, third
#' edition. CRC Press, Boca Raton, FL.
#'
"eight_schools"

#' @name wOBA
#' 
#' @title 2022 MLB wOBA Data
#'
#' @docType data
#' 
#' @description Weighted on-base average (wOBA) for hitters from the
#'   2022 Major League Baseball (MLB) season. Standard errors are
#'   calculated by modeling hitting outcomes as multinomially
#'   distributed and plugging in empirical proportions as the \dQuote{true}
#'   outcome probabilities. To handle small sample sizes, standard
#'   errors are lower bounded by the errors that would be obtained using
#'   league-wide proportions rather than the plug-in estimates.
#'
#' @format
#' \code{wOBA} is a data frame with 688 rows and 6 columns:
#' \describe{
#' 
#'   \item{FanGraphsID}{The hitter's FanGraphs identifier.}
#' 
#'   \item{Name}{The hitter's name.}
#'   \item{Team}{The hitter's MLB team (given as a three-letter code) or
#'   \code{NA} if the hitter played for multiple teams.}
#'   \item{PA}{The hitter's number of plate appearances over the 2022 season.}
#'   \item{x}{The hitter's wOBA over the 2022 season.}
#' 
#'   \item{s}{The standard error for the hitter's wOBA.}
#' }
#'
#' @source <https://fangraphs.com>
#'
#' @keywords data
#'
#' @examples
#' data(wOBA)
#' summary(wOBA)
#' 
NULL
