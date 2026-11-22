#' PaPeR: parentage assignment for microhaplotypes and SNPs
#'
#' PaPeR assigns genotyped dams and sires to genotyped individuals from
#' multiallelic microhaplotype or biallelic SNP genotypes and life-history
#' data (sex and birth year).
#'
#' PaPeR is inspired by the parentage step of Sequoia (Huisman 2017) and
#' follows its general heuristic approach: rounds of greedy assignment, each
#' parent tested against alternative relationships, age limits from birth
#' years, and Sequoia-style inputs, settings and output columns, so that it is
#' easy to pick up for people who already use Sequoia. Its likelihood
#' machinery differs substantially. Each decision is an exact local
#' pedigree-likelihood calculation that conditions on the parents already
#' assigned, multiallelic markers are modelled directly with a dropout/foreign
#' error model, and the error rate is calibrated from the inferred pedigree.
#' PaPeR is an independent implementation and does not include Sequoia code.
#' It finds genotyped parents only: unlike Sequoia, it does not cluster
#' sibships or add dummy parents.
#'
#' Start with \code{\link{assign_parents}} and the simulated
#' \code{\link{pap_example}} data. Other functions: \code{\link{pap_control}}
#' (settings), \code{\link{MHStats}} (allele frequencies and marker QC),
#' \code{\link{find_duplicates}}, \code{\link{get_maybe_rel}},
#' \code{\link{calc_llr}} (score a fixed pedigree), \code{\link{calc_pair_LL}},
#' \code{\link{calc_trio_LL}} and \code{\link{score_pedigree}}.
#'
#' @references
#' Huisman, J. (2017). Pedigree reconstruction from SNP data: parentage
#' assignment, sibship clustering and beyond. \emph{Molecular Ecology
#' Resources}, 17(5), 1009--1024. \doi{10.1111/1755-0998.12665}
#'
#' @seealso The Sequoia R package:
#'   \url{https://CRAN.R-project.org/package=sequoia} and
#'   \url{https://jiscah.github.io/}.
#' @keywords internal
"_PACKAGE"
