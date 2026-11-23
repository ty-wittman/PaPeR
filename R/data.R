#' Simulated example population
#'
#' A simulated population with overlapping generations, genotyped on both
#' microhaplotypes and SNPs, with the true pedigree for checking results.
#' 150 founders were born in 2015, then 110 offspring a year from 2016 to
#' 2020. Parents breed from age 1 to 4, and sires can have several mates a
#' year. 595 of the 700 individuals (85\%) are genotyped, so some true parents
#' cannot be found. Genotypes include dropout, miscalls and 2\% missing calls.
#' The script that builds it is \code{data-raw/make_example_data.R} in the
#' source repository.
#'
#' @format A list with four elements:
#' \describe{
#'   \item{\code{GenoM}}{character matrix, 595 individuals by 100
#'     microhaplotype loci with 3 to 8 alleles. Calls are \code{"a/b"};
#'     \code{"0"} is missing. About 0.6\% of true heterozygotes read as
#'     homozygotes, and 0.05\% of allele copies are miscalled.}
#'   \item{\code{GenoM_snp}}{character matrix, the same individuals by 150
#'     SNPs (alleles 1 and 2, minor allele frequency 0.1 to 0.5). About 0.5\%
#'     of genotypes are wrong.}
#'   \item{\code{LifeHistData}}{data frame with \code{ID}, \code{Sex} (1
#'     female, 2 male) and \code{BirthYear} for the genotyped individuals.}
#'   \item{\code{Pedigree}}{the true pedigree of all 700 individuals:
#'     \code{id}, \code{dam}, \code{sire}, \code{Sex}, \code{BirthYear} and
#'     \code{genotyped}.}
#' }
#' @examples
#' str(pap_example, max.level = 1)
#' pap_example$GenoM[1:3, 1:6]
#' table(pap_example$LifeHistData$BirthYear)
"pap_example"
