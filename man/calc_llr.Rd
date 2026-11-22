\name{calc_llr}
\alias{calc_llr}
\title{Score the assigned parents of a fixed pedigree}
\usage{
calc_llr(Pedigree, GenoM, LifeHistData, Dropout = NULL, Foreign = NULL,
         control = list())
}
\arguments{
\item{Pedigree}{Data frame with \code{id}, \code{dam}, and \code{sire}.
All IDs and assigned parents must be genotyped. Omitted genotyped rows are founders.}
\item{GenoM}{Character genotype matrix, one column per locus, row names as IDs,
calls \code{"a/b"}, \code{"0"} missing.}
\item{LifeHistData}{Data frame with \code{ID}, \code{Sex}, \code{BirthYear},
optionally \code{BY.min}, \code{BY.max}, \code{Year.last}; see \code{\link{assign_parents}}.}
\item{Dropout, Foreign}{Positive panel-wide error rates; see \code{\link{assign_parents}}.}
\item{control}{Settings from \code{\link{pap_control}}. \code{evidence} selects
\code{"ancestral"} (default) or \code{"full"}; \code{hypotheses} selects
\code{"S7"}, \code{"S37"}, or \code{"S45"}. Error calibration and search
settings do not initiate a search or fit new rates.}
}
\description{
Computes read-only pedigree likelihood relationship LLRs with the supplied error rates.
Each assigned parent is compared with the suite's nonparent alternatives,
keeping the other assigned parent. Full evidence includes descendants through
forest messages; forest cuts approximate external correlations in loopy pedigrees.
}
\value{
The supplied pedigree with base-10 \code{LLRdam}, \code{LLRsire},
\code{LLRpair} (the smaller side, only with both parents), \code{OHdam},
\code{OHsire}, \code{MEpair}, \code{DamLLRStatus}, \code{SireLLRStatus}, and
\code{PairLLRStatus}. Status is \code{"available"}, \code{"not_assigned"},
or \code{"unavailable"}. Unassigned and unavailable scores are \code{NA}.
Rows and parent links are preserved. Attributes \code{final_llr_evidence},
\code{final_llr_hypotheses}, and \code{final_llr_unavailable} record the settings
and native availability counters.
}
