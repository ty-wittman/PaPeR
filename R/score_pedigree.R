# Truth scoring for typed parentage pedigrees.

#' Score an inferred parentage pedigree against the truth
#'
#' @param inferred data.frame(id, dam, sire) from \code{assign_parents}.
#' @param truth data.frame(id, dam, sire) ground truth.
#' @param typed_ids optional IDs of all genotyped individuals. When supplied,
#'   true parents that are not genotyped count as no parent, so assigning a
#'   parent on such a side is a false assignment.
#' @return list with per-parent (dam, sire) counts: assigned, correct, wrong,
#'   missed, precision, recall, and true-null assignment accuracy.
#' @examples
#' ex <- pap_example
#' first <- ex$LifeHistData$ID[ex$LifeHistData$BirthYear <= 2016]
#' ped <- assign_parents(ex$GenoM[first, ], ex$LifeHistData,
#'                       Dropout = 0.005, Foreign = 5e-4)
#' str(score_pedigree(ped, ex$Pedigree, typed_ids = first))
#' @export
score_pedigree <- function(inferred, truth, typed_ids = NULL) {
  if (!is.null(typed_ids)) truth <- .typed_truth(truth, typed_ids)
  tr <- truth[match(inferred$id, truth$id), , drop = FALSE]
  one <- function(inf, tru) {
    truth <- tru
    genotyped_truth <- !is.na(truth)
    assigned <- !is.na(inf)
    correct  <- assigned & genotyped_truth & inf == truth
    correct_null <- !assigned & is.na(truth)
    wrong    <- assigned & (is.na(truth) | inf != truth)
    missed   <- genotyped_truth & !correct
    true_null <- is.na(truth)
    list(assigned = sum(assigned), correct = sum(correct),
         wrong = sum(wrong), missed = sum(missed),
         precision = sum(correct) / max(1, sum(assigned)),
         recall = sum(correct) / max(1, sum(genotyped_truth)),
         true_null = sum(true_null), correct_null = sum(correct_null),
         false_positive_true_null = sum(assigned & true_null),
         null_accuracy = sum(correct_null) / max(1, sum(true_null)))
  }
  list(dam = one(inferred$dam, tr$dam), sire = one(inferred$sire, tr$sire))
}

# The truth restricted to genotyped individuals: a true parent that is not
# genotyped becomes NA.
.typed_truth <- function(truth, typed_ids) {
  if (!is.data.frame(truth) || !all(c("id", "dam", "sire") %in% names(truth)))
    stop("truth must contain id, dam, and sire", call. = FALSE)
  typed_ids <- unique(as.character(typed_ids))
  if (anyNA(typed_ids) || any(!nzchar(typed_ids)))
    stop("typed_ids must be nonmissing nonempty IDs", call. = FALSE)
  out <- truth[, c("id", "dam", "sire"), drop = FALSE]
  out$id <- as.character(out$id)
  out <- out[out$id %in% typed_ids, , drop = FALSE]
  for (side in c("dam", "sire")) {
    value <- as.character(out[[side]])
    value[is.na(value) | !nzchar(value) | value == "0" | !value %in% typed_ids] <- NA_character_
    out[[side]] <- value
  }
  out
}
