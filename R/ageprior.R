# Sequoia-style first-pass and scaffold-rebuilt age priors.
# PaPeR keeps likelihoods in natural logs; ratios remain on the ordinary
# probability scale until the scoring boundary.

.parent_bounds <- function(lh, MinAgeParent = c(1L, 1L),
                           MaxAgeParent = NULL) {
  minp <- .normalize_age_vector(MinAgeParent, "MinAgeParent")
  if (any(minp < 1L))
    stop("MinAgeParent must contain one or two positive whole numbers", call. = FALSE)
  yrs <- c(lh$BirthYear, lh$BY.min, lh$BY.max)
  yrs <- yrs[!is.na(yrs)]
  span <- if (length(yrs)) max(1L, diff(range(yrs)) + 1L) else 1L
  if (is.null(MaxAgeParent)) maxp <- rep(span, 2L)
  else maxp <- .normalize_age_vector(MaxAgeParent, "MaxAgeParent")
  if (any(maxp < minp))
    stop("MaxAgeParent must contain one or two whole numbers not below MinAgeParent",
         call. = FALSE)
  names(minp) <- names(maxp) <- c("M", "P")
  list(min = minp, max = maxp, inferred = is.null(MaxAgeParent))
}

