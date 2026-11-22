# Central life-history validation and discrete birth-year representation.

.whole_or_na <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  is.na(y) | (is.finite(y) & y == floor(y))
}

#' Validate and normalize life-history data
#'
#' Sex is normalized to 1=female, 2=male, and NA=unknown; both input 3 and NA
#' map to the same unknown state. Birth-year bounds are inclusive integer bounds.
#' @keywords internal
#' @noRd
.normalize_life_history <- function(LifeHistData, ids = NULL) {
  if (!is.data.frame(LifeHistData))
    stop("LifeHistData must be a data.frame", call. = FALSE)
  need <- c("ID", "Sex", "BirthYear")
  if (!all(need %in% names(LifeHistData)))
    stop("LifeHistData must contain ID, Sex, and BirthYear", call. = FALSE)
  lh <- LifeHistData
  lh$ID <- as.character(lh$ID)
  if (anyNA(lh$ID) || any(!nzchar(lh$ID)) || anyDuplicated(lh$ID))
    stop("LifeHistData IDs must be unique, non-missing strings", call. = FALSE)
  if (!"BY.min" %in% names(lh)) lh$BY.min <- NA_integer_
  if (!"BY.max" %in% names(lh)) lh$BY.max <- NA_integer_
  if (!"Year.last" %in% names(lh)) lh$Year.last <- NA_integer_

  sx <- suppressWarnings(as.integer(as.character(lh$Sex)))
  bad.sex <- !is.na(lh$Sex) & (is.na(sx) | !sx %in% 1:3)
  if (any(bad.sex))
    stop("Sex must be 1=female, 2=male, or 3/NA=unknown", call. = FALSE)
  sx[sx == 3L] <- NA_integer_
  lh$Sex <- sx

  for (z in c("BirthYear", "BY.min", "BY.max", "Year.last")) {
    if (any(!.whole_or_na(lh[[z]])))
      stop(z, " must contain whole numbers or NA", call. = FALSE)
    lh[[z]] <- suppressWarnings(as.integer(as.character(lh[[z]])))
  }
  one.bound <- xor(is.na(lh$BY.min), is.na(lh$BY.max))
  if (any(one.bound))
    stop("BY.min and BY.max must either both be supplied or both be NA", call. = FALSE)
  if (any(lh$BY.min > lh$BY.max, na.rm = TRUE))
    stop("BY.min must be less than or equal to BY.max", call. = FALSE)
  bounded <- !is.na(lh$BirthYear) & !is.na(lh$BY.min)
  if (any(bounded & (lh$BirthYear < lh$BY.min | lh$BirthYear > lh$BY.max)))
    stop("BirthYear must lie inside the inclusive BY.min/BY.max interval", call. = FALSE)
  exact <- !is.na(lh$BirthYear)
  lh$BY.min[exact] <- lh$BirthYear[exact]
  lh$BY.max[exact] <- lh$BirthYear[exact]

  if (!is.null(ids)) {
    ids <- as.character(ids)
    extra <- setdiff(ids, lh$ID)
    if (length(extra))
      lh <- rbind(lh, data.frame(ID = extra, Sex = NA_integer_,
        BirthYear = NA_integer_, BY.min = NA_integer_, BY.max = NA_integer_,
        Year.last = NA_integer_, stringsAsFactors = FALSE))
    lh <- lh[match(ids, lh$ID), , drop = FALSE]
  }
  rownames(lh) <- NULL
  lh
}

.by_mass <- function(birth, lo, hi) {
  if (!is.na(birth)) return(setNames(1, as.character(birth)))
  if (is.na(lo) || is.na(hi)) return(NULL)
  yy <- seq.int(lo, hi)
  setNames(rep(1 / length(yy), length(yy)), as.character(yy))
}

.require_birth_bounds <- function(lh) {
  missing <- is.na(lh$BirthYear) & (is.na(lh$BY.min) | is.na(lh$BY.max))
  if (any(missing)) stop(
    "Each typed individual requires BirthYear or both BY.min and BY.max; missing for: ",
    paste(head(lh$ID[missing], 10L), collapse = ", "), call. = FALSE)
  invisible(lh)
}

.life_history_masses <- function(lh) {
  out <- lapply(seq_len(nrow(lh)), function(i)
    .by_mass(lh$BirthYear[i], lh$BY.min[i], lh$BY.max[i]))
  names(out) <- lh$ID
  out
}

.normalize_age_vector <- function(x, name) {
  if (is.null(x) || !length(x) || length(x) > 2L)
    stop(name, " must contain one or two whole numbers", call. = FALSE)
  y <- suppressWarnings(as.numeric(as.character(x)))
  if (anyNA(y) || any(!is.finite(y)) || any(y != floor(y)))
    stop(name, " must contain whole numbers", call. = FALSE)
  as.integer(rep(y, length.out = 2L))
}


.birth_year_support <- function(lh, id) {
  i <- match(as.character(id), lh$ID)
  if (is.na(i)) stop("birth-year support requested for an unknown ID", call. = FALSE)
  if (!is.na(lh$BirthYear[i])) return(as.integer(lh$BirthYear[i]))
  if (is.na(lh$BY.min[i]) || is.na(lh$BY.max[i])) return(NULL)
  seq.int(lh$BY.min[i], lh$BY.max[i])
}

.parent_age_feasible <- function(lh, offspring, parent, role,
                                 min_age_parent, max_age_parent) {
  role <- match.arg(role, c("M", "P"))
  offspring <- as.character(offspring)[1L]
  parent <- as.character(parent)[1L]
  if (is.na(offspring) || is.na(parent) || offspring == parent) return(FALSE)
  io <- match(offspring, lh$ID); ip <- match(parent, lh$ID)
  if (is.na(io) || is.na(ip)) return(FALSE)

  minp <- setNames(.normalize_age_vector(min_age_parent, "MinAgeParent"), c("M", "P"))
  maxp <- setNames(.normalize_age_vector(max_age_parent, "MaxAgeParent"), c("M", "P"))
  if (any(maxp < minp))
    stop("MaxAgeParent must not be below MinAgeParent", call. = FALSE)

  child_years <- .birth_year_support(lh, offspring)
  parent_years <- .birth_year_support(lh, parent)
  if (is.null(child_years)) return(TRUE)
  if (is.null(parent_years)) return(TRUE)

  ok <- outer(parent_years, child_years, function(year_parent, year_child) {
    age <- year_child - year_parent
    age >= minp[role] & age <= maxp[role]
  })
  any(ok)
}

.parent_year_last_feasible <- function(lh, offspring, parent) {
  offspring <- as.character(offspring)[1L]
  parent <- as.character(parent)[1L]
  if (is.na(offspring) || is.na(parent) || offspring == parent) return(FALSE)
  io <- match(offspring, lh$ID); ip <- match(parent, lh$ID)
  if (is.na(io) || is.na(ip)) return(FALSE)
  last <- lh$Year.last[ip]
  if (is.na(last)) return(TRUE)
  child_years <- .birth_year_support(lh, offspring)
  is.null(child_years) || any(child_years <= last)
}

.parent_age_year_feasible <- function(lh, offspring, parent, role,
                                      min_age_parent, max_age_parent) {
  role <- match.arg(role, c("M", "P"))
  offspring <- as.character(offspring)[1L]
  parent <- as.character(parent)[1L]
  if (is.na(offspring) || is.na(parent) || offspring == parent) return(FALSE)
  io <- match(offspring, lh$ID); ip <- match(parent, lh$ID)
  if (is.na(io) || is.na(ip)) return(FALSE)

  minp <- setNames(.normalize_age_vector(min_age_parent, "MinAgeParent"), c("M", "P"))
  maxp <- setNames(.normalize_age_vector(max_age_parent, "MaxAgeParent"), c("M", "P"))
  if (any(maxp < minp))
    stop("MaxAgeParent must not be below MinAgeParent", call. = FALSE)
  child_years <- .birth_year_support(lh, offspring)
  parent_years <- .birth_year_support(lh, parent)
  last <- lh$Year.last[ip]
  if (is.null(child_years)) return(TRUE)
  if (is.null(parent_years))
    return(is.na(last) || any(child_years <= last))
  ok <- outer(parent_years, child_years, function(year_parent, year_child) {
    age <- year_child - year_parent
    possible <- age >= minp[role] & age <= maxp[role]
    if (!is.na(last)) possible <- possible & year_child <= last
    possible
  })
  any(ok)
}

.log_probability_ln <- function(probability, context = "probability") {
  if (!is.numeric(probability) || is.object(probability) || anyNA(probability) ||
      any(!is.finite(probability)) || any(probability < 0))
    stop(context, " must be nonnegative finite", call. = FALSE)
  out <- rep.int(-Inf, length(probability))
  positive <- probability > 0
  out[positive] <- log(probability[positive])
  out
}

.log_likelihood_contrast_ln <- function(focal_ln, alternative_ln,
                                         context = "log likelihood contrast") {
  invalid <- function(x) {
    !is.numeric(x) || is.object(x) || anyNA(x) || any(is.nan(x)) ||
      any(is.infinite(x) & x > 0)
  }
  if (invalid(focal_ln) || invalid(alternative_ln))
    stop(context, " inputs must be finite or -Inf", call. = FALSE)
  if (length(focal_ln) != length(alternative_ln) &&
      length(focal_ln) != 1L && length(alternative_ln) != 1L)
    stop(context, " inputs must have equal lengths or length one", call. = FALSE)
  out <- focal_ln - alternative_ln
  # A scored zero against positive support is a defined -Inf contrast.
  # Two scored zeros supply no ratio: preserve unknown evidence as NA rather
  # than turning it into a rejection, and do not expose arithmetic NaN.
  both_zero <- (is.infinite(focal_ln) & focal_ln < 0) &
    (is.infinite(alternative_ln) & alternative_ln < 0)
  out[both_zero] <- NA_real_
  out
}

.logsumexp <- function(x) {
  if (!is.numeric(x) || is.object(x) || anyNA(x) || any(is.nan(x)) ||
      any(is.infinite(x) & x > 0))
    stop("log-sum-exp inputs must be finite or -Inf", call. = FALSE)
  if (!length(x) || all(is.infinite(x) & x < 0)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}
