#' @describeIn doseSimulation Find the dose required to hit a target.
#'
#' @param fit the tdmorefit object
#' @param regimen the regimen
#' @param doseRows which rows of the regimen to adapt when searching for a new dose, or NULL to take the last one
#' @param interval which interval to search a dose in. Defaults to a ridiculously high range
#' @param target target value, as a data.frame
#' @param ... extra arguments passed to stats::uniroot
#'
#' @return a recommendation object
#' @export
#' @importFrom stats uniroot
findDoseCautiously <- function(fit, regimen=fit$regimen, doseRows=NULL, interval=c(0, 1E10), target, ...) {
  if(is.null(doseRows))
    doseRows <- nrow(regimen)

  if(nrow(target) > 1) {
    stop("Cannot find the dose to hit multiple targets. Split the treatment regimen and perform findDose separately per target.")
  }
  if(ncol(target[, colnames(target) != "TIME", drop=FALSE]) > 1) {
    stop("Cannot find the dose to hit multiple targets. Split the treatment regimen and perform findDose separately per target.")
  }

  # Find the best dose for the estimated parameters
  rootFunction <- function(AMT) {
    myRegimen <- updateRegimen(regimen = regimen, doseRows = doseRows, newDose = AMT)
    obs <- stats::predict(fit, newdata = target, regimen = myRegimen)
    result <- obs[, colnames(obs) != "TIME", drop=TRUE] - target[, colnames(target) != "TIME", drop=TRUE]
    if(length(result) > 1) stop("Cannot use findDose to hit multiple targets!")
    result
  }
  result <- runUniroot(rootFunction, interval, ...)
  return(convertResultToRecommendation(fit, result, regimen, doseRows, target))
}

runUniroot <- function(rootFunction, interval, ...) {
  iValues <- c(rootFunction(interval[1]), rootFunction(interval[2]))

  if(sign(iValues[1]) == sign(iValues[2])) {
    warning("Predicted values at edges of interval both ",
            switch(sign(iValues[1]), "0"="equal to", "1"="above", "2"="below"),
            " target, returning closest value...")
    i <- which.min(abs(iValues))

    result <- list(
      root=interval[i],
      f.root=iValues[i],
      iter=0,
      estim.prec=Inf
    )
    return(result)
  } else {
    return(uniroot(f=rootFunction, interval=interval, ...))
  }
}

#' Update a regimen with the specified dose.
#'
#' @param regimen the regimen to update
#' @param doseRows which rows of the regimen to adapt, if not specified, the last dose will be adapted
#' @param newDose the specified new dose
#'
#' @return the updated regimen
#' @noRd
updateRegimen <- function(regimen, doseRows = NULL, newDose) {
  if (is.null(doseRows))
    doseRows <- nrow(regimen)

  dose <- numeric(length = length(doseRows)) + as.numeric(newDose)
  names(dose) <- paste0(doseRows, ".AMT")
  updatedRegimen <- transformRegimen(regimen, dose)

  return(updatedRegimen)
}

#' @describeIn doseSimulation Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization of multiple doses.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#' @param targetMetadata defined target troughs as list(min=X, max=Y). If NULL or all NA, taken from the model metadata.
#'
#' @export
findDosesWithCaution <- function(fit, regimen=fit$regimen, targetMetadata=NULL) {
  if(! "FIX" %in% colnames(regimen) ) regimen$FIX <- FALSE
  if(is.null(targetMetadata) || all(is.na(targetMetadata))) {
    targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
    if(is.null(targetMetadata)) stop("No target defined in model metadata")
  }
  stopifnot( all( c("min", "max") %in% names(targetMetadata) ) )

  rowNumber <- regimen %>%
    dplyr::mutate(ROW_NO=dplyr::row_number()) %>%
    dplyr::filter(!FIX) %>%
    dplyr::filter(OCC==OCC[1]) %>%
    dplyr::pull(ROW_NO) %>% dplyr::last()

  target <- list(
    TIME=getTroughs(fit$tdmore, regimen[rowNumber, ], adj=TRUE)
  )

  outputVar <- fit$tdmore$res_var[[1]]$var
  targetValue <- mean( c(targetMetadata$min, targetMetadata$max) )
  if(is.na(targetValue)) stop("Target not defined, cannot optimize treatment...")
  target[outputVar] <- targetValue
  target <- tibble::as_tibble(target)
  print(target)

  modified <- rep(FALSE, nrow(regimen))
  result <- list()
  #step-wise: row per row
  for(i in seq_along(target$TIME)) {
    row <- target[i, ]
    iterationRows <- which( regimen$FIX == FALSE & regimen$TIME < row$TIME & !modified )
    if(length(iterationRows) == 0) next
    rec <- findDose(fit, regimen, iterationRows, target=row)
    regimen <- rec$regimen
    roundedAmt <- purrr::pmap_dbl(list(regimen$AMT, regimen$FORM), function(amt, form) {
      form <- tdmore::getMetadataByName(fit$tdmore, form)
      form$round_function(amt)
    })
    regimen$AMT[iterationRows] <- roundedAmt[iterationRows] #rounded amounts only
    modified[ iterationRows ] <- TRUE
    result <- c( result, rec$result )
  }

  return(structure(
    list(
      tdmorefit=fit,
      dose=regimen$AMT[ modified ],
      regimen = regimen,
      target=target,
      result=result
    ),
    class = c("recommendation")
  ))
}
