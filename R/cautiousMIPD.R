
isAdvagraf <- function(form) {
  values <- c("1", "2", "3", "4")
  if (!all(form %in% values)) stop("formulation must be 1, 2, 3 or 4")
  return(form %in% c("3"))
}

getMaxAmountFromPast <-function(regimen, fit, min_dose) {
  fixedRegimen <- regimen %>%
    dplyr::filter(FIX)

  if(nrow(fixedRegimen)==0) {
    return(min_dose)
  }

  normalisedAmounts <- purrr::pmap_dbl(list(fixedRegimen$AMT, fixedRegimen$FORM), function(amt, form) {
    if (isAdvagraf(form)) {
      return(amt/2)
    } else {
      return(amt)
    }
  })
  return(max(normalisedAmounts))
}

#' @describeIn doseSimulation Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization of multiple doses.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#' @param targetMetadata defined target troughs as list(min=X, max=Y). If NULL or all NA, taken from the model metadata.
#' @param span user-given time span (in hours) in which cautious dosing is applied
#' @param cap dose cap that cannot be exceeded (unit: prograf dose)
#' @param ff fold factor threshold (recommended doses compared with max dose already given)
#' @param ff_min_dose min dose to use together with the fold factor check
#' @export
findDosesCautiously <- function(fit, regimen=fit$regimen, targetMetadata=NULL, span, cap, ff, ff_min_dose=3) {
  if(! "FIX" %in% colnames(regimen) ) regimen$FIX <- FALSE
  if(is.null(targetMetadata) || all(is.na(targetMetadata))) {
    targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
    if(is.null(targetMetadata)) stop("No target defined in model metadata")
  }
  stopifnot( all( c("min", "max") %in% names(targetMetadata) ) )

  # Retrieving target value
  targetValue <- mean( c(targetMetadata$min, targetMetadata$max) )
  if(is.na(targetValue)) stop("Target not defined, cannot optimize treatment...")

  obsConc <- fit$observed$CONC

  # Retrieving the last observed concentration
  lastConc <- obsConc[obsConc %>% length()]

  # Getting list of troughs
  unfixRegimen <- regimen[regimen$FIX==FALSE, ]
  target <- list(
    TIME=getTroughs(fit$tdmore, unfixRegimen, adj=TRUE)
  )

  # Apply cautious dosing if last concentration is lower than the target value
  if (lastConc < targetValue && nrow(unfixRegimen) > 0 && length(target$TIME) > 0) {
    targetedTime <- unfixRegimen$TIME[[1]] + span
    # Overwrite target by finding the closest trough value to the targeted time
    target <- list(
      TIME=target$TIME[which.min(abs(target$TIME - targetedTime))]
    )
  }

  outputVar <- fit$tdmore$res_var[[1]]$var
  target[outputVar] <- targetValue
  target <- tibble::as_tibble(target)
  print(target)

  modified <- rep(FALSE, nrow(regimen))
  result <- list()
  #step-wise: row per row
  for(i in seq_along(target$TIME)) {
    row <- target[i, ]
    iterationRows <- which( regimen$FIX == FALSE & regimen$TIME < row$TIME & !modified )
    maxAllowedAmount <- getMaxAmountFromPast(regimen=regimen, fit=fit, min_dose=ff_min_dose)*ff
    if(maxAllowedAmount > cap) {
      maxAllowedAmount <- cap # Cannot exceed dose cap
    }
    if(length(iterationRows) == 0) next
    # browser()
    # Add WEIGHT column for Advagraf
    regimen <- regimen %>%
      dplyr::mutate(WEIGHT=dplyr::if_else(isAdvagraf(FORM), 2, 1))

    rec <- findDose(fit=fit, regimen=regimen, doseRows=iterationRows, weightForm=TRUE, target=row)
    regimen <- rec$regimen

    roundedAmt <- purrr::pmap_dbl(list(regimen$AMT, regimen$FORM), function(amt, form) {
      formulation <- tdmore::getMetadataByName(fit$tdmore, form)
      if (isAdvagraf(form)) {
        doseCap_ <- formulation$round_function(min(cap*2, maxAllowedAmount*2))
      } else {
        doseCap_ <- formulation$round_function(min(cap, maxAllowedAmount))
      }
      retValue <- formulation$round_function(amt)
      if (retValue > doseCap_) {
        retValue <- doseCap_
      }
      return(retValue)
    })
    regimen$AMT[iterationRows] <- roundedAmt[iterationRows] # Rounded amounts only
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
