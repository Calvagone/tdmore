
isAdvagraf <- function(form) {
  values <- c("1", "2", "3", "4")
  if (!all(form %in% values)) stop("formulation must be 1, 2, 3 or 4")
  return(form %in% c("3"))
}

getMaxAmountFromPast <-function(regimen, fit) {
  fixedRegimen <- regimen %>%
    dplyr::filter(FIX)

  if(nrow(fixedRegimen)==0) {
    return(3000)
  }

  normalisedAmounts <- purrr::pmap_dbl(list(fixedRegimen$AMT, fixedRegimen$FORM), function(amt, form) {
    formulation <- tdmore::getMetadataByName(fit$tdmore, form)
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
#' @param doseCap dose cap that cannot be exceeded (unit: prograf dose)
#' @param foldFactor fold factor threshold (recommended doses compared with max dose already given)
#' @export
findDosesCautiously <- function(fit, regimen=fit$regimen, targetMetadata=NULL, doseCap, foldFactor) {
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

  if (lastConc < targetValue) {
    # smoother <- 0L # target=morning dose trough of current occasion
    smoother <- 1L   # target=morning dose trough of next occasion
  } else {
    smoother <- NULL # Free-style
  }

  # Unfixed part of regimen
  unfixedRegimen <- regimen %>%
    dplyr::mutate(ROW_NO=dplyr::row_number()) %>%
    dplyr::filter(!FIX)

  if (nrow(unfixedRegimen) > 0 && !is.null(smoother)) {
    currentOcc <- unfixedRegimen %>% dplyr::pull(OCC) %>% dplyr::first()
    occasions <- unfixedRegimen$OCC
    # Find latest occasion available from data
    for (occ in (((0:smoother) + currentOcc) %>% rev())) {
      if (occ %in% occasions) break
    }
    rowNumber <- unfixedRegimen %>%
      dplyr::filter(OCC==occ) %>%
      dplyr::pull(ROW_NO) %>% dplyr::last()
    target <- list(
      TIME=getTroughs(fit$tdmore, regimen[rowNumber, ], adj=TRUE)
    )
  } else {
    target <- list(
      TIME=getTroughs(fit$tdmore, regimen[regimen$FIX==FALSE, ], adj=TRUE)
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
    maxAllowedAmount <- getMaxAmountFromPast(regimen=regimen, fit=fit)*foldFactor
    if(maxAllowedAmount > doseCap) {
      maxAllowedAmount <- doseCap # Cannot exceed dose cap
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
        doseCap_ <- formulation$round_function(min(doseCap*2, maxAllowedAmount*2))
      } else {
        doseCap_ <- formulation$round_function(min(doseCap, maxAllowedAmount))
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
