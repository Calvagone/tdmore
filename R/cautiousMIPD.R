
#' @describeIn doseSimulation Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization of multiple doses.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#' @param targetMetadata defined target troughs as list(min=X, max=Y). If NULL or all NA, taken from the model metadata.
#'
#' @export
findDosesCautiously <- function(fit, regimen=fit$regimen, targetMetadata=NULL) {
  if(! "FIX" %in% colnames(regimen) ) regimen$FIX <- FALSE
  if(is.null(targetMetadata) || all(is.na(targetMetadata))) {
    targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
    if(is.null(targetMetadata)) stop("No target defined in model metadata")
  }
  stopifnot( all( c("min", "max") %in% names(targetMetadata) ) )
  # stopifnot( is.integer(smoother) && smoother >= 0L )

  obsConc <- fit$observed$CONC

  # Retrieving the last observed concentration
  lastConc <- obsConc[obsConc %>% length()]

  if (lastConc < 12.5) {
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
    iRoundedAmt <- roundedAmt[iterationRows] # Rounded amounts only

    # Capped dose when smoother is used (max 20mg/occasion)
    # Number of occasions = smoother + 1
    if (!is.null(smoother)) {
      maxAmount <- 20000*(smoother + 1)
      if (sum(iRoundedAmt) > maxAmount) {
        iRoundedAmt <- rep(maxAmount/(iRoundedAmt %>% length()), iRoundedAmt %>% length())
      }
    }
    regimen$AMT[iterationRows] <- iRoundedAmt
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
