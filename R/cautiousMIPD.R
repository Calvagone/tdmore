
getFormulationWeight <- function(form) {
  values <- c("1", "2", "3", "4")
  if (!all(form %in% values)) stop("formulation must be 1, 2, 3 or 4")
  weight <- rep(0, length(form))
  weight[form=="1"] <- 1
  weight[form=="2"] <- 1
  weight[form=="3"] <- 2
  weight[form=="4"] <- 10 # 20% bioavailability on average
  return(weight)
}

getMaxAmountFromPast <- function(regimen, fit, min_dose, loading_dose) {
  fixedRegimen <- regimen %>%
    dplyr::filter(FIX)

  if(nrow(fixedRegimen)==0) {
    return(min_dose)
  }

  # Process loading dose
  if (loading_dose) {
    fixedRegimen <- fixedRegimen %>%
      dplyr::mutate(AMT=dplyr::if_else(dplyr::row_number()==1, AMT/2, AMT))
  }

  # Retrieve normalised amounts
  normalisedAmounts <- getNormalisedAmounts(fixedRegimen)

  return(max(normalisedAmounts))
}

getNormalisedAmounts <- function(regimen) {
  retValue <- purrr::pmap_dbl(list(regimen$AMT, regimen$FORM), function(amt, form) {
    return(amt/getFormulationWeight(form))
  })
  return(retValue)
}

getLastAmountInFuture <- function(regimen) {
  unfixedRegimen <- regimen %>%
    dplyr::filter(!FIX)

  if (nrow(unfixedRegimen)==0) {
    return(NULL)
  }

  # Retrieve normalized amounts, reverse vector
  normalisedAmounts <- getNormalisedAmounts(unfixedRegimen) %>% rev()

  # Return last positive amount
  index <- which(normalisedAmounts > 0)
  if (length(index) > 0) {
    return(normalisedAmounts[index[1]])
  } else {
    return(NULL)
  }
}

getTrailingZeroIndexes <- function(x) {
  xRev <- rev(x)
  xRevIndexes <- rev(seq_along(xRev))

  index <- which(xRev != 0)
  if (length(index) > 0) {
    firstIndex <- index[1]
    if (firstIndex > 1) {
      indexes <- seq_len(firstIndex - 1)
    } else {
      indexes <- integer(0)
    }
  } else {
    indexes <- xRev==0
  }
  return(rev(xRevIndexes[indexes]))
}

# getTrailingZeroIndexes(c(1,2,3,0,0,0,0))
# getTrailingZeroIndexes(c(1,2,3,0,0,0,1))
# getTrailingZeroIndexes(c(0,0,0,0,0,0,0))
# getTrailingZeroIndexes(c(0,0,0,0,0,-9,0))

fillTrailingZeroesInFuture <- function(regimen, replacement) {
  if (nrow(regimen)==0 || is.null(replacement)) {
    return(regimen)
  }

  indexes <- getTrailingZeroIndexes(regimen$AMT)

  regimen <- regimen %>%
    dplyr::mutate(TRAILING_ZERO=(dplyr::row_number() %in% indexes) & !FIX) %>%
    dplyr::mutate(AMT=dplyr::if_else(TRAILING_ZERO, replacement*getFormulationWeight(FORM), AMT)) %>%
    dplyr::select(-TRAILING_ZERO)

  return(regimen)
}

#' @describeIn doseSimulation Optimize the regimen, using the metadata defined by the model
#'
#' This performs a step-wise optimization of multiple doses.
#'
#' @param fit tdmorefit object
#' @param regimen the treatment regimen to optimize
#' @param targetMetadata defined target troughs as list(min=X, max=Y). If NULL or all NA, taken from the model metadata.
#' @param cautious cautious mode, default is TRUE (=enabled). FALSE will disable the cautious mode.
#' @param ttt time-to-target (in hours) when cautious mode is used. Corresponds to the  time interval between the first dose to adjust and a trough concentration in the future.
#' @param cap dose cap that cannot be exceeded (unit: prograf dose)
#' @param ff fold factor threshold (recommended doses compared with max dose already given)
#' @param ff_min_dose min dose to use together with the fold factor check
#' @param ff_loading_dose loading dose is used, logical value. Default is FALSE.
#' @export
findDosesCautiously <- function(fit, regimen=fit$regimen, targetMetadata=NULL, cautious=TRUE, ttt, cap, ff, ff_min_dose=3, ff_loading_dose=FALSE) {
  if(! "FIX" %in% colnames(regimen) ) regimen$FIX <- FALSE
  if(is.null(targetMetadata) || all(is.na(targetMetadata))) {
    targetMetadata <- tdmore::getMetadataByClass(fit$tdmore, "tdmore_target")
    if(is.null(targetMetadata)) stop("No target defined in model metadata")
  }
  stopifnot( all( c("min", "max") %in% names(targetMetadata) ) )

  # All doses in future init to 0
  regimen <- regimen %>%
    dplyr::mutate(AMT=dplyr::if_else(!FIX, 0, AMT))

  # Retrieving target value
  targetValue <- mean( c(targetMetadata$min, targetMetadata$max) )
  if(is.na(targetValue)) stop("Target not defined, cannot optimize treatment...")

  # Getting list of troughs
  unfixRegimen <- regimen[regimen$FIX==FALSE, ]
  target <- list(
    TIME=getTroughs(fit$tdmore, unfixRegimen, adj=TRUE)
  )

  # Apply cautious dosing if last concentration is lower than the target value
  if (cautious && nrow(unfixRegimen) > 0 && length(target$TIME) > 0) {
    targetedTime <- unfixRegimen$TIME[[1]] + ttt
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
    maxAllowedAmount <- getMaxAmountFromPast(regimen=regimen, fit=fit, min_dose=ff_min_dose, loading_dose=ff_loading_dose)*ff
    if(maxAllowedAmount > cap) {
      maxAllowedAmount <- cap # Cannot exceed dose cap
    }
    if(length(iterationRows) == 0) next
    # browser()
    # Add WEIGHT column for Advagraf
    regimen <- regimen %>%
      dplyr::mutate(WEIGHT=getFormulationWeight(FORM))

    tacIVRateFun <-function(amt, form) {
      return(ifelse(form=="4", amt/24*1000, 0))
    }
    rec <- findDose(fit=fit, regimen=regimen, doseRows=iterationRows, weightForm=TRUE, rateFun=tacIVRateFun, target=row)

    regimen <- rec$regimen

    roundedAmt <- purrr::pmap_dbl(list(regimen$AMT, regimen$FORM), function(amt, form) {
      formulation <- tdmore::getMetadataByName(fit$tdmore, form)
      doseCap_ <- formulation$round_function(min(cap*getFormulationWeight(form), maxAllowedAmount*getFormulationWeight(form)))
      retValue <- formulation$round_function(amt)
      if (retValue > doseCap_) {
        retValue <- doseCap_
      }
      return(retValue)
    })
    regimen$AMT[iterationRows] <- roundedAmt[iterationRows] # Rounded amounts only

    # Filling zeroes
    if (cautious) {
      lastAmount <- getLastAmountInFuture(regimen=regimen)
      regimen <- fillTrailingZeroesInFuture(regimen=regimen, replacement=lastAmount)
    }

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
