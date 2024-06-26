

#' @rdname posthoc
#' @param .obs Either a string or NULL. If a string, the output will contain a column
#' with that name, describing how many observations were taken into account
#' @param .models List of models to use per occasion (index 1=No OCC, index 2=OCC1, index 3=OCC2, etc.)
#' @section Proseval:
#' The proseval tool calculates a prospective evaluation, and therefore creates multiple outputs per input row.
#' If N is the number of observations, it will calculate N+1 fits.
#' The first fit will be the population fit (not using any observations).
#' The next fit will only use one observation.
#' The next fit will use two observations.
#' And so on.
#'
#' Please note that ipred incorporates all time-varying covariates, even if these may be considered to be "in the future".
#' The 'par' argument in dataTibble can optionally be a list, with a starting value for each iteration.
#' @export
proseval2 <- function(x, ..., .fit="fit", .prediction="ipred", .elapsed="elapsed", .obs="OBS", .models, .window) {
  if(is.null(.fit)) stop(".fit column needs to be present for proseval to work")
  if(tibble::is_tibble(x)) {
    res <- dplyr::rowwise(x) %>%
      dplyr::do({proseval2(., ..., .fit=.fit, .prediction=.prediction, .elapsed=.elapsed, .obs=.obs, .models=.models, .window=.window)})
    return(res)
  } else {
    stopifnot(is.list( x ))

    output <- list()
    for(i in c(0, seq_len(nrow(x$observed))) ){
      argsPosthoc <- x
      if(!is.null(.models)) {
        # tmp <- tar_read(base_model_ctrl_iov_data_tibble)[1, ]
        # x <- list()
        # x$observed <- tmp$observed[[1]]
        # x$regimen <- tmp$regimen[[1]]
        # x$covariates <- tmp$covariates[[1]]
        if (i==0) {
          lastObsTime <- 0
        } else {
          lastObsTime <- x$observed[i, "TIME"]
        }
        # Covariates: DAY either corresponds to an observation or an evid-2 row (both are at time 6am)
        # So we need to look at the day just before
        day <- x$covariates %>%
          dplyr::filter(TIME < lastObsTime) %>%
          dplyr::pull(DAY) %>%
          dplyr::last()
        if (i==0) {
          day <- 0
        }
        modelIndex <- day + 1
        if (modelIndex > .models %>% length()) {
          modelIndex <- .models %>% length()
        }
        cat(paste0("Using model no.=", modelIndex))
        argsPosthoc$object <- .models[[modelIndex]] # Overriding default model
      }
      if(i == 0) {
        argsPosthoc['observed'] <- list(NULL)
      } else {
        if (is.null(.window)) {
          indices <- seq_len(i)
        } else {
          indices <- seq(i - .window + 1, i)
          indices <- indices[indices >= 1]
        }
        argsPosthoc$observed <- x$observed[indices, ]
      }
      if(is.list(x$par)) {
        argsPosthoc$par <- if(i>0) x$par[[i]] else NULL
      }
      res <- posthoc(argsPosthoc, ..., .fit=.fit, .prediction=NULL, .elapsed=.elapsed)
      fit <- res[[.fit]][[1]]
      res[["observed"]] <- list(x$observed) #always use the same OBSERVED data.frame
      if(!is.null(.prediction)) res[[.prediction]] <- list( stats::predict(fit, x$observed) )
      if(!is.null(.obs)) res[[.obs]] <- i
      output[[i+1]] <- res
    }
    return(bind_rows(output))
  }
}
