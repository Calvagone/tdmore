library(tdmore)

rxode2::rxode2({
  TVV1 = 24.4
  TVV2 = 7.01
  TVQ = 4.97
  TVCL = 9.87
  ECL = 0 # ECL Fixed to 0. Not considered as a parameter in TDMore.
  CL = TVCL * exp(ECL)
  V1 = TVV1 * exp(EV1)
  V2 = TVV2
  Q = TVQ
  K12 = Q/V1
  K21 = Q/V2

  d/dt(center) = - CL/V1 * center - K12*center + K21 * periph
  d/dt(periph) = K12*center - K21 * periph

  CONC = center / V1
}) %>%
  # EV1 corresponds to OMEGA_V1 (54% SD)
  tdmore(omega=c(EV1=0.287),
         res_var=list(errorModel(var="CONC", prop=0.371))) # Proportional error (37% SD)
