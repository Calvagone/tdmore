library(tdmore)

rxode2::rxode2({
  TVKA = 3.7
  TVV1 = 61
  TVQ = 10
  TVCL = 3.7

  KA = TVKA
  CL = TVCL * (WT/70)^0.75 * exp(ECL)
  V1 = TVV1 * exp(EV1)
  V2 = V1
  Q = TVQ
  K12 = Q/V1
  K21 = Q/V2

  d/dt(depot) = -KA*depot
  d/dt(center) = KA*depot - CL/V1 * center - K12*center + K21 * periph
  d/dt(periph) = K12*center - K21 * periph

  CONC = center / V1
}) %>%
tdmore(omega=c(ECL=0.0784, EV1=0.0361),
       res_var=list(errorModel(var="CONC", prop=0.23))) %>%
metadata(covariate(name="WT", label="Weight", unit="kg", min=15, max=300)) %>%
metadata(output(name="CONC", label="Concentration", unit="mg/L")) %>%
metadata(formulation(name="Drug", unit="mg", dosing_interval=8, default_value=5, round_function=function(x){round(x*2)/2})) %>%
metadata(target(min=12, max=20))


