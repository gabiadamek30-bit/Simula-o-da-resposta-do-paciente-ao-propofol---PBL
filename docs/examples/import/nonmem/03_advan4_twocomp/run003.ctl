$PROBLEM Two-compartment IV bolus PK model
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  TVQ  = THETA(3)
  TVV2 = THETA(4)

  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  Q  = TVQ
  V2 = TVV2

  S1 = V1

$ERROR
  IPRED = F
  W = IPRED * THETA(5)
  IF(W.EQ.0) W = 1
  IRES = DV - IPRED
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)

$THETA
  (0, 5.0)    ; CL (L/h)
  (0, 10.0)   ; V1 (L)
  (0, 2.0)    ; Q (L/h)
  (0, 20.0)   ; V2 (L)
  (0, 0.15)   ; Proportional error

$OMEGA
  0.09        ; IIV CL (CV ~30%)
  0.0625      ; IIV V1 (CV ~25%)

$SIGMA
  1 FIX       ; Proportional residual

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=10 NOABORT
$COVARIANCE PRINT=E
