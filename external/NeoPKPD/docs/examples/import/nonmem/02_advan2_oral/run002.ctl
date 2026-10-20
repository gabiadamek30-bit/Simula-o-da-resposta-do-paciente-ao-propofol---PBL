$PROBLEM One-compartment oral PK model with first-order absorption
$INPUT ID TIME DV AMT MDV EVID
$DATA ../data.csv IGNORE=@

$SUBROUTINES ADVAN2 TRANS2

$PK
  TVKA = THETA(1)
  TVCL = THETA(2)
  TVV  = THETA(3)

  KA = TVKA * EXP(ETA(1))
  CL = TVCL * EXP(ETA(2))
  V  = TVV  * EXP(ETA(3))

  S2 = V

$ERROR
  IPRED = F
  W = SQRT(THETA(4)**2 + (THETA(5)*IPRED)**2)
  IF(W.EQ.0) W = 1
  IRES = DV - IPRED
  IWRES = IRES / W
  Y = IPRED + W * ERR(1)

$THETA
  (0, 1.5)    ; Ka (1/h)
  (0, 5.0)    ; CL (L/h)
  (0, 50.0)   ; V (L)
  (0, 0.5)    ; Additive error (mg/L)
  (0, 0.10)   ; Proportional error

$OMEGA
  0.16        ; IIV Ka (CV ~40%)
  0.09        ; IIV CL (CV ~30%)
  0.0625      ; IIV V (CV ~25%)

$SIGMA
  1 FIX       ; Combined residual

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 PRINT=10 NOABORT
$COVARIANCE PRINT=E
