$PROBLEM THEOPHYLLINE PK - ORAL ONE COMPARTMENT
; Migration example for NeoPKPD

$DATA theo_sd.csv IGNORE=@
$INPUT ID TIME DV AMT EVID CMT WT

$SUBROUTINE ADVAN2 TRANS2

$PK
; Weight-based covariate model (allometric scaling)
TVCL = THETA(1) * (WT/70)**0.75
TVV  = THETA(2) * (WT/70)
TVKA = THETA(3)

; Individual parameters with IIV
CL = TVCL * EXP(ETA(1))
V  = TVV  * EXP(ETA(2))
KA = TVKA * EXP(ETA(3))

; Scaling
S2 = V

$ERROR
IPRED = F
Y = F * (1 + ERR(1))

$THETA
(0, 2.8)    ; CL (L/h) - typical for 70 kg
(0, 35)     ; V (L) - typical for 70 kg
(0, 1.5)    ; KA (1/h)

$OMEGA
0.09        ; IIV CL (CV ~30%)
0.04        ; IIV V (CV ~20%)
0.16        ; IIV KA (CV ~40%)

$SIGMA
0.04        ; Proportional residual error (CV ~20%)

$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=10
$COVARIANCE
$TABLE ID TIME IPRED PRED DV CWRES ETA1 ETA2 ETA3
       NOPRINT ONEHEADER FILE=sdtab001
