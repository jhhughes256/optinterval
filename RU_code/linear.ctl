$PROBLEM An example of linear interpolation of observed concentrations

$INPUT ID	TIME	COBS	INT	SLOP DV

$DATA lineardata.csv IGNORE=C

$SUBROUTINES ADVAN6 TOL=5

$MODEL
COMP=(EFFECT)

$PK
B = INT                     ;Reassign intercept
M = SLOP                    ;Reassign slope
KEO = THETA(1)*EXP(ETA(1))  ;Effect compartment rate constant

$DES
CTMP = B+(M*T) ; temporary concentration interpolated from supplied slopes and intercepts
DADT(1) = KEO*(CTMP-A(1)) ;differential equation for effect compartment


$ERROR
 CINT = B+(M*TIME)  ;recalculated as differential equations carry last value over into next subject
 CEFF = A(1)        ;recalculated as differential equations carry last value over into next subject
 Y=CEFF+EPS(1)

$THETA  0.2 ;KEOPOP
$OMEGA  0.05
$SIGMA  0.01

$SIMULATION (9215690) ; seed 1-7 digits

$TABLE ID TIME COBS INT SLOP CINT CEFF DV NOPRINT FILE= NOHEADER

