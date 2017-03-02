$PROBLEM An example of cubic spline interpolation of observed concentrations

$INPUT ID	TIME	COBS	KNOT	COF1	COF2	COF3	COF4	DV

$DATA cubicdata.int.csv IGNORE=C

$SUBROUTINES ADVAN6 TOL=5

$MODEL
COMP=(EFFECT)

$PK
CO1 = COF1    ;Reassign spline coefficients
CO2 = COF2
CO3 = COF3
CO4 = COF4

KEO = THETA(1)*EXP(ETA(1)) ;Effect compartment rate constant


$DES
TINT = T - TIME     ;TINT is time in interval between knots (starting at zero)
                    ;to be consistent with the way the coefficients are calculated by interpSpline in R
CCUB = CO1+CO2*TINT+CO3*(TINT**2)+CO4*(TINT**3) ; time varying cubic spline interpolated from supplied coefficients

DADT(1) = KEO*(CCUB-A(1))    ;differential equation for effect compartment


$ERROR
  CEFF=A(1)
  Y=CEFF + EPS(1)

$THETA  0.2  ;KEOPOP
$OMEGA  0.05
$SIGMA  0.01

$SIMULATION (9215690)

$TABLE ID TIME TINT COBS CCUB CEFF DV NOPRINT FILE= NOHEADER

