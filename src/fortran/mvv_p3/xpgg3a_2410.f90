module xpgg3a_2410
  use iso_c_binding
!  character(len=*), parameter :: name_xpgg3 = "xpgg3a"
contains
!
! ..File: xPgg3a.f   
!
! ..The 4-loop MSbar gluon-gluon splitting function P_{gg}^(3) 
!    for the  evolution of the flavour-singlet parton distributions.
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
! ..These are approximations for fixed nf = 3, 4 and 5 based on the 
!    first 10 even moments together with small-x/large-x constraints.
!    The two sets providing the error estimate are called via IMOD = 1 
!    and IMOD = 2.  Any other value of IMOD invokes their average.
!
! ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to 
!    the kernel superscripts [2], [3], and [1] in that equation.
!   
! ..Reference: Four-loop splitting functions in QCD 
!                - The gluon-gluon case -
!
!              G. Falcioni, F. Herzog, S. Moch, A. Pelloni and A. Vogt
!              ZU-TH 47/24, DESY-24-144, LTH 1384
!
! =====================================================================
!
! ..The regular part
!
       FUNCTION P3GGA_2410 (Y, NF, IMOD) bind(c, name="p3gga_2410")
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf,nf2,nf3
!
       YM  = 1.D0/Y
       Y1  = 1.D0-Y
       DL  = DLOG(Y)
       DL1 = DLOG(1.D0-Y)
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
!  Known large-x coefficients [except delta(1-x)]
!
       A4gluon =  40880.330D0     - 11714.246D0*nf &
     &          + 440.04876D0*nf2 + 7.3627750D0*nf3
       Ccoeff  = 8.5814120D+4 - 1.3880515D+4*nf + 1.3511111D+2*nf2
       Dcoeff  = 5.4482808D+4 - 4.3411337D+3*nf - 2.1333333D+1*nf2 
!
       x1L4cff = 5.6460905D+1*nf - 3.6213992D+0*nf2
       x1L3cff = 2.4755054D+2*nf - 4.0559671D+1*nf2 + 1.5802469D+0*nf3
!    
!  Known small-x coefficients 
!
       bfkl0  = - 8.3086173D+3
       bfkl1  = - 1.0691199D+5 - 9.9638304D+2*nf
!
       x0L6cff =  1.44D+2 - 2.7786008D+1*nf + 7.9012346D-1*nf2
       x0L5cff = -1.44D+2 - 1.6208066D+2*nf + 1.4380247D+1*nf2
       x0L4cff =  2.6165784D+4     - 3.3447551D+3*nf &
     &          + 9.1522635D+1*nf2 - 1.9753086D-1*nf3
!
!  The resulting part of the function
!
       P3gg01 =                 &
     &      + bfkl0* DL**3*YM   &
     &      + bfkl1* DL**2*YM   &
     &      + x0L6cff* DL**6    &
     &      + x0L5cff* DL**5    &
     &      + x0L4cff* DL**4    &
     &      + Ccoeff* DL1       &
     &      + Dcoeff - A4gluon  &
     &      + x1L4cff* Y1*DL1**4&
     &      + x1L3cff* Y1*DL1**3
!
!  The selected approximations for nf = 3, 4, 5
!
       IF ( NF .EQ. 3 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 421311.0D0 * Y1*DL*YM  &
     &      - 325557.0D0 * Y1*YM     &
     &      + 1679790.0D0* Y1        &
     &      - 1456863.0D0* Y1*Y      &
     &      + 3246307.0D0* Y1*DL     &
     &      + 2026324.0D0* DL*DL     &
     &      + 549188.0D0 * DL**3     &
     &      +   8337.0D0 * Y1*DL1    &
     &      +  26718.0D0 * Y1*DL1*DL1&
     &      -  27049.0D0 * Y1*Y1*DL1**3
         P3ggApp2 = P3gg01            &   
     &      - 700113.0D0 * Y1*DL*YM      &
     &      - 2300581.0D0* Y1*YM         &
     &      + 896407.0D0 * Y1*(1.0D0+2.0D0*Y)  &
     &      - 162733.0D0 * Y1*Y*Y        &
     &      - 2661862.0D0* Y1*DL         &
     &      + 196759.0D0 * DL*DL         &
     &      - 260607.0D0 * DL**3         &
     &      +  84068.0D0 * Y1*DL1        &
     &      + 346318.0D0 * Y1*DL1*DL1    &
     &      + 315725.0D0 * DL*DL1*DL1
       ELSE IF ( NF .EQ. 4 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 437084.0D0 * Y1*DL*YM  &
     &      - 361570.0D0 * Y1*YM     &
     &      + 1696070.0D0* Y1        &
     &      - 1457385.0D0* Y1*Y      &
     &      + 3195104.0D0* Y1*DL     &
     &      + 2009021.0D0* DL*DL     &
     &      + 544380.0D0 * DL**3     &
     &      +  9938.0D0  * Y1*DL1    &
     &      +  24376.0D0 * Y1*DL1*DL1&
     &      -  22143.0D0 * Y1*Y1*DL1**3
         P3ggApp2 = P3gg01           &
     &      - 706649.0D0 * Y1*DL*YM     &
     &      - 2274637.0D0* Y1*YM        &
     &      + 836544.0D0 * Y1*(1.0D0+2.0D0*Y) &
     &      - 199929.0D0 * Y1*Y*Y       &
     &      - 2683760.0D0* Y1*DL        &
     &      + 168802.0D0 * DL*DL        &
     &      - 250799.0D0 * DL**3        &
     &      +  36967.0D0 * Y1*DL1       &
     &      +  24530.0D0 * Y1*DL1*DL1   &
     &      -  71470.0D0 * Y1*Y1*DL1*DL1
       ELSE IF ( NF .EQ. 5 ) THEN
         P3ggApp1 = P3gg01        &
     &      - 439426.0D0 * Y1*DL*YM  &
     &      - 293679.0D0 * Y1*YM     &
     &      + 1916281.0D0* Y1        &
     &      - 1615883.0D0* Y1*Y      &
     &      + 3648786.0D0* Y1*DL     &
     &      + 2166231.0D0* DL*DL     &
     &      + 594588.0D0 * DL**3     &
     &      +  50406.0D0 * Y1*DL1    &
     &      +  24692.0D0 * Y1*DL1*DL1&
     &      + 174067.0D0 * Y1*Y1*DL1
         P3ggApp2 = P3gg01           &
     &      - 705978.0D0 *  Y1*DL*YM    &
     &      - 2192234.0D0* Y1*YM        &
     &      + 1730508.0D0* Y1*Y         &
     &      + 353143.0D0 * Y1*(2.0D0-Y*Y)  &
     &      - 2602682.0D0* Y1*DL        &
     &      + 178960.0D0 * DL*DL        &
     &      - 218133.0D0 * DL**3        &
     &      +   2285.0D0 * Y1*DL1       &
     &      +  19295.0D0 * Y1*DL1*DL1   &
     &      -  13719.0D0 * Y1*Y1*DL1*DL1
      ELSE IF ( NF .EQ. 6 ) THEN
         P3ggApp1 = P3gg01              &
     &      - 476018.0D0 * Y1*DL*YM        &
     &      - 469289.0D0 * Y1*YM           &
     &      + 2049351.0D0 * Y1             &
     &      - 1589000.0D0 * Y1*Y           &
     &      + 3185549.0D0 * Y1*DL          &
     &      + 1994521.0D0 * DL*DL          &
     &      + 527723.0D0 * DL**3           &
     &      - 340674.0D0 * Y1*DL1          &
     &      +  22460.0D0 * Y1*DL1*DL1      &
     &      - 394556.0D0 * DL*DL1             
         P3ggApp2 = P3gg01              &
     &      - 709863.0D0 *  Y1*DL*YM       &
     &      - 2134347.0D0* Y1*YM           &
     &      + 1605315.0D0* Y1*Y            &
     &      + 360743.0D0 * Y1*(2.0D0-Y*Y)     &
     &      - 2426250.0D0* Y1*DL           &
     &      + 230631.0D0 * DL*DL           &
     &      - 185804.0D0 * DL**3           &
     &      - 7992.9D0 * Y1*DL1           &
     &      + 15918.0D0 * Y1*DL1*DL1       &
     &      - 32771.0D0 * Y1*Y1*DL1        
       ELSE
         WRITE(6,*) '  Error in function P3ggA: choice of nf   '
         CALL ABORT
       END IF
!
!  We return one of the two error-band representatives or the present 
!  best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GGA_2410 = P3ggApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GGA_2410 = P3ggApp2
       ELSE
         P3GGA_2410 = 0.5D0 * ( P3ggApp1 + P3ggApp2 )
       END IF
!
       RETURN
       END
!
! ---------------------------------------------------------------------
!
! ..The singular (soft) piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGB_2410 (Y, NF, IMOD) bind(c, name="p3ggb_2410")
       IMPLICIT REAL*8 (A - Z)
       INTEGER nf, IMOD
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
       A4gluon =  40880.330D0     - 11714.246D0*nf &
            &          + 440.04876D0*nf2 + 7.3627750D0*nf3
       P3GGB_2410  = A4gluon
!
       RETURN
       END
!
! ---------------------------------------------------------------------
!
!
! ..The 'local' piece of P_gg^(3).
!   Note: A4gluon is provided by a common block set in P3GGA
!
       FUNCTION P3GGC_2410 (Y, NF, IMOD) bind(c, name="p3ggc_2410")
!
       IMPLICIT REAL*8 (A - Z)
       INTEGER IMOD, nf,nf2,nf3
!
       nf2 = nf*nf
       nf3 = nf*nf2
!
!  The coefficient of delta(1-x), also called the virtual anomalous
!  dimension. nf^0 and nf^1 are still approximate, but the error at 
       !  nf^1 is far too small to be relevant in this context.

       B4gluon = 68587.64D0 - 18143.983D0*nf + 423.81135D0*nf2 + 9.0672154D-1*nf3
       IF ( IMOD .EQ. 1 ) B4gluon = B4gluon - 0.2
       IF ( IMOD .EQ. 2 ) B4gluon = B4gluon + 0.2
!
       P3GGC_2410 = B4gluon
!
       RETURN
       END
!
! =================================================================av==
end module xpgg3a_2410
