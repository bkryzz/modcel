!!!http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.35.shtml

        subroutine voron2()

        common / block1 / rx, ry

!    *******************************************************************
!    ** CONSTRUCTION OF THE VORONOI POLYGON IN 2D.                    **
!    **                                                               **
!    ** THIS PROGRAM TAKES IN A CONFIGURATION IN A SQUARE BOX WITH    **
!    ** CONVENTIONAL PERIODIC BOUNDARY CONDITIONS AND FOR EACH ATOM   **
!    ** OBTAINS THE SURROUNDING VORONOI POLYGON, DEFINED AS THAT      **
!    ** REGION OF SPACE CLOSER TO THE CHOSEN ATOM THAN TO ANY OTHER.  **
!    ** NEIGHBOURING POLYGONS DEFINE NEIGHBOURING ATOMS.              **
!    ** THE PROGRAM IS SLOW BUT ESSENTIALLY FOOLPROOF.                **
!    ** WE USE THE MINIMUM IMAGE CONVENTION AND SET A CUTOFF BEYOND   **
!    ** WHICH ATOMS ARE ASSUMED NOT TO BE NEIGHBOURS: BOTH OF THESE   **
!    ** MEASURES ARE DANGEROUS FOR SMALL AND/OR RANDOM SYSTEMS.       **
!    ** WE DELIBERATELY DO NOT USE PREVIOUSLY-FOUND NEIGHBOURS IN     **
!    ** CONSTRUCTING NEIGHBOUR LISTS, SO THAT AN INDEPENDENT CHECK    **
!    ** MAY BE MADE AT THE END.                                       **
!    ** HERE WE SIMPLY PRINT OUT THE GEOMETRICAL INFORMATION AT THE   **
!    ** END.  THE OUTPUT IS QUITE LENGTHY.  IN PRACTICE, IT WOULD     **
!    ** PROBABLY BE ANALYZED DIRECTLY WITHOUT PRINTING OUT.           **
!    ** NB: BEWARE DEGENERATE CONFIGURATIONS, I.E. ONES IN WHICH MORE **
!    ** THAN THREE VORONOI DOMAINS SHARE A VERTEX. THE SQUARE LATTICE **
!    ** IS AN EXAMPLE.                                                **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                        NUMBER OF ATOMS              **
!    ** REAL    RX(N),RY(N)              POSITIONS                    **
!    ** REAL    PX(MAXCAN),PY(MAXCAN)    CANDIDATE RELATIVE POSITIONS **
!    ** REAL    PS(MAXCAN)               SQUARED RELATIVE DISTANCES   **
!    ** INTEGER NVER                     NUMBER OF VERTICES FOUND     **
!    ** INTEGER NEDGE                    NUMBER OF EDGES FOUND        **
!    ** INTEGER VERTS(MAXCAN)            VERTICES FOR EACH CANDIDATE  **
!    **                                  = 0 IF NOT A NEIGHBOUR       **
!    **                                  = 2 ( 1 EDGE ) IF NEIGHBOUR  **
!    ** REAL    RXVER(MAXVER)            VERTEX RELATIVE X-COORD      **
!    ** REAL    RYVER(MAXVER)            VERTEX RELATIVE Y-COORD      **
!    ** INTEGER IVER(MAXVER)             ATOMIC INDICES TAGGING       **
!    ** INTEGER JVER(MAXVER)             .. EACH VERTEX OF POLYGON    **
!    **                                                               **
!    ** ROUTINES REFERENCED:                                          **
!    **                                                               **
!    ** SUBROUTINE READCN ( CNFILE, N, BOX )                          **
!    **    READS IN CONFIGURATION, NUMBER OF ATOMS, BOX SIZE          **
!    ** SUBROUTINE SORT2D ( MAXCAN, PX, PY, PS, TAG, NCAN )           **
!    **    SORTS NEIGHBOUR DETAILS INTO ASCENDING DISTANCE ORDER      **
!    ** SUBROUTINE WORK2D ( MAXCAN, MAXVER, NCAN, NVER, NEDGE,        **
!    **     PX, PY, PS, VERTS, RXVER, RYVER, IVER, JVER )             **
!    **    CARRIES OUT THE VORONOI CONSTRUCTION                       **
!    *******************************************************************

        INTEGER     MAXN, MAXCAN, MAXVER
        PARAMETER ( MAXN = 108, MAXCAN = 50, MAXVER = 50 )

        REAL        RX(MAXN), RY(MAXN)

        REAL        PX(MAXCAN), PY(MAXCAN), PS(MAXCAN)
        INTEGER     TAG(MAXCAN), VERTS(MAXCAN)

        REAL        RXVER(MAXVER), RYVER(MAXVER)
        INTEGER     IVER(MAXVER), JVER(MAXVER)
        INTEGER     NABLST(MAXVER,MAXN), NNAB(MAXN), INAB, JNAB

        INTEGER     NCAN, NVER, NCOORD, NEDGE
        INTEGER     I, J, CAN, VER, N
        REAL        BOX, BOXINV, RCUT, RCUTSQ, COORD
        REAL        RXJ, RYJ, RZJ, RXIJ, RYIJ, RZIJ, RIJSQ
        CHARACTER   CNFILE*30
        LOGICAL     OK

!    *******************************************************************

        WRITE(*,'(1H1,'' **** PROGRAM VORON2 ****                  '')')
        WRITE(*,'(//1X,''VORONOI CONSTRUCTION IN 2D                '')')

!    ** BASIC PARAMETERS **

        WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
        READ (*,'(A)') CNFILE
        WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE

!    ** READCN MUST READ IN INITIAL CONFIGURATION  **

        CALL READCN ( CNFILE, N, BOX )

        WRITE(*,'(1X,I5,''-ATOM CONFIGURATION'')') N
        WRITE(*,'('' BOX LENGTH = '',F10.5)') BOX
        WRITE(*,'('' ENTER NEIGHBOUR CUTOFF IN SAME UNITS '')')
        READ (*,*) RCUT
        WRITE(*,'('' NEIGHBOUR CUTOFF = '',F10.5)') RCUT

        RCUTSQ = RCUT ** 2
        BOXINV = 1.0 / BOX

!    ** ZERO ACCUMULATORS **

        DO 100 J = 1, N

           NNAB(J) = 0

           DO 90 INAB = 1, NVER

              NABLST(INAB,J) = 0

 90         CONTINUE

 100     CONTINUE

!    *******************************************************************
!    ** MAIN LOOP STARTS                                              **
!    *******************************************************************

        DO 1000 J = 1, N

           IF ( MOD ( J, 2 ) .EQ. 0 ) THEN

              WRITE(*,'(///1X,''RESULTS FOR ATOM '',I5)') J

           ELSE

              WRITE(*,'(1H1,''RESULTS FOR ATOM '',I5)') J

           ENDIF

           RXJ = RX(J)
           RYJ = RY(J)
           CAN = 0

!       ** SELECT CANDIDATES **

           DO 500 I = 1, N

              IF ( I .NE. J ) THEN

                 RXIJ = RX(I) - RXJ
                 RYIJ = RY(I) - RYJ
                 RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                 RIJSQ  = RXIJ ** 2 + RYIJ ** 2

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    CAN = CAN + 1

                    IF ( CAN .GT. MAXCAN ) THEN

                       WRITE(*,'('' TOO MANY CANDIDATES '')')
                       STOP

                    ENDIF

                    PX(CAN)  = RXIJ
                    PY(CAN)  = RYIJ
                    PS(CAN)  = RIJSQ
                    TAG(CAN) = I

                 ENDIF

              ENDIF

 500        CONTINUE

!       ** CANDIDATES HAVE BEEN SELECTED **

           NCAN = CAN

!       ** SORT INTO INCREASING DISTANCE ORDER **
!       ** THIS SHOULD IMPROVE EFFICIENCY      **

           CALL SORT2D ( MAXCAN, PX, PY, PS, TAG, NCAN )

!       ** PERFORM VORONOI CONSTRUCTION **

           CALL WORK2D ( MAXCAN, MAXVER, NCAN, NVER, NEDGE,
     :                 PX, PY, PS, VERTS,
     :                 RXVER, RYVER, IVER, JVER )

!       ** WRITE OUT RESULTS **

           WRITE(*,'(/1X,''NUMBER OF NEIGHBOURS '',I5)') NEDGE
           WRITE(*,'(/1X,''NEIGHBOUR LIST '')')
           WRITE(*,10001)

           DO 800 CAN = 1, NCAN

              IF ( VERTS(CAN) .NE. 0 ) THEN

                 PS(CAN) = SQRT ( PS(CAN) )
                 WRITE(*,'(1X,I5,3X,I5,3X,2F12.5,3X,F12.5)')
     :              TAG(CAN), VERTS(CAN), PX(CAN), PY(CAN), PS(CAN)
                 NNAB(J) = NNAB(J) + 1
                 NABLST(NNAB(J),J) = TAG(CAN)

              ENDIF

 800        CONTINUE

           WRITE(*,'(/1X,''NUMBER OF VERTICES   '',I5)') NVER
           WRITE(*,'(/1X,''VERTEX LIST '')')
           WRITE(*,10002)

           DO 900 VER = 1, NVER

              WRITE(*,'(1X,2I5,3X,2F12.5)')
     :           TAG(IVER(VER)), TAG(JVER(VER)),
     :           RXVER(VER), RYVER(VER)

 900        CONTINUE

 1000    CONTINUE

!    *******************************************************************
!    ** MAIN LOOP ENDS                                                **
!    *******************************************************************

        WRITE(*,'(1H1,''FINAL SUMMARY'')')
        WRITE(*,10003)

        NCOORD = 0

        DO 2000 J = 1, N

           NCOORD = NCOORD + NNAB(J)

           WRITE(*,'(1X,I5,3X,I5,3X,30I3)') J, NNAB(J),
     :       ( NABLST(INAB,J), INAB = 1, NNAB(J) )

!       ** CHECK THAT IF I IS A NEIGHBOUR OF J **
!       ** THEN J IS ALSO A NEIGHBOUR OF I     **

           DO 1500 INAB = 1, NNAB(J)

              I = NABLST(INAB,J)

              OK = .FALSE.
              JNAB = 0

 1200          IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN

                 OK = ( J .EQ. NABLST(JNAB,I) )
                 JNAB = JNAB + 1
                 GOTO 1200

              ENDIF

              IF ( .NOT. OK ) THEN

                 WRITE(*,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') J, I

              ENDIF

 1500       CONTINUE

 2000    CONTINUE

        COORD = REAL ( NCOORD ) / REAL ( N )

        WRITE(*,'(/1X,'' AVERAGE COORDINATION NUMBER = '',F10.5)') COORD

        STOP

 10001   FORMAT(/1X,'ATOM ',3X,'EDGE ',
     :         /1X,'INDEX',3X,'VERTS',3X,
     :         '      RELATIVE POSITION   ',3X,'  DISTANCE  ')
 10002   FORMAT(/1X,'   INDICES         RELATIVE POSITION ')
 10003   FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')

        END

!===================================================================

        SUBROUTINE READCN ( CNFILE, N, BOX )

        COMMON / BLOCK1 / RX, RY

!    *******************************************************************
!    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION                   **
!    *******************************************************************

        INTEGER     MAXN
        PARAMETER ( MAXN = 108 )

        REAL        RX(MAXN), RY(MAXN), BOX
        INTEGER     N

        CHARACTER   CNFILE*(*)

        INTEGER     CNUNIT, I
        PARAMETER ( CNUNIT = 10 )

!    *******************************************************************

        OPEN ( UNIT = CNUNIT, FILE = CNFILE,
     :         STATUS = 'OLD', FORM = 'UNFORMATTED' )

        READ ( CNUNIT ) N, BOX
        IF ( N .GT. MAXN ) STOP ' N TOO LARGE '
        READ ( CNUNIT ) ( RX(I), I = 1, N ), ( RY(I), I = 1, N )

        CLOSE ( UNIT = CNUNIT )

        RETURN
        END

!=======================================================================


        SUBROUTINE WORK2D ( MAXCAN, MAXV, NN, NV, NE, RX, RY, RS, VERTS,
     :                    VX, VY, IV, JV )

!    *******************************************************************
!    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
!    **                                                               **
!    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
!    ** I.E. THE EDGES OF THE POLYGON GO THROUGH THE POINTS.          **
!    *******************************************************************

        INTEGER     MAXCAN, NN, MAXV, NV, NE
        INTEGER     VERTS(MAXCAN)
        REAL        RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
        REAL        VX(MAXV), VY(MAXV)
        INTEGER     IV(MAXV), JV(MAXV)

        LOGICAL     OK
        INTEGER     I, J, L, NN1, N, V
        REAL        AI, BI, CI, AJ, BJ, CJ, DET, DETINV
        REAL        VXIJ, VYIJ
        REAL        TOL
        PARAMETER ( TOL = 1.E-6 )

!    *******************************************************************

!    ** IF THERE ARE LESS THAN 3 POINTS GIVEN **
!    ** WE CANNOT CONSTRUCT A POLYGON         **

        IF ( NN .LT. 3 ) THEN

           WRITE(*,'('' LESS THAN 3 POINTS GIVEN TO WORK '',I5)') NN
           STOP

        ENDIF

        NN1 = NN - 1
        V = 0

!    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
!    ** DEFINED BY THE INTERSECTION OF 2 EDGES  **
!    ** EACH EDGE IS DEFINED BY RX,RY,RS.       **

        DO 400 I = 1, NN1

           AI =  RX(I)
           BI =  RY(I)
           CI = -RS(I)

           DO 300 J = I + 1, NN

              AJ =  RX(J)
              BJ =  RY(J)
              CJ = -RS(J)

              DET = AI * BJ - AJ * BI

              IF ( ABS ( DET ) .GT. TOL ) THEN

!             ** THE EDGES INTERSECT **

                 DETINV = 1.0 / DET

                 VXIJ = ( BI * CJ - BJ * CI ) * DETINV
                 VYIJ = ( AJ * CI - AI * CJ ) * DETINV

!             ** NOW WE TAKE SHOTS AT THE VERTEX **
!             ** USING THE REMAINING EDGES ..... **

                 OK = .TRUE.
                 L  = 1

 100              IF ( OK .AND. ( L .LE. NN ) ) THEN

                    IF ( ( L .NE. I ) .AND. ( L .NE. J ) ) THEN

                       OK = ( RX(L) * VXIJ + RY(L) * VYIJ ) .LE. RS(L)

                    ENDIF

                    L = L + 1
                    GOTO 100

                 ENDIF

!             ** IF THE VERTEX MADE IT      **
!             ** ADD IT TO THE HALL OF FAME **
!             ** CONVERT TO CORRECT SCALE   **

                 IF ( OK ) THEN

                    V = V + 1
                    IF ( V .GT. MAXV ) STOP 'TOO MANY VERTICES'
                    IV(V)  = I
                    JV(V)  = J
                    VX(V) = 0.5 * VXIJ
                    VY(V) = 0.5 * VYIJ

                 ENDIF

              ENDIF

 300        CONTINUE

 400     CONTINUE

!    ** THE SURVIVING VERTICES DEFINE THE VORONOI POLYGON **

        NV = V

        IF ( NV .LT. 3 ) THEN

           WRITE(*,'('' LESS THAN 3 VERTICES FOUND IN WORK '',I5)') NV
           STOP

        ENDIF

!    ** IDENTIFY NEIGHBOURING POINTS **

        DO 500 N = 1, NN

           VERTS(N) = 0

500     CONTINUE

        DO 600 V = 1, NV

           VERTS(IV(V)) = VERTS(IV(V)) + 1
           VERTS(JV(V)) = VERTS(JV(V)) + 1

600     CONTINUE

!    ** POINTS WITH NONZERO VERTS ARE NEIGHBOURS **
!    ** IF NONZERO, VERTS SHOULD BE EQUAL TO 2   **

!    ** CHECK RESULT AND COUNT EDGES **

        OK = .TRUE.
        NE = 0

        DO 700 N = 1, NN

           IF ( VERTS(N) .GT. 0 ) THEN

              NE = NE + 1

              IF ( VERTS(N) .NE. 2 ) THEN

                 OK = .FALSE.

              ENDIF

           ENDIF

 700     CONTINUE

        IF ( .NOT. OK ) THEN

           WRITE (*,'('' **** VERTEX ERROR: DEGENERACY ? **** '')')

        ENDIF

        IF ( NE .NE. NV ) THEN

           WRITE(*,'('' **** EDGE   ERROR: DEGENERACY ? ****  '')')

        ENDIF

        RETURN
        END

!=========================================================================

        SUBROUTINE SORT2D ( MAXCAN, RX, RY, RS, TAG, NN )

!    *********************************************************************
!    ** ROUTINE TO SORT2D NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
!    **                                                                 **
!    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.      **
!    *********************************************************************

        INTEGER MAXCAN, NN
        REAL    RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
        INTEGER TAG(MAXCAN)

        LOGICAL CHANGE
        INTEGER I, ITOP, I1, TAGI
        REAL    RXI, RYI, RSI

!    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NN - 1

 1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN

                 RXI = RX(I)
                 RYI = RY(I)
                 RSI = RS(I)
                 TAGI = TAG(I)

                 RX(I) = RX(I1)
                 RY(I) = RY(I1)
                 RS(I) = RS(I1)
                 TAG(I) = TAG(I1)

                 RX(I1) = RXI
                 RY(I1) = RYI
                 RS(I1) = RSI
                 TAG(I1) = TAGI

                 CHANGE = .TRUE.

              ENDIF

 100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

        RETURN
        END


