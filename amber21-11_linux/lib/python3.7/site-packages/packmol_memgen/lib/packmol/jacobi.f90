!
! JACOBI
! Jacobi diagonalizer with sorted output.  Same calling sequence as
! EISPACK routine, but must specify nrot!
!
      SUBROUTINE jacobi (a, n, np, d, v, nrot)
      IMPLICIT CHARACTER (A-Z)
!
      INTEGER n, np, nrot
      DOUBLEPRECISION a (np, n)
      DOUBLEPRECISION d (n)
      DOUBLEPRECISION v (np, n)
!
      DOUBLEPRECISION onorm, dnorm
      DOUBLEPRECISION b, dma, q, t, c, s
      DOUBLEPRECISION atemp, vtemp, dtemp
      INTEGER i, j, k, l
!
      DO 10000 j = 1, n
        DO 10010 i = 1, n
          v (i, j) = 0.0D0
10010     CONTINUE
        v (j, j) = 1.0D0
        d (j) = a (j, j)
10000   CONTINUE
!
      DO 20000 l = 1, nrot
        dnorm = 0.0D0
        onorm = 0.0D0
        DO 20100 j = 1, n
          dnorm = dnorm + ABS (d (j))
          DO 20110 i = 1, j - 1
            onorm = onorm + ABS (a (i, j))
20110       CONTINUE
20100     CONTINUE
        IF (onorm / dnorm .LE. 0.0D0) GOTO 19999
        DO 21000 j = 2, n
          DO 21010 i = 1, j - 1
            b = a (i, j)
            IF (ABS (b) .GT. 0.0D0) THEN
              dma = d (j) - d (i)
              IF (ABS (dma) + ABS (b) .LE. ABS (dma)) THEN
                t = b / dma
              ELSE
                q = 0.5D0 * dma / b
                t = SIGN (1.0D0 / (ABS (q) + SQRT (1.0D0 + q * q)), q)
              ENDIF
              c = 1.0D0 / SQRT (t * t + 1.0D0)
              s = t * c
              a (i, j) = 0.0D0
              DO 21110 k = 1, i - 1
                atemp    = c * a (k, i) - s * a (k, j)
                a (k, j) = s * a (k, i) + c * a (k, j)
                a (k, i) = atemp
21110           CONTINUE
              DO 21120 k = i + 1, j - 1
                atemp    = c * a (i, k) - s * a (k, j)
                a (k, j) = s * a (i, k) + c * a (k, j)
                a (i, k) = atemp
21120           CONTINUE
              DO 21130 k = j + 1, n
                atemp    = c * a (i, k) - s * a (j, k)
                a (j, k) = s * a (i, k) + c * a (j, k)
                a (i, k) = atemp
21130           CONTINUE
              DO 21140 k = 1, n
                vtemp    = c * v (k, i) - s * v (k, j)
                v (k, j) = s * v (k, i) + c * v (k, j)
                v (k, i) = vtemp
21140           CONTINUE
              dtemp = c * c * d (i) + s * s * d (j) -&
                      2.0D0 * c * s * b
              d (j) = s * s * d (i) + c * c * d (j) +&
                      2.0D0 * c * s * b
              d (i) = dtemp
              ENDIF
21010       CONTINUE
21000     CONTINUE
20000   CONTINUE
19999 CONTINUE
      nrot = l
!
      DO 30000 j = 1, n - 1
        k = j
        dtemp = d (k)
        DO 30100 i = j + 1, n
          IF (d (i) .LT. dtemp) THEN
            k = i
            dtemp = d (k)
            ENDIF
30100     CONTINUE
        IF (k .GT. j) THEN
          d (k) = d (j)
          d (j) = dtemp
          DO 30200 i = 1, n
            dtemp    = v (i, k)
            v (i, k) = v (i, j)
            v (i, j) = dtemp
30200       CONTINUE
          ENDIF
30000   CONTINUE
!
RETURN
END subroutine jacobi        

 
