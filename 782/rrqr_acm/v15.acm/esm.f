c     *---------------------------------------------------------------*
c     i                                                               i
c     i                   entrada/salida matricial                    i
c     i                                                               i
c     *---------------------------------------------------------------*


c     *****************************************************************
c     *                           escmatriz                           *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar una matriz de *
c     *  dimension m x n cuyos componentes son reales en simple       *
c     *  precision. se escribe cada elemento en una linea.            *
c     *****************************************************************
      subroutine escmatriz (a,lda,m,n)
      integer lda
      integer m,n
      real    a(lda,n)

      integer f,c
      do f=1,m
          do c=1,n
              print *,' a(',f,' , ',c,' ) : ',a(f,c)
          end do
      end do
      return
      end


c     *****************************************************************
c     *                            escmaf                             *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar una matriz de *
c     *  dimension m x n cuyos componentes son reales en simple       *
c     *  precision. todos los elementos de una fila se escriben en    *
c     *  una misma fila de pantalla, si caben.                        *
c     *****************************************************************
      subroutine escmaf (a,lda,m,n)
      integer lda
      integer m,n
      real    a(lda,n)

      integer f,c
      do f=1,m
          WRITE (*,10) (a(f,c), c=1,n)
 10       format (1x, 100(e12.5))
      end do
      return
      end


c     *****************************************************************
c     *                           escvector                           *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar un vector de  *
c     *  n numeros reales en simple precision, situando uno en cada   *
c     *  linea.                                                       *
c     *****************************************************************
      subroutine escvector (v,n)
      integer n
      real    v(n)

      integer f
      do f=1,n
          print *,'el elemento v( ',f,'  ) es: ',v(f)
      end do
      return
      end


c     *****************************************************************
c     *                           escvef                              *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar un vector de  *
c     *  n numeros reales en simple precision. en este caso se situan *
c     *  todos los elementos en una sola linea.                       *
c     *****************************************************************
      subroutine escvef (v,n)
      integer n
      real    v(n)

      integer f
      WRITE (*,10) (v(f), f=1,n)
 10   format (1x, 100(e12.5))
      return
      end


c     *****************************************************************
c     *                           escvint                             *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar un vector de  *
c     *  n numeros enteros, situando uno en cada linea.               *
c     *****************************************************************
      subroutine escvint (v,n)
      integer n
      integer v(n)

      integer f
      do f=1,n
          print *,'el elemento v( ',f,'  ) es: ',v(f)
      end do
      return
      end


c     *****************************************************************
c     *                           escvif                              *
c     *****************************************************************
c     *    esta subrutina escribe en la salida estandar un vector de  *
c     *  n numeros reales en simple precision. en este caso se situan *
c     *  todos los elementos en una sola linea.                       *
c     *****************************************************************
      subroutine escvif (v,n)
      integer n
      integer v(n)

      integer f
      WRITE (*,10) (v(f), f=1,n)
 10   format (1x, 100(i8))
      return
      end


c     *****************************************************************
c     *                           espera                              *
c     *****************************************************************
c     *    esta subrutina lee un caracter del teclado, para lo cual   *
c     *  detiene momentaneamente la ejecucion del programa.           *
c     *****************************************************************
      subroutine espera
      character*1 c

      read (*,10) c
  10  format (a1)
      return
      end


c     *****************************************************************
c     *                             leed                              *
c     *****************************************************************
c     *    esta subrutina lee un entero de la entrada estandar        *
c     *  comprendido entre 1 y "max". si el numero tecleado no esta   *
c     *  dentro de este rango, se vuelve a demandar otro.             *
c     *****************************************************************
      subroutine leed ( max,n )
      integer max,n

      n = 0
      do while ( (n .lt. 1) .or. (n .gt. max) )
          print *,'dame un numero  1 .. ',max,'  : '
          read (*,*) n
      end do
      return
      end

