C--**--CH2470--734--P:RW--22:9:1999
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> rosenbrockfunc.f90
MODULE RosenbrockFunc

    USE Precisions, only    : stnd, short
    USE Supp_Codes, only    : justf, justg, both, OK
    USE Reals,      only    : one, two, c100, c200, c400

CONTAINS

    SUBROUTINE Rosenbrock ( x, f, g, job )

        Implicit None

        Real(stnd),     intent(IN)     :: x(:)
        Real(stnd),     intent(OUT)    :: f, g(size(x))

        Integer(short), intent(INOUT)  :: job

        Real(stnd)  :: w1, w2       ! temporaries

        Logical     :: dof, dog

        dof =    job == justf  .or. job == both
        dog =    job == justg  .or. job == both

        w1 = x(2) - x(1)**2;    w2 =  one - x(1)

        IF ( dof ) then
             f = c100*w1**2 + w2**2
        end if

        IF ( dog ) then
             g(1) = -c400 * w1 * x(1) - two*w2
             g(2) =  c200 * w1
        end if

        job =  OK
    return
    end Subroutine Rosenbrock

end Module RosenbrockFunc
! <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> <<=>> example.f90
Program Min_Rosenbrock

    USE Precisions, only    : stnd, short
    USE Minimize
    USE RosenbrockFunc
    USE Reals,          only :  one

    Implicit None

    Integer(short),PARAMETER :: n = 2
    Integer(short)           :: state
    Integer                  :: nout, i1mach

    Real(stnd),   PARAMETER  :: acc = 0.001
    Real(stnd)               :: f, x(n), g(n)

    Type(MinimizeState)      :: C

    x = (/ -1.2_stnd, one /)  ! set the initial guess.

    state = Normal
    nout = i1mach(2)

    Call Minimize_f (Rosenbrock, x, f, g, acc, state, 0, C, &
                     Frequency = 10, EvalLimit = 200 )

    IF ( state == Done ) then
       write(nout,'(1x,a,e16.8,a)') &
           'Least function value:      ',f, '.'
       write(nout,'(1x,a,i5,a)') &
           'Function evaluation count: ',C%EvalCts%FEvals,'.'
       write(nout,'(1x,a,i5,a)') &
           'Iteration count:           ',C%Shared%it,'.'
    else
       write(nout,'(1x,a,i5)') &
           'Failure: state = ', state
    end if
    stop
end
