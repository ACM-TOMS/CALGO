#include "SpecializedECSpaces.h"

#include <Core/Exceptions.h>
#include <Core/Math/Constants.h>

// References
//
// [1] M. Brilleaud and M.-L. Mazure. 2012. Mixed hyperbolic/trigonometric spaces for design.
//     Computers and Mathematics with Applications 64, 8 (2012), 2459--2477.
//     DOI: http://dx.doi.org/10.1016/j.camwa.2012.05.019
//
// [2] J.-M. Carnicer and J. M. Pena. 1993. Shape preserving representations and optimality of the
//     Bernstein basis. Advances in Computational Mathematics 1, 2 (1993), 173--196.
//     DOI: http://dx.doi.org/10.1007/BF02071384
//
// [3] J.-M. Carnicer, E. Mainar, and J. M. Pena. 2004. Critical length for design purposes and extended
//     Chebyshev spaces. Constructive Approximation 20, 1 (2004), 55--71.
//     DOI: http://dx.doi.org/10.1007/s00365-002-0530-1
//
// [4] J.-M. Carnicer, E. Mainar, and J. M. Pena. 2014. On the critical length of cycloidal spaces.
//     Constructive Approximation 39, 3 (2014), 573--583.
//     DOI: http://dx.doi.org/10.1007/s00365-013-9223-1
//
// [5] J. Sanchez-Reyes. 1998. Harmonic rational Bezier curves, p-Bezier curves and trigonometric polynomials.
//     Computer Aided Geometric Design 15, 9 (1998), 909--923.
//     DOI: http://dx.doi.org/10.1016/S0167-8396(98)00031-4
//
// [6] W.-Q. Shen and G.-Z. Wang. 2005. A class of quasi Bezier curves based on hyperbolic polynomials.
//     Journal of Zhejiang University SCIENCE 6A (Suppl. I), 9 (2005), 116--123.
//     DOI: http://dx.doi.org/10.1007/BF02887226

namespace cagd
{
    // Constructs the EC space $\mathbb{P}_n^{\alpha, \beta} = \left\langle \left\{\varphi_{n,k}\left(u\right) = u^k : u \in \left[\alpha,\beta\right] \right\}_{k=0}^n \right\rangle$ of polynomials of degree at most $n$,
    // where $-\infty < \alpha < \beta < \infty$. It was shown in \citep{Carnicer1993} that its normalized B-basis is the system
    // $\bigg\{b_{n,i}\left(u\right) = \binom{n}{i}\left(\frac{u-\alpha}{\beta - \alpha}\right)^i \left(\frac{\beta - u}{\beta - \alpha}\right)^{n-i} : u \in \left[\alpha,\beta\right]\bigg\}_{i=0}^n$ of Bernstein-polynomials of degree $n$.
    // Note that $\mathbb{P}_n^{\alpha, \beta}$ can be identified with the solution space of the constant-coefficient homogeneous linear
    // differential equation determined by the characteristic polynomial $p_{n+1}\left(z\right) = z^{n+1},~z\in\mathbb{C}.$
    PolynomialECSpace::PolynomialECSpace(
        double alpha, double beta, int degree,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _degree(degree)
    {
        if (_degree < 0)
        {
            throw Exception("The degree of the polynomial EC space should be a "
                            "non-negative integer!");
        }

        // we store the single zero $z=0$ of order $n+1$ of the characteristic polynomial $p_{n+1}\left(z\right) = z^{n+1},~z\in\mathbb{C}$
        if (!insertZero(0.0, 0.0, _degree + 1, true,
                        check_for_ill_conditioned_matrices,
                        expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "vector space of polynomials of finite degree could not "
                            "be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    PolynomialECSpace* PolynomialECSpace::clone() const
    {
        return new (std::nothrow) PolynomialECSpace(*this);
    }

    // Builds the EC space $\mathbb{T}_{2n}^{\alpha, \beta} = \left\langle \left\{\varphi_{2n,0}\equiv 1, \left\{\varphi_{2n,2k-1}\left(u\right)=\cos\left(ku\right), \varphi_{2n,2k}\left(u\right) = \sin\left(ku\right) \right\}_{k=1}^n  : u \in \left[\alpha,\beta\right] \right\} \right\rangle$
    // of trigonometric polynomials of order at most $n$ (i.e., of degree at most $2n$), where $0 < \beta-\alpha < \ell^{\prime}\left(\mathbb{T}_{2n}^{\alpha,\beta}\right)=\pi$,
    // $\forall n \geq 1$. Based on \citep{Sanchez1998}, it can be shown that its normalized B-basis is the system
    // $\left\{b_{n,i}\left(u\right) = c_{2n,i} \sin^{2n-i}\left(\frac{\beta - u}{2}\right) \sin^i \left(\frac{u - \alpha}{2}\right) : u \in \left[\alpha,\beta\right]\right\}_{i=0}^{2n}$, where the normalizing coefficients $\left\{c_{2n,i}\right\}_{i=0}^{2n}$
    // fulfill the symmetry $c_{2n,i} = c_{2n,2n-i} = \frac{1}{\sin^{2n}\left(\frac{\beta - \alpha}{2}\right)}\sum_{r = 0}^{\left\lfloor \frac{i}{2} \right\rfloor} \binom{n}{i-r} \binom{i-r}{r} \left(2\cos\left(\frac{\beta - \alpha}{2}\right)\right)^{i-2r},~i=0,1,\ldots,n$.
    // Observe that $\mathbb{T}_{2n}^{\alpha, \beta}$ coincides with the solution space of the constant-coefficient homogeneous linear differential
    // equation determined by the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^{n}\left(z^2+k^2\right),~z\in\mathbb{C}.$
    TrigonometricECSpace::TrigonometricECSpace(
        double alpha, double beta, int order,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _order(order)
    {
        if (_order < 0)
        {
            throw Exception("The order of the vector space of pure trigonometric "
                            "polynomials should be a non-negative integer!");
        }

        if (beta - alpha >= PI)
        {
            throw Exception("The vector space of trigonometric polynomials of finite "
                            "order is EC provided that the length of its definition "
                            "domain is strictly less than pi!");
        }

        // we store the first order zeros $z=\pm\mathbf{i}k$ of the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^n\left(z^2+k^2\right),~z\in\mathbb{C}$
        // note that, in case of complex roots it is sufficient to store one of the conjugately equivalent zeros and the
        // first order real zero $z=0$ was already stored by the constructor of the base class ECSpace
        for (int k = 1; k <= _order; k++)
        {
            insertZero(0.0, k, 1, false);
        }

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "vector space of trigonometric polynomials of finite order "
                            "could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    TrigonometricECSpace* TrigonometricECSpace::clone() const
    {
        return new (std::nothrow) TrigonometricECSpace(*this);
    }

    // Constructs the EC space $\mathbb{H}_{2n}^{\alpha, \beta} = \left\langle \left\{\varphi_{2n,0}\equiv 1, \left\{\varphi_{2n,2k-1}\left(u\right)=e^{ku}, \varphi_{2n,2k}\left(u\right) = e^{-ku} \right\}_{k=1}^n  : u \in \left[\alpha,\beta\right] \right\} \right\rangle$ of
    // hyperbolic polynomials of order at most $n$ (i.e., of degree at most $2n$), where $0<\beta-\alpha<\infty$.
    // Based on \citep{ShenWang2005}, it can be shown that its normalized B-basis is the system
    // $\left\{b_{n,i}\left(u\right) = c_{2n,i} \sinh^{2n-i}\left(\frac{\beta - u}{2}\right) \sinh^i \left(\frac{u - \alpha}{2}\right) : u \in \left[\alpha,\beta\right]\right\}_{i=0}^{2n}$, where the normalizing coefficients $\left\{c_{2n,i}\right\}_{i=0}^{2n}$
    // fulfill the symmetry $c_{2n,i} = c_{2n,2n-i} = \frac{1}{\sinh^{2n}\left(\frac{\beta - \alpha}{2}\right)}\sum_{r = 0}^{\left\lfloor \frac{i}{2} \right\rfloor} \binom{n}{i-r} \binom{i-r}{r} \left(2\cosh\left(\frac{\beta - \alpha}{2}\right)\right)^{i-2r},~i=0,1,\ldots,n$.
    // (In \citep{ShenWang2005}, these unique normalizing coefficients were expressed in a different way.)
    // Note that $\mathbb{H}_{2n}^{\alpha, \beta}$ is in fact the solution space of the constant-coefficient homogeneous linear differential
    // equation determined by the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^{n}\left(z^2-k^2\right),~z\in\mathbb{C}.$
    HyperbolicECSpace::HyperbolicECSpace(
        double alpha, double beta, int order,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _order(order)
    {
        if (_order < 0)
        {
            throw Exception("The order of the vector space of pure hyperbolic "
                            "polynomials should be a non-negative integer!");
        }

        // we store the first order real zeros $z=\pm k$ of the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^n\left(z^2-k^2\right)$
        // note that, first order real zero $z=0$ was already stored by the constructor of the base class ECSpace
        for (int k = 1; k <= _order; k++)
        {
            insertZero( k, 0.0, 1, false);
            insertZero(-k, 0.0, 1, false);
        }

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "vector space of hyperbolic polynomials of finite order "
                            "could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    HyperbolicECSpace* HyperbolicECSpace::clone() const
    {
        return new (std::nothrow) HyperbolicECSpace(*this);
    }

    // Constructs the algebraic-trigonometric EC space
    // $\mathbb{AT}_{n\left(n+2\right)}^{\alpha,\beta}=\langle\{1,u,\ldots,u^n,\left\{u^r \cos\left(ku\right),u^r \sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n}:u \in \left[\alpha,\beta\right]\}\rangle$
    // of order $n$ and dimension $\left(n+1\right)^2$. Its critical length $\ell^{\prime}\left(\mathbb{AT}_{n\left(n+2\right)}^{\alpha,\beta}\right)$ for design and the closed form of its
    // normalized B-basis, in general, were not investigated in the literature. However, some special cycloidal
    // subspaces of this space were studied e.g.\ in \citep{CarnicerMainarPena2004, CarnicerMainarPena2014} and references therein.
    // Note that $\mathbb{AT}_{n\left(n+2\right)}^{\alpha, \beta}$ can be identified with the solution space of the constant-coefficient homogeneous linear
    // differential equation defined by the characteristic polynomial $p_{\left(n+1\right)  ^{2}}\left(  z\right)=z^{n+1}\prod_{k=1}^{n}\left(  z^{2}+k^{2}\right)^{n+1-k},~z\in\mathbb{C}.$
    ATECSpace::ATECSpace(
        double alpha, double beta, int order,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _order(order)
    {
        if (_order < 0)
        {
            throw Exception("The order of the mixed algebraic-trigonometric vector "
                            "space should be a non-negative integer!");
        }

        // we store the (higher order) zeros of the characteristic polynomial $p_{\left(n+1\right)^2}\left(z\right) = z^{n+1}\prod_{k=1}^n\left(z^2+k^2\right)^{n+1-k}$
        insertZero(0.0, 0.0, _order + 1, false);

        // note that, in case of complex roots it is sufficient to store one of the conjugately equivalent zeros
        for (int k = 1; k <= _order; k++)
        {
            insertZero(0.0, (double)k, _order + 1 - k, false);
        }

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "algebraic-trigonometric vector space could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    ATECSpace* ATECSpace::clone() const
    {
        return new (std::nothrow) ATECSpace(*this);
    }

    // Constructs the algebraic-exponential-trigonometric EC space
    // $\mathbb{AET}_{n\left(2n+3\right)}^{\alpha,\beta}=\left\langle\left\{\left\{u^r\right\}_{r=0}^n,\left\{u^r e^{ku} \cos\left(ku\right),u^r e^{ku}\sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n},\right.\right.$
    // \hspace{2.0cm}$\left.\left.\left\{u^r e^{-ku} \cos\left(ku\right),u^r e^{-ku}\sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n}:u \in \left[\alpha,\beta\right]\right\}\right\rangle$
    // of order $n$ and dimension $2n^2+3n+1$. Its critical length $\ell^{\prime}\left(\mathbb{AET}_{n\left(2n+3\right)}^{\alpha,\beta}\right)$ for design and the explicit form of
    // its normalized B-basis, in general, was not investigated in the literature. Note that $\mathbb{AET}_{n\left(2n+3\right)}^{\alpha, \beta}$ coincides with
    // the solution space of the constant-coefficient homogeneous linear differential equation determined by the
    // characteristic polynomial $p_{2n^2+3n+1}\left(z\right)=z^{n+1}\prod_{k=1}^{n}\left[\left(z^{2}-2kz+2k^{2}\right)\cdot\left(z^{2}+2kz+2k^{2}\right)\right]^{n+1-k},~z\in\mathbb{C}.$
    AETECSpace::AETECSpace(
        double alpha, double beta, int order,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _order(order)
    {
        if (_order < 0)
        {
            throw Exception("The order of the mixed algebraic-exponential-trigonometric "
                            "vector space should be a non-negative integer!");
        }

        // we store the (higher order) zeros of the characteristic polynomial
        // $p_{2n^2+3n+1}\left(z\right) = z^{n+1}\prod_{k=1}^n\left[\left(z^2-2kz+2k^2\right)\cdot\left(z^2+2kz+2k^2\right)\right]^{n+1-k},~z\in\mathbb{C}$
        insertZero(0.0, 0.0, _order + 1, false);

        // note that, it is sufficient to store one of the conjugately equivalent zeros $z=k\pm\mathbf{i}k$ and $z=-k\pm\mathbf{i}k$
        // of order $n+1-k$
        for (int k = 1; k <= _order; k++)
        {
            insertZero( k, k, _order + 1 - k, false);
            insertZero(-k, k, _order + 1 - k, false);
        }

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "algebraic-exponential-trigonometric vector space "
                            "could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    AETECSpace* AETECSpace::clone() const
    {
        return new (std::nothrow) AETECSpace(*this);
    }

    // The following two special EC spaces will be used for the B-representation of the ordinary exponential-
    // trigonometric integral surface $\def\arraystretch{1.3}\mathbf{s}\left(u, v\right) =\left[\begin{array}{c}s^0\left(u, v\right)\\s^1\left(u, v\right)\\s^2\left(u, v\right)\end{array}\right]=\left[\begin{array}{c}\left(1-e^{\omega_0 u}\right) \cos\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\\left(e^{\omega_0 u}-1\right) \sin\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\7-e^{\omega_1 u} - \sin\left(v\right) + e^{\omega_0 u} \sin\left(v\right)\end{array}\right],\def\arraystretch{1.0}$
    // where $\left(u, v\right)\in\left[\frac{7\pi}{2},\frac{49\pi}{8}\right]\times\left[-\frac{\pi}{3},\frac{5\pi}{3}\right]$, $\omega_0 = \frac{1}{6\pi}$ and $\omega_1 = \frac{1}{3\pi}$.

    // The critical length for design and the explicit form of the normalized B-basis of the EC space
    // $\mathbb{ET}_{6}^{\alpha, \beta}=\big\langle \big\{1,\cos\left(u\right), \sin\left(u\right),\allowbreak{}  e^{\omega_0 u},\allowbreak{} e^{\omega_1 u}, \allowbreak{}e^{\omega_0 u}\allowbreak{} \cos\left(u\right), \allowbreak{}e^{\omega_0 u} \sin\left(u\right) : u\in\left[\alpha,\beta\right]\big\} \big\rangle$ were not
    // studied in the literature. Observe that $\mathbb{ET}_{6}^{\alpha, \beta}$ is in fact the solution space of the constant-coefficient
    // homogeneous linear differential equation determined by the characteristic polynomial
    // $p_7\left(z\right)=z\left(z-\mathbf{i}\right)\left(z+\mathbf{i}\right)\left(z-\omega_0\right)\left(z-\omega_1\right)\left(z-\left(\omega_0-\mathbf{i}\right)\right)\left(z-\left(\omega_0+\mathbf{i}\right)\right),~z\in\mathbb{C}.$

    // special constructor
    SnailUECSpace::SnailUECSpace(
        double alpha, double beta,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits)
    {
        // we store one of the conjugately equivalent first order zeros of the characteristic polynomial
        // $p_7\left(z\right)=z\left(z-\mathbf{i}\right)\left(z+\mathbf{i}\right)\left(z-\omega_0\right)\left(z-\omega_1\right)\left(z-\left(\omega_0-\mathbf{i}\right)\right)\left(z-\left(\omega_0+\mathbf{i}\right)\right),~z\in\mathbb{C}$
        insertZero(0.0,            1.0, 1, false);
        insertZero(1.0 / 6.0 / PI, 0.0, 1, false);
        insertZero(1.0 / 3.0 / PI, 0.0, 1, false);
        insertZero(1.0 / 6.0 / PI, 1.0, 1, false);

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "SnailUECSpace could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    SnailUECSpace* SnailUECSpace::clone() const
    {
        return new (std::nothrow) SnailUECSpace(*this);
    }

    // Constructs the first order special case of the pure trigonometric EC space represented by
    // the class TrigonometricECSpace.
    SnailVECSpace::SnailVECSpace(
        double alpha, double beta,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        TrigonometricECSpace(alpha, beta, 1,
                             check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits)
    {
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    SnailVECSpace* SnailVECSpace::clone() const
    {
        return new (std::nothrow) SnailVECSpace(*this);
    }

    // The following algebraic-hyperbolic-trigonometric EC space $\mathbb{M}_{n+4,a,b}^{\alpha,\beta} = \big\langle\big\{1,u,\ldots,u^n,\cosh\left(au\right)\cos\left(bu\right),\cosh\left(au\right)\cdot$}
    // $\cdot\sin\left(bu\right),\sinh\left(au\right)\cos\left(bu\right),\sinh\left(au\right)\sin\left(bu\right) : u \in \left[\alpha,\beta\right]\big\}\big\rangle$ was investigated in \citep{BrilleaudMazure2012},}
    // where, for $n=0$, it was shown that $\ell^{\prime}\left(\mathbb{M}_{4,a,b}^{\alpha,\beta}\right)$ coincides with the only solution of the equation}
    // $b\tanh\left(au\right)=a \tan\left(bu\right)$, where $u \in \left(\frac{\pi}{b},\frac{3\pi}{2b}\right)$. For example, in case of parameters $n=0$, $a=1$ and $b = 0.2$ one}
    // obtains that $\ell^{\prime}\left(\mathbb{M}_{4,1,0.2}^{\alpha,\beta}\right) \approx 16.694941067922716$, i.e., the space $\mathbb{M}_{4,1,0.2}^{\alpha,\beta}$ possesses a unique normalized B-basis}
    // provided that $\beta - \alpha \in \left(0,16.694941067922716\right)$. Moreover, $\ell^{\prime}\left(\mathbb{M}_{n+4,a,b}^{\alpha,\beta}\right) \geq \ell^{\prime}\left(\mathbb{M}_{4,a,b}^{\alpha,\beta}\right),~\forall n\geq 1$, $\forall a,b>0$.}
    // Note that $\mathbb{M}_{n+4,a,b}^{\alpha,\beta}$ coincides with the solution space of a constant-coefficient homogeneous linear differential}
    // equation determined by the characteristic polynomial}
    // $p_{n+5}\left(z\right) = z^{n+1}\left(z-\left(a-\mathbf{i}b\right)\right)\left(z-\left(a+\mathbf{i}b\right)\right)\left(z-\left(-a-\mathbf{i}b\right)\right)\left(z-\left(-a+\mathbf{i}b\right)\right),~z\in\mathbb{C}$.}

    // special constructor
    BrilleaudMazureSpace::BrilleaudMazureSpace(
        double alpha, double beta, int order, double a, double b,
        bool check_for_ill_conditioned_matrices,
        int expected_correct_significant_digits):
        ECSpace(alpha, beta,
                check_for_ill_conditioned_matrices, expected_correct_significant_digits),
        _order(order),
        _a(a), _b(b)
    {
        if (_order < 0)
        {
            throw Exception("The order of the mixed algebraic-hyperbolic-trigonometric "
                            "vector space should be a non-negative integer!");
        }

        if (_a <= 0.0 || _b <= 0.0)
        {
            throw Exception("Parameters a and b should be positive reals!");
        }

        // we store the (first and higher order) zeros of the characteristic polynomial
        // $p_{n+5}\left(z\right) = z^{n+1}\left(z-\left(a-\mathbf{i}b\right)\right)\left(z-\left(a+\mathbf{i}b\right)\right)\left(z-\left(-a-\mathbf{i}b\right)\right)\left(z-\left(-a+\mathbf{i}b\right)\right),~a,b>0,~z\in\mathbb{C}$
        insertZero(0.0, 0.0, _order + 1, false);

        insertZero( a, b, 1, false);
        insertZero(-a, b, 1, false);

        if (!updateBothBases(check_for_ill_conditioned_matrices,
                             expected_correct_significant_digits))
        {
            throw Exception("The ordinary basis and the normalized B-basis of the "
                            "BrilleaudMazureSpace could not be updated!");
        }
    }

    // clone function required by smart pointers based on the deep copy ownership policy
    BrilleaudMazureSpace* BrilleaudMazureSpace::clone() const
    {
        return new (std::nothrow) BrilleaudMazureSpace(*this);
    }
}
