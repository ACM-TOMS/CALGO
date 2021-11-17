#pragma once

#include <EC/ECSpaces.h>

namespace cagd
{
    (*@\Green{// Represents the EC space $\mathbb{P}_n^{\alpha, \beta} = \left\langle \left\{\varphi_{n,k}\left(u\right) = u^k : u \in \left[\alpha,\beta\right] \right\}_{k=0}^n \right\rangle$ of polynomials of degree at most $n$,}@*)
    (*@\Green{// where $-\infty < \alpha < \beta < \infty$. It was shown in \citep{Carnicer1993} that its normalized B-basis is the system}@*)
    (*@\Green{// $\bigg\{b_{n,i}\left(u\right) = \binom{n}{i}\left(\frac{u-\alpha}{\beta - \alpha}\right)^i \left(\frac{\beta - u}{\beta - \alpha}\right)^{n-i} : u \in \left[\alpha,\beta\right]\bigg\}_{i=0}^n$ of Bernstein-polynomials of degree $n$.}@*)
    (*@\Green{// Note that $\mathbb{P}_n^{\alpha, \beta}$ can be identified with the solution space of the constant-coefficient homogeneous linear}@*)
    (*@\Green{// differential equation determined by the characteristic polynomial $p_{n+1}\left(z\right) = z^{n+1},~z\in\mathbb{C}.$}@*)
    class PolynomialECSpace: public ECSpace
    {
    protected:
        int _degree; (*@\Green{// $n$}@*)

    public:
        (*@\Green{// special constructor}@*)
        PolynomialECSpace(double alpha, double beta, int degree,
                          bool check_for_ill_conditioned_matrices = false,
                          int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        PolynomialECSpace* clone() const;
    };

    (*@\Green{// Represents the EC space $\mathbb{T}_{2n}^{\alpha, \beta} = \left\langle \left\{\varphi_{2n,0}\equiv 1, \left\{\varphi_{2n,2k-1}\left(u\right)=\cos\left(ku\right), \varphi_{2n,2k}\left(u\right) = \sin\left(ku\right) \right\}_{k=1}^n  : u \in \left[\alpha,\beta\right] \right\} \right\rangle$}@*)
    (*@\Green{// of trigonometric polynomials of order at most $n$ (i.e., of degree at most $2n$), where $0 < \beta-\alpha < \ell^{\prime}\left(\mathbb{T}_{2n}^{\alpha,\beta}\right)=\pi$,}@*)
    (*@\Green{// $\forall n \geq 1$. Based on \citep{Sanchez1998}, it can be shown that its normalized B-basis is the system}@*)
    (*@\Green{// $\left\{b_{n,i}\left(u\right) = c_{2n,i} \sin^{2n-i}\left(\frac{\beta - u}{2}\right) \sin^i \left(\frac{u - \alpha}{2}\right) : u \in \left[\alpha,\beta\right]\right\}_{i=0}^{2n}$, where the normalizing coefficients $\left\{c_{2n,i}\right\}_{i=0}^{2n}$}@*)
    (*@\Green{// fulfill the symmetry $c_{2n,i} = c_{2n,2n-i} = \frac{1}{\sin^{2n}\left(\frac{\beta - \alpha}{2}\right)}\sum_{r = 0}^{\left\lfloor \frac{i}{2} \right\rfloor} \binom{n}{i-r} \binom{i-r}{r} \left(2\cos\left(\frac{\beta - \alpha}{2}\right)\right)^{i-2r},~i=0,1,\ldots,n$.}@*)
    (*@\Green{// Observe that $\mathbb{T}_{2n}^{\alpha, \beta}$ coincides with the solution space of the constant-coefficient homogeneous linear differential}@*)
    (*@\Green{// equation determined by the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^{n}\left(z^2+k^2\right),~z\in\mathbb{C}.$}@*)
    class TrigonometricECSpace: public ECSpace
    {
    protected:
        int _order; (*@\Green{// $n$}@*)

    public:
        (*@\Green{// special constructor}@*)
        TrigonometricECSpace(double alpha, double beta, int order,
                             bool check_for_ill_conditioned_matrices = false,
                             int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        TrigonometricECSpace* clone() const;
    };

    (*@\Green{// Represents the EC space $\mathbb{H}_{2n}^{\alpha, \beta} = \left\langle \left\{\varphi_{2n,0}\equiv 1, \left\{\varphi_{2n,2k-1}\left(u\right)=e^{ku}, \varphi_{2n,2k}\left(u\right) = e^{-ku} \right\}_{k=1}^n  : u \in \left[\alpha,\beta\right] \right\} \right\rangle$ of}@*)
    (*@\Green{// hyperbolic polynomials of order at most $n$ (i.e., of degree at most $2n$), where $0<\beta-\alpha<\infty$.}@*)
    (*@\Green{// Based on \citep{ShenWang2005}, it can be shown that its normalized B-basis is the system}@*)
    (*@\Green{// $\left\{b_{n,i}\left(u\right) = c_{2n,i} \sinh^{2n-i}\left(\frac{\beta - u}{2}\right) \sinh^i \left(\frac{u - \alpha}{2}\right) : u \in \left[\alpha,\beta\right]\right\}_{i=0}^{2n}$, where the normalizing coefficients $\left\{c_{2n,i}\right\}_{i=0}^{2n}$}@*)
    (*@\Green{// fulfill the symmetry $c_{2n,i} = c_{2n,2n-i} = \frac{1}{\sinh^{2n}\left(\frac{\beta - \alpha}{2}\right)}\sum_{r = 0}^{\left\lfloor \frac{i}{2} \right\rfloor} \binom{n}{i-r} \binom{i-r}{r} \left(2\cosh\left(\frac{\beta - \alpha}{2}\right)\right)^{i-2r},~i=0,1,\ldots,n$.}@*)
    (*@\Green{// (In \citep{ShenWang2005}, these unique normalizing coefficients were expressed in a different way.)}@*)
    (*@\Green{// Note that $\mathbb{H}_{2n}^{\alpha, \beta}$ is in fact the solution space of the constant-coefficient homogeneous linear differential}@*)
    (*@\Green{// equation determined by the characteristic polynomial $p_{2n+1}\left(z\right) = z\prod_{k=1}^{n}\left(z^2-k^2\right),~z\in\mathbb{C}.$}@*)
    class HyperbolicECSpace: public ECSpace
    {
    protected:
        int _order; (*@\Green{// $n$}@*)

    public:
        (*@\Green{// special constructor}@*)
        HyperbolicECSpace(double alpha, double beta, int order,
                          bool check_for_ill_conditioned_matrices = false,
                          int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        HyperbolicECSpace* clone() const;
    };

    (*@\Green{// An algebraic-trigonometric EC space $\mathbb{AT}_{n\left(n+2\right)}^{\alpha,\beta}=\langle\{1,u,\ldots,u^n,\left\{u^r \cos\left(ku\right),u^r \sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n}:u \in \left[\alpha,\beta\right]\}\rangle$}@*)
    (*@\Green{// of order $n$ and dimension $\left(n+1\right)^2$. Its critical length $\ell^{\prime}\left(\mathbb{AT}_{n\left(n+2\right)}^{\alpha,\beta}\right)$ for design and the closed form of its}@*)
    (*@\Green{// normalized B-basis, in general, were not investigated in the literature. However, some special cycloidal}@*)
    (*@\Green{// subspaces of this space were studied e.g.\ in \citep{CarnicerMainarPena2004, CarnicerMainarPena2014} and references therein.}@*)
    (*@\Green{// Note that $\mathbb{AT}_{n\left(n+2\right)}^{\alpha, \beta}$ can be identified with the solution space of the constant-coefficient homogeneous linear}@*)
    (*@\Green{// differential equation defined by the characteristic polynomial $p_{\left(n+1\right)  ^{2}}\left(  z\right)=z^{n+1}\prod_{k=1}^{n}\left(  z^{2}+k^{2}\right)^{n+1-k},~z\in\mathbb{C}.$}@*)
    class ATECSpace: public ECSpace
    {
    protected:
        int _order; (*@\Green{// $n$}@*)

    public:
        (*@\Green{// special constructor}@*)
        ATECSpace(double alpha, double beta, int order,
                  bool check_for_ill_conditioned_matrices = false,
                  int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        ATECSpace* clone() const;
    };

    (*@\Green{// An algebraic-exponential-trigonometric EC space}@*)
    (*@\Green{// $\mathbb{AET}_{n\left(2n+3\right)}^{\alpha,\beta}=\left\langle\left\{\left\{u^r\right\}_{r=0}^n,\left\{u^r e^{ku} \cos\left(ku\right),u^r e^{ku}\sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n},\right.\right.$}@*)
    (*@\Green{// \hspace{2.0cm}$\left.\left.\left\{u^r e^{-ku} \cos\left(ku\right),u^r e^{-ku}\sin\left(ku\right)\right\}_{r=0,\,k=1}^{n-k,\,n}:u \in \left[\alpha,\beta\right]\right\}\right\rangle$}@*)
    (*@\Green{// of order $n$ and dimension $2n^2+3n+1$. Its critical length $\ell^{\prime}\left(\mathbb{AET}_{n\left(2n+3\right)}^{\alpha,\beta}\right)$ for design and the explicit form of}@*)
    (*@\Green{// its normalized B-basis, in general, was not investigated in the literature. Note that $\mathbb{AET}_{n\left(2n+3\right)}^{\alpha, \beta}$ coincides with}@*)
    (*@\Green{// the solution space of the constant-coefficient homogeneous linear differential equation determined by the}@*)
    (*@\Green{// characteristic polynomial $p_{2n^2+3n+1}\left(z\right)=z^{n+1}\prod_{k=1}^{n}\left[\left(z^{2}-2kz+2k^{2}\right)\cdot\left(z^{2}+2kz+2k^{2}\right)\right]^{n+1-k},~z\in\mathbb{C}.$}@*)
    class AETECSpace: public ECSpace
    {
    protected:
        int _order; (*@\Green{// $n$}@*)

    public:
        (*@\Green{// special constructor}@*)
        AETECSpace(double alpha, double beta, int order,
                   bool check_for_ill_conditioned_matrices = false,
                   int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        AETECSpace* clone() const;
    };

    (*@\Green{// The following two special EC spaces will be used for the B-representation of the ordinary exponential-}@*)
    (*@\Green{// trigonometric integral surface $\def\arraystretch{1.3}\mathbf{s}\left(u, v\right) =\left[\begin{array}{c}s^0\left(u, v\right)\\s^1\left(u, v\right)\\s^2\left(u, v\right)\end{array}\right]=\left[\begin{array}{c}\left(1-e^{\omega_0 u}\right) \cos\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\\left(e^{\omega_0 u}-1\right) \sin\left(u\right) \left(\frac{5}{4}+\cos\left(v\right)\right)\\7-e^{\omega_1 u} - \sin\left(v\right) + e^{\omega_0 u} \sin\left(v\right)\end{array}\right],\def\arraystretch{1.0}$}@*)
    (*@\Green{// where $\left(u, v\right)\in\left[\frac{7\pi}{2},\frac{49\pi}{8}\right]\times\left[-\frac{\pi}{3},\frac{5\pi}{3}\right]$, $\omega_0 = \frac{1}{6\pi}$ and $\omega_1 = \frac{1}{3\pi}$.}@*)

    (*@\Green{// The critical length for design and the explicit form of the normalized B-basis of the EC space}@*)
    (*@\Green{// $\mathbb{ET}_{6}^{\alpha, \beta}=\big\langle \big\{1,\cos\left(u\right), \sin\left(u\right),\allowbreak{}  e^{\omega_0 u},\allowbreak{} e^{\omega_1 u}, \allowbreak{}e^{\omega_0 u}\allowbreak{} \cos\left(u\right), \allowbreak{}e^{\omega_0 u} \sin\left(u\right) : u\in\left[\alpha,\beta\right]\big\} \big\rangle$ were not}@*)
    (*@\Green{// studied in the literature. Observe that $\mathbb{ET}_{6}^{\alpha, \beta}$ is in fact the solution space of the constant-coefficient}@*)
    (*@\Green{// homogeneous linear differential equation determined by the characteristic polynomial}@*)
    (*@\Green{// $p_7\left(z\right)=z\left(z-\mathbf{i}\right)\left(z+\mathbf{i}\right)\left(z-\omega_0\right)\left(z-\omega_1\right)\left(z-\left(\omega_0-\mathbf{i}\right)\right)\left(z-\left(\omega_0+\mathbf{i}\right)\right),~z\in\mathbb{C}.$}@*)
    class SnailUECSpace: public ECSpace
    {
    public:
        (*@\Green{// special constructor}@*)
        SnailUECSpace(double alpha, double beta,
                      bool check_for_ill_conditioned_matrices = false,
                      int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        SnailUECSpace* clone() const;
    };

    (*@\Green{// This is only the first order special case of the pure trigonometric EC space represented by}@*)
    (*@\Green{// the class TrigonometricECSpace.}@*)
    class SnailVECSpace: public TrigonometricECSpace
    {
    public:
        (*@\Green{// special constructor}@*)
        SnailVECSpace(double alpha, double beta,
                      bool check_for_ill_conditioned_matrices = false,
                      int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        SnailVECSpace* clone() const;
    };

    (*@\Green{// The following algebraic-hyperbolic-trigonometric EC space $\mathbb{M}_{n+4,a,b}^{\alpha,\beta} = \big\langle\big\{1,u,\ldots,u^n,\cosh\left(au\right)\cos\left(bu\right),\cosh\left(au\right)\cdot$}@*)
    (*@\Green{// $\cdot\sin\left(bu\right),\sinh\left(au\right)\cos\left(bu\right),\sinh\left(au\right)\sin\left(bu\right) : u \in \left[\alpha,\beta\right]\big\}\big\rangle$ was investigated in \citep{BrilleaudMazure2012},}@*)
    (*@\Green{// where, for $n=0$, it was shown that $\ell^{\prime}\left(\mathbb{M}_{4,a,b}^{\alpha,\beta}\right)$ coincides with the only solution of the equation}@*)
    (*@\Green{// $b\tanh\left(au\right)=a \tan\left(bu\right)$, where $u \in \left(\frac{\pi}{b},\frac{3\pi}{2b}\right)$. For example, in case of parameters $n=0$, $a=1$ and $b = 0.2$ one}@*)
    (*@\Green{// obtains that $\ell^{\prime}\left(\mathbb{M}_{4,1,0.2}^{\alpha,\beta}\right) \approx 16.694941067922716$, i.e., the space $\mathbb{M}_{4,1,0.2}^{\alpha,\beta}$ possesses a unique normalized B-basis}@*)
    (*@\Green{// provided that $\beta - \alpha \in \left(0,16.694941067922716\right)$. Moreover, $\ell^{\prime}\left(\mathbb{M}_{n+4,a,b}^{\alpha,\beta}\right) \geq \ell^{\prime}\left(\mathbb{M}_{4,a,b}^{\alpha,\beta}\right),~\forall n\geq 1$, $\forall a,b>0$.}@*)
    (*@\Green{// Note that $\mathbb{M}_{n+4,a,b}^{\alpha,\beta}$ coincides with the solution space of a constant-coefficient homogeneous linear differential}@*)
    (*@\Green{// equation determined by the characteristic polynomial}@*)
    (*@\Green{// $p_{n+5}\left(z\right) = z^{n+1}\left(z-\left(a-\mathbf{i}b\right)\right)\left(z-\left(a+\mathbf{i}b\right)\right)\left(z-\left(-a-\mathbf{i}b\right)\right)\left(z-\left(-a+\mathbf{i}b\right)\right),~z\in\mathbb{C}$.}@*)
    class BrilleaudMazureSpace: public ECSpace
    {
    protected:
        int     _order; (*@\Green{// $n \in \mathbb{N}$}@*)
        double  _a, _b; (*@\Green{// $a, b > 0$}@*)

    public:
        (*@\Green{// special constructor}@*)
        BrilleaudMazureSpace(double alpha, double beta, int order, double a, double b,
                             bool check_for_ill_conditioned_matrices = false,
                             int expected_correct_significant_digits = 5);

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        BrilleaudMazureSpace* clone() const;
    };
}
