//----------------------------------------------------------------------------------
// File:        EC/CharacteristicPolynomials.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef CHARACTERISTICPOLYNOMIALS_H
#define CHARACTERISTICPOLYNOMIALS_H

#include <iostream>
#include <vector>

namespace cagd
{
    (*@\Green{// represents the characteristic polynomial (\mref{eq:characteristic_polynomial}) of the constant-coefficient linear homogeneous differential}@*)
    (*@\Green{// equation (\mref{eq:differential_equation})}@*)
    class CharacteristicPolynomial
    {
        (*@\Green{// friend class that represents an extended Chebyshev (EC) vector space of functions}@*)
        friend class ECSpace;

        (*@\Green{// friend output to stream operator}@*)
        friend std::ostream& operator <<(std::ostream &lhs,
                                         const CharacteristicPolynomial &rhs);

    public:
        (*@\Green{// the nested public class Zero represents a higher order complex root of the characteristic polynomial (\mref{eq:characteristic_polynomial})}@*)
        class Zero
        {
        public:
            double  real;      (*@\Green{// the real part of the complex root}@*)
            double  imaginary; (*@\Green{// the imaginary part of the complex root}@*)
            int     order;     (*@\Green{// the multiplicity or the order of the complex root, its value does}@*)
                               (*@\Green{// not affect the comparison operator == declared below in line \mref{src:CharacteristicPolynomial::comparison_operator}}@*)

            (*@\Green{// special/default constructor}@*)
            explicit Zero(double real = 0.0, double imaginary = 0.0, int order = 1);

            (*@\Green{// binary logical operator for possible comparisons (conjugate roots are considered to be equivalent)}@*)
            bool operator ==(const Zero &rhs) const;(*@\label{src:CharacteristicPolynomial::comparison_operator}@*)

            (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
            Zero* clone() const;
        };

    protected:
        (*@\Green{// the factorization of the characteristic polynomial (\mref{eq:characteristic_polynomial}) (we only store conjugately pairwise different zeros)}@*)
        std::vector<Zero> _factorization;

        (*@\Green{// a flag that will be set to true whenever the factorization of the characteristic polynomial (\mref{eq:characteristic_polynomial}) or the}@*)
        (*@\Green{// definition domain of the underlying EC space is changed}@*)
        bool              _factorization_changed;

    public:
        (*@\Green{// default constructor}@*)
        CharacteristicPolynomial();

        (*@\Green{// If during insertion the given zero (or its conjugate) already exists, the value of its order will be changed}@*)
        (*@\Green{// to the maximum of its formerly and newly given orders.}@*)
        void insertZero(double real, double imaginary, int order);
        void insertZero(const Zero &zero);

        (*@\Green{// If exists, deletes the given zero (or its conjugate) independently of its order.}@*)
        bool deleteZero(double real, double imaginary);
        bool deleteZero(const Zero &zero);

        bool factorizationChanged() const;    (*@\Green{// returns the state of the flag \_factorization\_changed}@*)
        void flipFactorizationChangedState(); (*@\Green{// flips the state of the flag \_factorization\_changed}@*)

        (*@\Green{// overloaded function operator that evaluates the characteristic polynomial (\mref{eq:characteristic_polynomial}) at the parameter value $u$}@*)
        double operator ()(double u) const;

        (*@\Green{// determines whether the characteristic polynomial is an either odd or even function}@*)
        bool isEvenOrOddFunction() const;

        (*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
        CharacteristicPolynomial* clone() const;
    };

    (*@\Green{// overloaded output to stream operators}@*)
    std::ostream& operator <<(std::ostream &lhs,
                              const CharacteristicPolynomial::Zero &rhs);
    std::ostream& operator <<(std::ostream &lhs, const CharacteristicPolynomial &rhs);
}

#endif (*@\Green{// CHARACTERISTICPOLYNOMIALS\_H}@*)
