//----------------------------------------------------------------------------------
// File:        EC/CharacteristicPolynomials.cpp
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#include "CharacteristicPolynomials.h"
#include "../Core/Math/Constants.h"
#include <cmath>
#include <algorithm>

using namespace std;

namespace cagd
{
    //(*@\Green{// special/default constructor}@*)
    CharacteristicPolynomial::Zero::Zero(double real, double imaginary, int order):
            real(real), imaginary(imaginary), order(order >= 1 ? order : 1)
    {
        assert("The multiplicity (or the order) of the given complex zero should be "
               "at least 1!" && order >= 1);
    }

    //(*@\Green{// binary logical operator for possible comparisons (conjugate roots are considered to be equivalent)}@*)
    bool CharacteristicPolynomial::Zero::operator ==(
            const CharacteristicPolynomial::Zero &rhs) const
    {
        return ((real == rhs.real) && (abs(imaginary) == abs(rhs.imaginary)));
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    CharacteristicPolynomial::Zero* CharacteristicPolynomial::Zero::clone() const
    {
        return new (nothrow) Zero(*this);
    }

    //(*@\Green{// default constructor}@*)
    CharacteristicPolynomial::CharacteristicPolynomial():
            _factorization(0), _factorization_changed(false)
    {
    }

    //(*@\Green{// If during insertion the given zero (or its conjugate) already exists, the value of its order will be changed}@*)
    //(*@\Green{// to the maximum of its formerly and newly given orders.}@*)
    void CharacteristicPolynomial::insertZero(double real, double imaginary, int order)
    {
        Zero z(real, imaginary, order);

        vector<Zero>::iterator i = find(_factorization.begin(), _factorization.end(), z);

        if (i == _factorization.end())
        {
            _factorization.push_back(z);
            _factorization_changed = true;
        }
        else
        {
            if (i->order < order)
            {
                i->order = order;
                _factorization_changed = true;
            }
            else
            {
                _factorization_changed = false;
            }
        }
    }

    void CharacteristicPolynomial::insertZero(const Zero &zero)
    {
        vector<Zero>::iterator i = find(_factorization.begin(), _factorization.end(),
                                        zero);

        if (i == _factorization.end())
        {
            _factorization.push_back(zero);
            _factorization_changed = true;
        }
        else
        {
            if (i->order < zero.order)
            {
                i->order = zero.order;
                _factorization_changed = true;
            }
            else
            {
                _factorization_changed = false;
            }
        }
    }

    //(*@\Green{// If exists, deletes the given zero (or its conjugate) independently of its order.}@*)
    bool CharacteristicPolynomial::deleteZero(double real, double imaginary)
    {
        Zero z(real, imaginary);

        vector<Zero>::iterator i = find(_factorization.begin(), _factorization.end(), z);

        if (i == _factorization.end())
        {
            _factorization_changed = false;

            return false;
        }

        _factorization.erase(i);

        _factorization_changed = true;

        return true;
    }

    bool CharacteristicPolynomial::deleteZero(const Zero &zero)
    {
        vector<Zero>::iterator i = find(_factorization.begin(), _factorization.end(),
                                        zero);

        if (i == _factorization.end())
        {
            _factorization_changed = false;

            return false;
        }

        _factorization.erase(i);

        _factorization_changed = true;

        return true;
    }

    //(*@\Green{// returns the state of the flag \_factorization\_changed}@*)
    bool CharacteristicPolynomial::factorizationChanged() const
    {
        return _factorization_changed;
    }

    //(*@\Green{// flips the state of the flag \_factorization\_changed}@*)
    void CharacteristicPolynomial::flipFactorizationChangedState()
    {
        _factorization_changed = !_factorization_changed;
    }

    //(*@\Green{// overloaded function operator that evaluates the characteristic polynomial (\mref{eq:characteristic_polynomial}) at the parameter value $u$}@*)
    double CharacteristicPolynomial::operator ()(double u) const
    {
        double product = 1.0;
        for (vector<Zero>::const_iterator z = _factorization.begin();
             z != _factorization.end(); z++)
        {
            if (z->imaginary)
            {
                product *= pow(z->real * z->real + z->imaginary * z->imaginary
                               - 2.0 * z->real * u + u * u, z->order);
            }
            else
            {
                product *= pow(u - z->real, z->order);
            }
        }

        return product;
    }

    //(*@\Green{// determines whether the characteristic polynomial is an either odd or even function}@*)
    bool CharacteristicPolynomial::isEvenOrOddFunction() const
    {
        if (_factorization.empty())
        {
            return false;
        }

        for (vector<Zero>::const_iterator z = _factorization.begin();
             z != _factorization.end(); z++)
        {
            if (z->imaginary == 0.0)
            {
                CharacteristicPolynomial::Zero symmetric(-z->real, 0.0);
                vector<Zero>::const_iterator s = find(_factorization.begin(),
                                                      _factorization.end(), symmetric);

                if (s == _factorization.end())
                {
                    return false;
                }
                else
                {
                    if (z->order != s->order)
                    {
                        return false;
                    }
                }

            }
        }
        return true;
    }

    //(*@\Green{// clone function required by smart pointers based on the deep copy ownership policy}@*)
    CharacteristicPolynomial* CharacteristicPolynomial::clone() const
    {
        return new (nothrow) CharacteristicPolynomial(*this);
    }

    //(*@\Green{// overloaded output to stream operators}@*)
    ostream& operator <<(ostream &lhs, const CharacteristicPolynomial::Zero &rhs)
    {
        return lhs << rhs.real << " " << rhs.imaginary << " " << rhs.order;
    }

    ostream& operator <<(ostream &lhs, const CharacteristicPolynomial &rhs)
    {
        lhs << rhs._factorization.size() << endl;

        for (vector<CharacteristicPolynomial::Zero>::const_iterator
             i = rhs._factorization.begin(); i != rhs._factorization.end(); i++)
        {
            lhs << *i << endl;
        }

        return lhs;
    }
}
