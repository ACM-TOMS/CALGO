//----------------------------------------------------------------------------------
// File:        Core/Exceptions.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
//----------------------------------------------------------------------------------

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <iostream>
#include <string>

namespace cagd
{
    class Exception
    {
        //(*@\Green{// overloaded friend output to stream operator}@*)
        friend std::ostream& operator <<(std::ostream &lhs, const Exception &rhs);

    protected:
        std::string _reason;

    public:
        //(*@\Green{// special constructor}@*)
        Exception(const std::string &reason): _reason(reason)
        {
        }

        const std::string& reason() const
        {
            return _reason;
        }
    };

    inline std::ostream& operator <<(std::ostream& lhs, const Exception &rhs)
    {
        return lhs << rhs._reason;
    }
}

#endif //(*@\Green{// EXCEPTIONS\_H}@*)
