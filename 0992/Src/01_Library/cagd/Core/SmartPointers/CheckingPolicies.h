//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/CheckingPolicies.h
// Library:     cagd{32|64}[d]
//              An OpenGL and C++ based function library for curve and surface
//              modeling in a large class of extended Chebyshev spaces
// Version:     1.0.0 | August 15, 2017 - October 11, 2018
// Author:      (c) Agoston Roth
// Affiliation: Babes-Bolyai University
//              Department of Mathematics and Computer Science of the Hungarian Line
//              RO-400084, Cluj-Napoca, Romania
// Reference:   A. Alexandrescu. 2001. Modern C++ Design: Generic Programming and Design
//              Patterns Applied, 1st edition. Addison-Wesley Professional, USA.
//--------------------------------------------------------------------------------------

#ifndef CHECKINGPOLICIES_H
#define CHECKINGPOLICIES_H

#include <cassert>
#include <stdexcept>

namespace cagd
{
    //(*@\Green{// The helper class NullPointerException is used by the RejectNullDereferenceOrIndirection and RejectNull}@*)
    //(*@\Green{// checking policies. (Each instantiation attempt will generate a run-time error.)}@*)
    class NullPointerException: public std::runtime_error
    {
    public:
        NullPointerException() throw(): std::runtime_error("Null pointer exception!") {}
    };

    //(*@\Green{// The user can choose from the:}@*)
    //(*@\Green{// \hspace{0.5cm}- RejectNullDereferenceOrIndirection;}@*)
    //(*@\Green{// \hspace{0.5cm}- RejectNull;}@*)
    //(*@\Green{// \hspace{0.5cm}- AssertNullDereferenceOrIndirection; and}@*)
    //(*@\Green{// \hspace{0.5cm}- AssertNull}@*)
    //(*@\Green{// checking policies.}@*)
    template <typename T>
    class CheckingPolicy
    {
    public:
        //(*@\Green{// The NoCheck policy will not perform checkings during initialization, dereference and indirection.}@*)
        class NoCheck
        {
        public:
            static void onInitialize(const T* ptr) {}
            static void onDereferenceOrIndirection(const T* ptr) {}
        };

        //(*@\Green{// The RejectNullDereferenceOrIndirection checking policy will throw a NullPointerException whenever}@*)
        //(*@\Green{// one attempts either to dereference or indirect a smart pointer that stores a null raw pointer.}@*)
        class RejectNullDereferenceOrIndirection
        {
        public:
            static void onInitialize(const T* /*ptr*/) {}

            static void onDereferenceOrIndirection(const T* ptr)
            {
                if (!ptr)
                {
                    throw NullPointerException();
                }
            }
        };

        //(*@\Green{// The RejectNull checking policy will throw a NullPointerException whenever one bumps into null raw}@*)
        //(*@\Green{// pointers during the initialization/assignment, dereference and indirection of smart pointers.}@*)
        class RejectNull
        {
        public:
            static void onInitialize(const T* ptr)
            {
                if (!ptr)
                {
                    throw NullPointerException();
                }
            }

            static void onDereferenceOrIndirection(const T* ptr)
            {
                onInitialize(ptr);
            }
        };

        //(*@\Green{// Before dereferencing and indirection, the AssertOnDereferenceOrIndirection checking policy verifies in}@*)
        //(*@\Green{// debug mode the value of the stored raw pointer \_ptr.}@*)
        class AssertNullDereferenceOrIndirection
        {
        public:
            static void onInitialize(const T* /*ptr*/) {}

            static void onDereferenceOrIndirection(const T* ptr)
            {
                assert("Null pointer dereference or indirection!" && ptr);
            }
        };

        //(*@\Green{// The AssertNull checking policy verifies in debug mode the value of the stored raw pointer \_ptr}@*)
        //(*@\Green{// before every attempt of initialization/assignment, dereferencing and indirection.}@*)
        class AssertNull
        {
        public:
            static void onInitialize(const T* ptr)
            {
                assert("Initializing with null pointer is prohibited!" && ptr);
            }

            static void onDereferenceOrIndirection(const T* ptr)
            {
                assert("Null pointer dereference or indirection!" && ptr);
            }
        };
    };
}

#endif //(*@\Green{// CHECKINGPOLICIES\_H}@*)
