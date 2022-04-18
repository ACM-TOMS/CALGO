//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/SpecializedSmartPointers.h
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

#ifndef SPECIALIZEDSMARTPOINTERS_H
#define SPECIALIZEDSMARTPOINTERS_H

#include "SmartPointers.h"

namespace cagd
{
    template <typename T>
    struct SP
    {
        typedef
        SmartPointer
        <
            T,
            typename StoragePolicy<T>::Default,
            typename OwnershipPolicy<T>::DeepPrimitiveCopy,
            ImplicitConversionPolicy::Disallowed,
            typename CheckingPolicy<T>::RejectNullDereferenceOrIndirection
        >
        DefaultPrimitive;

        typedef
        SmartPointer
        <
            T,
            typename StoragePolicy<T>::Default,
            typename OwnershipPolicy<T>::DeepCopy,
            ImplicitConversionPolicy::Disallowed,
            typename CheckingPolicy<T>::RejectNullDereferenceOrIndirection
        >
        Default;

        typedef
        SmartPointer
        <
            T,
            typename StoragePolicy<T>::Array,
            typename OwnershipPolicy<T>::NoCopy,
            ImplicitConversionPolicy::Disallowed,
            typename CheckingPolicy<T>::RejectNullDereferenceOrIndirection
        >
        Array;

        typedef
        SmartPointer
        <
            T,
            typename StoragePolicy<T>::Default,
            typename OwnershipPolicy<T>::DestructiveCopy,
            ImplicitConversionPolicy::Disallowed,
            typename CheckingPolicy<T>::RejectNullDereferenceOrIndirection
        >
        DestructiveCopy;

        typedef
        SmartPointer
        <
            T,
            typename StoragePolicy<T>::Default,
            typename OwnershipPolicy<T>::NonIntrusiveReferenceCounting,
            ImplicitConversionPolicy::Disallowed,
            typename CheckingPolicy<T>::RejectNullDereferenceOrIndirection
        >
        NonIntrusiveReferenceCounting;
    };
}

#endif (*@\Green{// SPECIALIZEDSMARTPOINTERS\_H}@*)
