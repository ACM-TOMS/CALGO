//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/SmartPointers.h
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

#ifndef SMARTPOINTERS_H
#define SMARTPOINTERS_H

#include "CheckingPolicies.h"
#include "ImplicitConversionPolicies.h"
#include "OwnershipPolicies.h"
#include "StoragePolicies.h"
#include "TypeSelectors.h"

namespace cagd
{
    template
    <
        typename T,
        class    TSP  = typename StoragePolicy<T>::Default,    //(*@\Green{// default storage policy}@*)
        class    TOP  = typename OwnershipPolicy<T>::DeepCopy, //(*@\Green{// default ownership policy}@*)
        class    TICP = ImplicitConversionPolicy::Disallowed,  //(*@\Green{// default implicit conversion policy}@*)
        class    TCP  = typename CheckingPolicy<T>::NoCheck    //(*@\Green{// default checking policy}@*)
    >
    class SmartPointer: public TSP, public TOP, public TICP, public TCP
    {
    private:
        //(*@\Green{// smart pointers with different typename parameters have to be friends, otherwise we cannot perform}@*)
        //(*@\Green{// operations through the corresponding private member \_ptr}@*)
        template <typename U, class USP, class UOP, class UICP, class UCP>
        friend class SmartPointer;

    //(*@\Green{// Handling the implicit conversion policy:}@*)
    private:
        //(*@\Green{// represents a private helper class for disallowing implicit conversion}@*)
        class _Disallowed
        {
        };

        //(*@\Green{// if the implicit conversion is enabled, the \_ImplicitConversionType below will be defined as T*,}@*)
        //(*@\Green{// otherwise as \_Disallowed}@*)
        typedef typename TypeSelector<TICP::ENABLED, T*, _Disallowed>::Result
                         _ImplicitConversionType;

    public:
        //(*@\Green{// if the \_ImplicitConversionType type above coincides with \_Disallowed, then the type conversion operator}@*)
        //(*@\Green{// below will generate a compile-time error, since there exists no conversion that transforms the type T*}@*)
        //(*@\Green{// to \_Disallowed}@*)
        operator _ImplicitConversionType() const
        {
            return TSP::_ptr;
        }

    //(*@\Green{// Handling direct null pointer testing:}@*)
    private:
        //(*@\Green{// represents a private helper class for direct null pointer testing}@*)
        class NullPointerTester
        {
        private:
            void operator delete(void*); //(*@\Green{// it cannot be called}@*)
        };

    public:
        operator NullPointerTester*() const
        {
            static class NullPointerTester test;

            if (!TSP::_ptr)
            {
                return nullptr;
            }

            return &test;
        }

    //(*@\Green{// Handling the selected ownership policy:}@*)
    private:
        typedef typename
                TypeSelector< TOP::DESTRUCTIVE_COPY,
                              SmartPointer<T, TSP, TOP, TICP, TCP>,
                              const SmartPointer<T, TSP, TOP, TICP, TCP> >::Result
                _CopiedType;

        void _swap(SmartPointer& sp);

    //(*@\Green{// Declaration of special/default/copy constructors, of overloaded operators and of the destructor:}@*)
    public:
        //(*@\Green{// explicit special/default constructor}@*)
        explicit SmartPointer(T* ptr = nullptr);

        //(*@\Green{// in order to enable assignments like SmartPointer$<$const T$>$ = SmartPointer$<$T$>$,}@*)
        //(*@\Green{// we have to enable the conversion from SmartPointer$<$T$>$ to SmartPointer$<$const T$>$}@*)
        operator SmartPointer<const T, TSP, TOP, TICP, TCP>() const;

        //(*@\Green{// copy constructor}@*)
        SmartPointer(_CopiedType& sp);

        //(*@\Green{// overloaded assignment operator}@*)
        SmartPointer<T, TSP, TOP, TICP, TCP>& operator =(_CopiedType& rhs);

        //(*@\Green{// overloaded dereferencing operators}@*)
        T& operator *();
        const T& operator *() const;

        //(*@\Green{// overloaded indirection operators}@*)
        T* operator ->();
        const T* operator ->() const;

        //(*@\Green{// overloaded member comparison operators}@*)
        bool operator ==(const SmartPointer<T, TSP, TOP, TICP, TCP>& rhs) const;
        bool operator !=(const SmartPointer<T, TSP, TOP, TICP, TCP>& rhs) const;
        bool operator !() const;

        template <typename U, class USP, class UOP, class UICP, class UCP>
        bool operator ==(const SmartPointer<U, USP, UOP, UICP, UCP>& rhs) const;

        template <typename U, class USP, class UOP, class UICP, class UCP>
        bool operator !=(const SmartPointer<U, USP, UOP, UICP, UCP>& rhs) const;

        //(*@\Green{// overloaded friend comparison operators}@*)
        template <typename U, typename V, class VSP, class VOP, class VICP, class VCP>
        friend bool operator ==(const U* lhs,
                                const SmartPointer<V, VSP, VOP, VICP, VCP>& rhs);

        template <typename U, class USP, class UOP, class UICP, class UCP, typename V>
        friend bool operator ==(const SmartPointer<U, USP, UOP, UICP, UCP>& lhs,
                                const V* rhs);

        template <typename U, typename V, class VSP, class VOP, class VICP, class VCP>
        friend bool operator !=(const U* lhs,
                                const SmartPointer<V, VSP, VOP, VICP, VCP>& rhs);

        template <typename U, class USP, class UOP, class UICP, class UCP, typename V>
        friend bool operator !=(const SmartPointer<U, USP, UOP, UICP, UCP>& lhs,
                                const V* rhs);

        //(*@\Green{// destructor}@*)
        ~SmartPointer();
    };

    //(*@\Green{// special/default constructor}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline SmartPointer<T, TSP, TOP, TICP, TCP>::SmartPointer(T* ptr):
           TSP(ptr), TOP(), TICP(), TCP()
    {
        TCP::onInitialize(ptr);
    }

    //(*@\Green{// in order to enable assignments like SmartPointer$<$const T$>$ = SmartPointer$<$T$>$,}@*)
    //(*@\Green{// we have to enable the conversion from SmartPointer$<$T$>$ to SmartPointer$<$const T$>$}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline SmartPointer<T, TSP, TOP, TICP, TCP>::operator
           SmartPointer<const T, TSP, TOP, TICP, TCP>() const
    {
        return SmartPointer<const T, TSP, TOP, TICP, TCP>(TOP::clone(TSP::_ptr));
    }

    //(*@\Green{// copy constructor}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline SmartPointer<T, TSP, TOP, TICP, TCP>::SmartPointer(_CopiedType& sp):
           TSP(sp), TOP(sp), TICP(), TCP()
    {
        TSP::_ptr = TOP::_clone(sp._ptr);
    }

    //(*@\Green{// exchanging storage and ownership policies}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline void SmartPointer<T, TSP, TOP, TICP, TCP>::_swap(
           SmartPointer<T, TSP, TOP, TICP, TCP>& sp)
    {
        TSP::_swap(sp);
        TOP::_swap(sp);
    }

    //(*@\Green{// overloaded assignment operator}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline SmartPointer<T, TSP, TOP, TICP, TCP>&
           SmartPointer<T, TSP, TOP, TICP, TCP>::operator =(_CopiedType& rhs)
    {
        if (this != &rhs)
        {
            if (TOP::_release(TSP::_ptr))
            {
                TSP::_destroy();
            }

            SmartPointer aux(rhs);
            aux._swap(*this);
        }

        return *this;
    }

    //(*@\Green{// overloaded dereferencing operators}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline T& SmartPointer<T, TSP, TOP, TICP, TCP>::operator *()
    {
        TCP::onDereferenceOrIndirection(TSP::_ptr);
        return TSP::operator *();
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline const T& SmartPointer<T, TSP, TOP, TICP, TCP>::operator *() const
    {
        TCP::onDereferenceOrIndirection(TSP::_ptr);
        return TSP::operator *();
    }

    //(*@\Green{// overloaded indirection operators}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline T* SmartPointer<T, TSP, TOP, TICP, TCP>::operator ->()
    {
        TCP::onDereferenceOrIndirection(TSP::_ptr);
        return TSP::operator ->();
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline const T* SmartPointer<T, TSP, TOP, TICP, TCP>::operator ->() const
    {
        TCP::onDereferenceOrIndirection(TSP::_ptr);
        return TSP::operator ->();
    }

    //(*@\Green{// overloaded member comparison operators}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline bool SmartPointer<T, TSP, TOP, TICP, TCP>::operator ==(
           const SmartPointer<T, TSP, TOP, TICP, TCP>& rhs) const
    {
        return (TSP::_ptr == rhs._ptr);
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline bool SmartPointer<T, TSP, TOP, TICP, TCP>::operator !=(
           const SmartPointer<T, TSP, TOP, TICP, TCP>& rhs) const
    {
        return (TSP::_ptr != rhs._ptr);
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline bool SmartPointer<T, TSP, TOP, TICP, TCP>::operator !() const
    {
        return !TSP::_ptr;
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    template <typename U, class USP, class UOP, class UICP, class UCP>
    inline bool SmartPointer<T, TSP, TOP, TICP, TCP>::operator ==(
           const SmartPointer<U, USP, UOP, UICP, UCP>& rhs) const
    {
        return ((void*)TSP::_ptr == (void*)rhs._ptr);
    }

    template <typename T, class TSP, class TOP, class TICP, class TCP>
    template <typename U, class USP, class UOP, class UICP, class UCP>
    inline bool SmartPointer<T, TSP, TOP, TICP, TCP>::operator !=(
           const SmartPointer<U, USP, UOP, UICP, UCP>& rhs) const
    {
        return ((void*)TSP::_ptr != (void*)rhs._ptr);
    }

    //(*@\Green{// overloaded friend comparison operators}@*)
    template <typename U, typename V, class VSP, class VOP, class VICP, class VCP>
    inline bool operator ==(
            const U* lhs, const SmartPointer<V, VSP, VOP, VICP, VCP>& rhs)
    {
        return ((void*)lhs == (void*)rhs._ptr);
    }

    template <typename U, class USP, class UOP, class UICP, class UCP, typename V>
    inline bool operator ==(
            const SmartPointer<U, USP, UOP, UICP, UCP>& lhs, const V* rhs)
    {
        return ((void*)lhs._ptr == (void*)rhs);
    }

    template <typename U, typename V, class VSP, class VOP, class VICP, class VCP>
    inline bool operator !=(
            const U* lhs, const SmartPointer<V, VSP, VOP, VICP, VCP>& rhs)
    {
        return ((void*)lhs != (void*)rhs._ptr);
    }

    template <typename U, class USP, class UOP, class UICP, class UCP, typename V>
    inline bool operator !=(
            const SmartPointer<U, USP, UOP, UICP, UCP>& lhs, const V* rhs)
    {
        return ((void*)lhs._ptr != (void*)rhs);
    }

    //(*@\Green{// destructor}@*)
    template <typename T, class TSP, class TOP, class TICP, class TCP>
    inline SmartPointer<T, TSP, TOP, TICP, TCP>::~SmartPointer()
    {
        if (TOP::_release(TSP::_ptr))
        {
            TSP::_destroy();
        }
    }
}

#endif //(*@\Green{// SMARTPOINTERS\_H}@*)
