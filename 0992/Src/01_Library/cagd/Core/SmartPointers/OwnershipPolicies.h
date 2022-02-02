//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/OwnershipPolicies.h
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

#ifndef OWNERSHIPPOLICIES_H
#define OWNERSHIPPOLICIES_H

#include "StaticChecks.h"
#include "TypeSelectors.h"
#include <new>

namespace cagd
{
    //(*@\Green{// Currently one can choose from the no-copy, deep copy, deep primitive copy, destructive copy and}@*)
    //(*@\Green{// non-intrusive reference counting ownership policies.}@*)
    template <typename T>
    class OwnershipPolicy
    {
    public:
        //(*@\Green{// no copy ownership policy}@*)
        class NoCopy
        {
        public:
            enum {DESTRUCTIVE_COPY = false};

        protected:
            T*   _clone(T* const& ptr);
            bool _release(T* const& ptr);
            void _swap(NoCopy& policy);
        };

        //(*@\Green{// deep copy ownership policy for custom types}@*)
        class DeepCopy
        {
        public:
            enum {DESTRUCTIVE_COPY = false};

        protected:
            T*   _clone(T* const& ptr);
            bool _release(T* const& ptr);
            void _swap(DeepCopy& policy);
        };

        //(*@\Green{// deep copy ownership policy for primitive types}@*)
        class DeepPrimitiveCopy
        {
        public:
            enum {DESTRUCTIVE_COPY = false};

        protected:
            T*   _clone(T* const& ptr);
            bool _release(T* const& ptr);
            void _swap(DeepPrimitiveCopy& policy);
        };

        //(*@\Green{// destructive copy ownership policy}@*)
        class DestructiveCopy
        {
        public:
            enum {DESTRUCTIVE_COPY = true};

        protected:
            T*   _clone(T* & ptr);
            bool _release(T* const& ptr);
            void _swap(DestructiveCopy& policy);
        };

        //(*@\Green{// non-intrusive reference counting ownership policy}@*)
        class NonIntrusiveReferenceCounting
        {
        private:
            size_t *_reference_count; //(*@\Green{// records the total number of references}@*)

        public:
            enum {DESTRUCTIVE_COPY = false};

            //(*@\Green{// default constructor}@*)
            NonIntrusiveReferenceCounting();

        protected:
            T*   _clone(T* const& ptr);
            bool _release(T* const& ptr);
            void _swap(NonIntrusiveReferenceCounting& policy);
        };
    };

    //(*@\Green{// implementation of the no-copy ownership policy}@*)
    template <typename T>
    inline T* OwnershipPolicy<T>::NoCopy::_clone(T* const& ptr)
    {
        STATIC_CHECK(false, The_no_copy_ownership_policy_disallows_copying);
    }

    template <typename T>
    inline bool OwnershipPolicy<T>::NoCopy::_release(T* const& ptr)
    {
        return true;
    }

    template <typename T>
    inline void OwnershipPolicy<T>::NoCopy::_swap(NoCopy& policy)
    {
    }

    //(*@\Green{// implementation of the deep copy ownership policy for custom types}@*)
    template <typename T>
    inline T* OwnershipPolicy<T>::DeepCopy::_clone(T* const& ptr)
    {
        return (ptr ? ptr->clone() : nullptr);
    }

    template <typename T>
    inline bool OwnershipPolicy<T>::DeepCopy::_release(T* const& /*ptr*/)
    {
        return true;
    }

    template <typename T>
    inline void OwnershipPolicy<T>::DeepCopy::_swap(DeepCopy& /*policy*/)
    {
    }

    //(*@\Green{// implementation of the deep copy ownership policy for primitive types}@*)
    template <typename T>
    inline T* OwnershipPolicy<T>::DeepPrimitiveCopy::_clone(T* const& ptr)
    {
        if (ptr && primitive(*ptr))
        {
            return new (std::nothrow) T(*ptr);
        }

        return nullptr;
    }

    template <typename T>
    inline bool OwnershipPolicy<T>::DeepPrimitiveCopy::_release(T* const& ptr)
    {
        return true;
    }

    template <typename T>
    inline void OwnershipPolicy<T>::DeepPrimitiveCopy::_swap(DeepPrimitiveCopy& policy)
    {
    }

    //(*@\Green{// implementation of the destructive copy ownership policy:}@*)
    //(*@\Green{// during assignment or copying the ownership is transfered to the target smart pointer and the raw pointer}@*)
    //(*@\Green{// stored by the copied pointer is set to null}@*)
    template <typename T>
    inline T* OwnershipPolicy<T>::DestructiveCopy::_clone(T* & ptr)
    {
        T* aux(ptr);
        ptr = nullptr;
        return aux;
    }

    template <typename T>
    inline bool OwnershipPolicy<T>::DestructiveCopy::_release(T* const& /*ptr*/)
    {
        return true;
    }

    template <typename T>
    inline void OwnershipPolicy<T>::DestructiveCopy::_swap(DestructiveCopy& /*policy*/)
    {
    }

    //(*@\Green{// implementation of the non-intrusive reference counting ownership policy:}@*)
    //(*@\Green{// this ownership policy allows the sharing of a dynamically allocated object's memory address by multiple}@*)
    //(*@\Green{// smart pointers, records the number of total references, and destroys the referenced object when the total}@*)
    //(*@\Green{// reference count becomes zero}@*)
    template <typename T>
    inline OwnershipPolicy<T>::NonIntrusiveReferenceCounting::
        NonIntrusiveReferenceCounting():
        _reference_count(new (std::nothrow) size_t(1))
    {
    }

    template <typename T>
    inline T* OwnershipPolicy<T>::NonIntrusiveReferenceCounting::_clone(T* const& ptr)
    {
        if (_reference_count)
        {
            *_reference_count += 1; //(*@\Green{// during copying the total reference count is increased}@*)
        }

        return ptr;
    }

    template <typename T>
    inline bool OwnershipPolicy<T>::NonIntrusiveReferenceCounting::_release(
                T* const& ptr)
    {
        if (_reference_count)
        {
            //(*@\Green{// when a reference counting smart pointer goes out of scope the total reference count is decreased}@*)
            *_reference_count -= 1;

            //(*@\Green{// if the total reference count becomes 0, we delete the reference tracking pointer and we indicate that}@*)
            //(*@\Green{// the stored raw pointer should also be deleted}@*)
            if (*_reference_count == 0)
            {
                delete _reference_count, _reference_count = nullptr;

                return true;
            }
        }

        //(*@\Green{// if the total reference count is not 0, the stored raw pointer should be not deleted}@*)
        return false;
    }

    template <typename T>
    inline void OwnershipPolicy<T>::NonIntrusiveReferenceCounting::_swap(
            NonIntrusiveReferenceCounting& policy)
    {
        //(*@\Green{// The function \_swap will be called by a temporary smart pointer in the assignment operator of the}@*)
        //(*@\Green{// class SmartPointer that will be discussed later in Listing \mref{src:SmartPointers.h}.}@*)
        //(*@\Green{// After the exchange of the policies, the aforementioned temporary smart pointer will go out of scope}@*)
        //(*@\Green{// and consequently the number of total references will be further decreased by $1$. In order to avoid this,}@*)
        //(*@\Green{// initially we increase the number of total referecenes by $1$.}@*)
        if (policy._reference_count)
        {
            *policy._reference_count += 1;
        }

        size_t *aux             = _reference_count;
        _reference_count        = policy._reference_count;
        policy._reference_count = aux;
    }
}

#endif //(*@\Green{// OWNERSHIPPOLICIES\_H}@*)
