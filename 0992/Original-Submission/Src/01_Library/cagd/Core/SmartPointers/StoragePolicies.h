//--------------------------------------------------------------------------------------
// File:        Core/SmartPointers/StoragePolicies.h
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

#ifndef STORAGEPOLICIES_H
#define STORAGEPOLICIES_H

namespace cagd
{
    //(*@\Green{// The user can choose between the default (i.e., non-array) and array storage policies.}@*)
    template <typename T>
    class StoragePolicy
    {
    public:

        //(*@\Green{// default storage policiy}@*)
        class Default
        {
        protected:
            T* _ptr;

        public:
            //(*@\Green{// explicit special/default constructor}@*)
            explicit Default(T* ptr = nullptr);

            //(*@\Green{// copy constructor}@*)
            Default(const Default& policy);

            //(*@\Green{// overloaded dereferencing operators}@*)
            T& operator *();
            const T& operator *() const;

            //(*@\Green{// overloaded indirection operators}@*)
            T* operator ->();
            const T* operator ->() const;

        protected:
            //(*@\Green{// it will be called during policy exchanges}@*)
            void _swap(Default& policy);

            //(*@\Green{// it is responsible for the deallocation of the object referenced by the stored pointer \_ptr}@*)
            void _destroy();
        };

        //(*@\Green{// array storage policy:}@*)
        //(*@\Green{// it can only be combined with the no-copy ownership policy, since we do not know the size of the array}@*)
        //(*@\Green{// referenced by the stored raw pointer \_ptr}@*)
        class Array
        {
        protected:
            T* _ptr;

        public:
            //(*@\Green{// explicit special/default constructor}@*)
            explicit Array(T* ptr = nullptr);

            //(*@\Green{// copy constructor}@*)
            Array(const Array& policy);

            //(*@\Green{// overloaded dereferencing operators}@*)
            T& operator *();
            const T& operator *() const;

            //(*@\Green{// overloaded indexing operators}@*)
            T& operator [](int index);
            const T& operator [](int index) const;

            //(*@\Green{// overloaded indirection operators}@*)
            T* operator ->();
            const T* operator ->() const;

        protected:
            //(*@\Green{// it will be called during policy exchanges}@*)
            void _swap(Array& policy);

            //(*@\Green{// it is responsible for the deallocation of the array referenced by the stored pointer \_ptr}@*)
            void _destroy();
        };
    };

    //(*@\Green{// Implementation of the default storage policy.}@*)

    //(*@\Green{// special/default constructor}@*)
    template <typename T>
    inline StoragePolicy<T>::Default::Default(T* ptr): _ptr(ptr)
    {
    }

    //(*@\Green{// copy constructor}@*)
    template <typename T>
    inline StoragePolicy<T>::Default::Default(const Default& /*policy*/)
    {
        //(*@\Green{// the stored raw pointer \_ptr will be initialized by the clone function of the applied ownership policy}@*)
    }

    //(*@\Green{// overloaded dereferencing operators}@*)
    template <typename T>
    inline T& StoragePolicy<T>::Default::operator *()
    {
        return *_ptr;
    }

    template <typename T>
    inline const T& StoragePolicy<T>::Default::operator *() const
    {
        return *_ptr;
    }

    //(*@\Green{// overloaded indirection operators}@*)
    template <typename T>
    inline T* StoragePolicy<T>::Default::operator ->()
    {
        return _ptr;
    }

    template <typename T>
    inline const T* StoragePolicy<T>::Default::operator ->() const
    {
        return _ptr;
    }

    //(*@\Green{// it will be called during policy exchanges}@*)
    template <typename T>
    inline void StoragePolicy<T>::Default::_swap(Default& policy)
    {
        T* aux(_ptr);
        _ptr = policy._ptr;
        policy._ptr = aux;
    }

    //(*@\Green{// deallocation of the object referenced by the raw pointer \_ptr}@*)
    template <typename T>
    inline void StoragePolicy<T>::Default::_destroy()
    {
        if (_ptr)
        {
            delete _ptr, _ptr = nullptr;
        }
    }

    //(*@\Green{// Implementation of the array storage policy.}@*)

    //(*@\Green{// special/default constructor}@*)
    template <typename T>
    inline StoragePolicy<T>::Array::Array(T* ptr): _ptr(ptr)
    {
    }

    //(*@\Green{// copy constructor}@*)
    template <typename T>
    inline StoragePolicy<T>::Array::Array(const Array& /*policy*/)
    {
        //(*@\Green{// the stored raw pointer \_ptr will be initialized by the clone function of the applied ownership policy}@*)

        //(*@\Green{// do not forget: the array storage policy can only be used together with the no-copy ownership policy,}@*)
        //(*@\Green{// since we do not know the size of the array that is referenced by the stored raw pointer \_ptr}@*)
    }

    //(*@\Green{// overloaded dereferencing operators}@*)
    template <typename T>
    inline T& StoragePolicy<T>::Array::operator *()
    {
        return *_ptr;
    }

    template <typename T>
    inline const T& StoragePolicy<T>::Array::operator *() const
    {
        return *_ptr;
    }

    //(*@\Green{// overloaded indexing operators}@*)
    template <typename T>
    T& StoragePolicy<T>::Array::operator [](int index)
    {
        return _ptr[index];
    }

    template <typename T>
    inline const T& StoragePolicy<T>::Array::operator [](int index) const
    {
        return _ptr[index];
    }

    //(*@\Green{// overloaded indirection operators}@*)
    template <typename T>
    inline T* StoragePolicy<T>::Array::operator ->()
    {
        return _ptr;
    }

    template <typename T>
    inline const T* StoragePolicy<T>::Array::operator ->() const
    {
        return _ptr;
    }

    //(*@\Green{// it will be called during policy exchanges}@*)
    template <typename T>
    inline void StoragePolicy<T>::Array::_swap(Array& policy)
    {
        T* aux(_ptr);
        _ptr = policy._ptr;
        policy._ptr = aux;
    }

    //(*@\Green{// deallocation of the array referenced by the raw pointer \_ptr}@*)
    template <typename T>
    inline void StoragePolicy<T>::Array::_destroy()
    {
        if (_ptr)
        {
            delete[] _ptr, _ptr = nullptr;
        }
    }
}

#endif //(*@\Green{// STORAGEPOLICIES\_H}@*)
