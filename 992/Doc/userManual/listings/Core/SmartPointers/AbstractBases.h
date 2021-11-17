#ifndef ABSTRACTBASES_H
#define ABSTRACTBASES_H

namespace cagd
{
    class AbstractBase
    {
    public:
        (*@\Green{// pure virtual method that has to redeclared and defined in all derived classes}@*)
        virtual AbstractBase* clone() const = 0;

        (*@\Green{// virtual destructor}@*)
        virtual ~AbstractBase() {}
    };
}

#endif (*@\Green{// ABSTRACTBASES\_H}@*)
