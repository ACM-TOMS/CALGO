/*******************************************************************************
   
   The 'Elrint3dSequence' class provides the functions necessary to generate 
   the embedded lattice sequence used by elrint3d.  The sequence is obtained
   by copying the seed lattice (1, 230, 311)/823 three times, with intermediate
   lattices being constructed along the way.  Including the seed lattice and 
   the intermediate lattices, there is a total of ten lattices in the sequence.
   Each lattice has twice the number of points of the preceding lattice in the
   sequence.  The maximum number of points generated is 421376.

   Note
   With the exception of the static method 'getInstance', which returns the
   singleton instance of the lattice sequence, this singleton class is sealed
   by making all the methods private, including the constructor.  This is 
   because the class is based on a specific seed lattice and copying order, 
   neither of which can be altered by the user.  Access to the class is 
   delegated up to its superclass 'EmbeddedSequence<double>' and the class is 
   kept invisible by locating it in the source file. 
 
 *******************************************************************************/

// Application includes
#include    <elrint3dSequence.h>
#include    <precomputedCoords.h>

// Defines
#define  N  823  // The number of points in the seed lattice

class Elrint3dSequence : public EmbeddedSequence<double>
{
public:
    static const Elrint3dSequence* getInstance()
    // Returns singleton instance
    {
        static const Elrint3dSequence elrint3dSequence;
        return &elrint3dSequence;
    }

private:
    void get(double * const, int) const;
    // Returns the coordinates of the point indexed by 'iPoint'
    // along with the weight associated with that point. 

    int iLat(int iPoint) const;
    // Returns the lattice index of the point indexed by 'iPoint'

    int indexOfStart(int iLat) const;
    // Returns the index of the first point of the lattice indexed by 'iLat'	

    int size() const;
    // Returns the total number of points in the last lattice in the sequence
    // - fixed at 421376 for elrint3d.

    int sizeOfSeq() const;
    // Returns the total number of lattices in the sequence
    // - fixed at 10 for elrint3d.

    int count(int iLat) const;
    // Returns the number of points in the lattice indexed by 'iLat'     

    Elrint3dSequence()
    // Private constructor to prevent any direct instantiation from constructor
    {
    }

};

void Elrint3dSequence::get(double * const pa, int iPoint) const
{
    get_coords_by_index(pa, iPoint);
}

int Elrint3dSequence::iLat(int iPoint) const
{
    int n = iPoint / N;
    int i = 0;

    while (n)
    {
        i++;
        n >>= 1;
    }
    return i;
}

int Elrint3dSequence::indexOfStart(int iLat) const
{
    if (iLat == 0)
    {
        return 0;
    }
    return (1 << (iLat - 1))* N;
}

int Elrint3dSequence::size() const
{
    return 421376;
}

int Elrint3dSequence::sizeOfSeq() const
{
    return 10;
}

int Elrint3dSequence::count(int iLat) const
{
    return (1 << iLat) * N;
}

extern const EmbeddedSequence<double>* singleton_elrint3dSequence()
// Returns a global singleton instance of the elrint3d sequence
// to be shared by all cubature rules using this sequence
{
    return Elrint3dSequence::getInstance();
}
