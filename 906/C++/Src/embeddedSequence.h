#ifndef _EMBEDDED_SEQUENCE_H_
#define _EMBEDDED_SEQUENCE_H_ 

/********************************************************************************
                                                              
   This class is a generic class for generating an embedded lattice sequence.  
   It is an abstract class that cannot be instantiated by the user; it must be 
   extended to provide a specific embedded sequence.                    
                                                                               
   Notes                                                                              
   1. Template functions and parameters are used to facilitate use of various
	  arithmetic precisions; for example, float, double, long double.  The
	  data type chosen must be supported by native operations (addition,
	  subtraction, multiplication and division) and by the standard library.
                                                                               
   2. The class is not limited in dimensionality. It may also provide
	  information about points that is additional to the coordinates of the
	  points, for example, the weights associated with the points in the case
	  that a periodizing transformation has been used.
                                                                                  
   3. Like STL, the class is accessed through iterators which act like a pointer 
	  to members of the sequence.
                                                                                               
   Implementation issues                                                      
   To extend this class for a specific embedded sequence, the following pure
   virtual functions must be implemented in the inherited class:                                                          
				1. virtual int size() const
				2. virtual int sizeOfSeq() const
				3. virtual int count(int iLat) const
				4. virtual void get(T const, int iPoint ) const
				5. virtual int iLat(int iPoint) const
				6. virtual int iStart(int iLat) const

 *********************************************************************************/

// Standard include
#include	<stdexcept>

#define MAX_SEQ_VECT  256   // Reasonably assume the maximum dimension
                            // of the generated vector is 256

template <typename T>
class EmbeddedSequence
// This class is expected to be extended for a specific sequence
{
public:
	class Iterator
	// Stores the index of the current lattice in the sequence and the indices
	// of all points in the lattice.  Indices start at 0.
	{
	public:
		Iterator(const EmbeddedSequence<T>* pSequence, int iP) : pSeq(pSequence), iPoint(iP)
		// Constructor.  Iterator of 'pSequence' indexed by 'iPoint'
		{
		}

		T * operator *()
		// Returns coordinates and other information associated with the current iterator
		{
			if (iPoint < 0 || iPoint > pSeq->size())
				// Throw an exception if the index goes out of bounds
			{
				throw std::range_error("Point index out of bounds");
			}
			pSeq->get(coords, iPoint);
			return coords;
		}

		Iterator operator ++(int) // Suffix increment operator
		// Move forwards; an exception will be thrown if the index goes out of bounds
		{
			iPoint++;
			if (iPoint >= pSeq->size())
				// In this situation, stop at the end of the sequence
			{
				iPoint = pSeq->size() - 1; // Assumption: empty sequence not allowed
			}

			return *this;
		}

		Iterator operator --(int) // Suffix decrement operator
		// Move backwards; an exception will be thrown if the index goes out of bounds
		{
			iPoint--;
			if (iPoint < 0)
				// Stop at the beginning
			{
				iPoint = 0;
			}

			return *this;
		}

		bool eod() const
		// True if the iterator is at the beginning or end of the sequence
		{
			return (iPoint <= 0) || (iPoint >= pSeq->size() - 1);
		}

	private:
		const EmbeddedSequence* pSeq; // Store the sequence owning this iterator
		int iPoint; // Point index
		T coords[MAX_SEQ_VECT]; // Store current point

	};

	virtual int size() const = 0;
	// The total number of points in the final lattice of the sequence
	// (= 'count(last)', where 'last' is the index of the final lattice)

	virtual int sizeOfSeq() const = 0;
	// The total number of lattices in the sequence

	virtual int indexOfStart(int iLat) const = 0;
	// Returns the start point index of the lattice indexed by 'iLat'

	Iterator iFirst() const
	// Returns the first iterator in the sequence; used to move forwards
	{
		return Iterator(this, 0);
	}

	Iterator iLast() const
	// Returns the last iterator in the sequence; used to move backwards
	{
		return Iterator(this, size() - 1);
	}

	Iterator iStartAtSeq(int iSeq) const
	// Returns the iterator at the start of the lattice indexed by 'iLat'
	{
		return Iterator(this, indexOfStart(iSeq));
	}

	Iterator iStartAtPoint(int iPoint) const
	// Returns the iterator at the start of the point indexed by 'iPoint'
	{
		return Iterator(this, iPoint);
	}

	virtual int count(int iLat) const = 0;
	// Returns the number of points in the lattice indexed by 'iLat'

	virtual ~EmbeddedSequence()
	// Virtual destructor; required for an abstract class
	{
	}

private:
	virtual void get(T * const, int iPoint) const = 0;
	// Returns the coordinates and other information associated with the point
	// indexed by 'iPoint'. The length of output may exceed the dimensionality
	// of the point if additional information is being returned, such as an
	// associated weight following use of a periodizing transformation.

	virtual int iLat(int iPoint) const = 0;
	// Returns the index of the first lattice containing the point indexed by 'iPoint'

};

#endif  // _EMBEDDED_SEQUENCE_H_

