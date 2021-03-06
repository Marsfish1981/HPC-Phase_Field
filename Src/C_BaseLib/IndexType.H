
#ifndef BL_INDEXTYPE_H
#define BL_INDEXTYPE_H

#include <iosfwd>

#include <ccse-mpi.H>
#include <IntVect.H>
#include <SPACE.H>

//
// Cell-Based or Node-Based Indices
//
// The class IndexType defines an index as being cell based or node (edge)
// based in each of the BL_SPACEDIM directions.  This class defines an
// enumerated type CellIndex to be either CELL or NODE; i.e. each of the
// BL_SPACEDIM dimensions must be either CELL or NODE.
//

class IndexType
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<IndexType>::type();

public:
    //
    // The cell index type: one of CELL or NODE.
    //
    enum CellIndex { CELL = 0, NODE = 1 };
    //
    // The default constructor
    //
    IndexType ()
        :
        itype(0) {}
    //
    // Construct an IndexType identical to an IntVect.
    //
    explicit IndexType (const IntVect& iv)
    {
        itype = D_TERM((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2));
    }
    //
    // Construct an IndexType given an explicit CellIndex for
    // each direction.  D_DECL is a macro that sets the constructor
    // to take BL_SPACEDIM arguments.
    //
    IndexType (D_DECL(CellIndex i, CellIndex j, CellIndex k))
    {
        itype = D_TERM(i, | (j<<1), | (k<<2));
    }
    //
    // Set IndexType to be NODE based in direction dir.
    //
    void set (int dir) { itype |= mask(dir); }
    //
    // Set IndexType to be CELL based in direction dir.
    //
    void unset (int dir) { itype &= ~mask(dir); }
    //
    // True if IndexType is NODE based in direction dir.
    //
    bool test (int dir) const { return (itype & mask(dir)) != 0; }
    //
    // Set NODE based in all directions.
    //
    void setall () { itype = (1 << BL_SPACEDIM) - 1; }
    //
    // Set CELL based in all directions.
    //
    void clear () { itype = 0; }
    //
    // True if this IndexType is NODE based in any direction.
    //
    bool any () const { return itype != 0; }
    //
    // True if IndexType is valid.
    //
    bool ok () const { return itype < (1 << BL_SPACEDIM); }
    //
    // Change from CELL to NODE or NODE to CELL in direction dir.
    //
    void flip (int i) { itype ^= mask(i); }
    //
    // True if IndexTypes are identical.
    //
    bool operator== (const IndexType& t) const { return t.itype == itype; }
    //
    // True if IndexTypes are not identical.
    //
    bool operator!= (const IndexType& t) const { return t.itype != itype; }
    //
    // True if the IndexType is CELL based in all directions.
    //
    bool cellCentered () const { return itype == 0; }
    //
    // True if the IndexType is CELL based in dir-direction.
    //
    bool cellCentered (int dir) const { return (itype & mask(dir)) == 0; }
    //
    // True if the IndexType is NODE based in all directions.
    //
    bool nodeCentered () const { return itype == (1<<BL_SPACEDIM)-1; }
    //
    // True if the IndexType is NODE based in dir-direction.
    //
    bool nodeCentered (int dir) const { return (itype & mask(dir)) != 0; }
    //
    // Set IndexType to CellIndex type t in direction dir.
    //
    void setType (int dir, CellIndex t);
    //
    // Returns  the CellIndex in direction dir.
    //
    CellIndex ixType (int dir) const { return (CellIndex) ((itype & (1<<dir)) >> dir); }
    //
    // Return an integer representing the IndexType in direction dir.
    //
    int operator[] (int dir) const { return test(dir); }
    //
    // Fill an IntVect of size BL_SPACEDIM with IndexTypes.
    //
    IntVect ixType () const { return IntVect(D_DECL(itype&1, (itype>>1)&1, (itype>>2)&1)); }
    //
    // This static member function returns an IndexType object of value
    // IndexType::CELL.  It is provided as a convenience to our users
    // when defining a Box all of whose faces should be of type
    // IndexType::CELL.
    //
    static IndexType TheCellType ();
    //
    // This static member function returns an IndexType object of value
    // IndexType::NODE.  It is provided as a convenience to our users
    // when defining a Box all of whose faces should be of type
    // IndexType::NODE.
    //
    static IndexType TheNodeType ();

private:
    //
    // Returns 1<<k.
    //
    static int mask (int k) { return 1<<k; }
    //
    // An integer holding the CellIndex in bits 0 - BL_SPACEDIM-1.
    //
    unsigned int itype;
};

//
// Write an IndexType to an ostream in ASCII.
//
std::ostream& operator<< (std::ostream& os, const IndexType& itype);
//
// Read an IndexType from an istream.
//
std::istream& operator>> (std::istream& is, IndexType& itype);

#endif /*BL_INDEXTYPE_H*/
