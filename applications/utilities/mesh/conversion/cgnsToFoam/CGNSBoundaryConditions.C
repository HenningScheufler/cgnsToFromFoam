/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright Hydro-Quebec - IREQ, 2008
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Authors 
  Robert Magnan,  Hydro-Quebec - IREQ, 2007

\*---------------------------------------------------------------------------*/

#include "CGNSBoundaryConditions.H"
#include "CGNSElementType.H"
#include "FoamReplaceChar.H"
#include "HashTable.H"
#include "error.H"
#include "mathematicalConstants.H"

#include <limits>
#include "libcgnsoo3/cgnsoo.H"

CGNSBoundaryConditions::CGNSBoundaryConditions( const ConnectivityMapper& m, 
const map<std::string,CGNSOO::BCType_t>& fammap,
bool debug ) : 
    mapper( m ), familymap(fammap), debug_(debug)
{
    identifyBoundaryFaces();
}

static std::vector<foamNodeIndex> connectivityToPointList( const std::vector<int>& connec )
{
    std::set<foamNodeIndex> s;
    for ( std::vector<int>::const_iterator i=connec.begin() ; i!=connec.end() ; i++ )
        s.insert( static_cast<foamNodeIndex>(*i) );
    std::vector<foamNodeIndex> v(s.size());
    int k = 0;
    for ( std::set<foamNodeIndex>::const_iterator i=s.begin() ; i!=s.end() ; i++, k++ )
        v[k] = *i;
    return v;
}

void CGNSBoundaryConditions::addPatchesFromElements( CGNSOO::BCType_t bctype  )
{
    std::vector<ConnectivityMapper::zone2d_data_t> sections = mapper.getSectionsInfo();
    for ( size_t i=0 ; i<sections.size() ; i++ )
    {
        std::string sectionname  = sections[i].name;
        CGNSOO::ElementType_t etype = sections[i].etype;
        int indexZone       = sections[i].zindex;

        // convert the connectivity data into a list of points
        // we will figure out the corresponding faces later
        std::vector<int> connec = sections[i].connec;
        int expected_faces = connec.size() / CGNSElementType::getNumberOfNodesPerElement(etype);

        std::vector<foamNodeIndex> bc_point_list = connectivityToPointList( connec );

        // Convert CGNS indices to Foam indices using the zone mapper
        int k = 0;
        std::vector<foamNodeIndex> foam_point_list( bc_point_list.size() );
        for ( std::vector<int>::iterator i  = bc_point_list.begin() ;
              i != bc_point_list.end() ;
              i++ )
        {
            foamNodeIndex index_local = (*i)-1; // CGNS indices start at 1
            foam_point_list[k++] = mapper.n_localToGlobal(indexZone,index_local);
        }

        // Remember the BC name in a std::set, so we can figure out how many DISTINCT one there are
        //Foam::Info << "Defining a new patch of " << foam_point_list.size() 
        //     << " points named " << sectionname << endl;
        s_bcname.insert( sectionname );

        // store the patch information (type & pointlist) 
        // in a multimap using the section name as a key
        bcpatches.insert( mm_bcp::value_type(sectionname, patchinfo(bctype,foam_point_list,expected_faces) ) );
    }
}

void CGNSBoundaryConditions::addBoundaryPatch( CGNSOO::Base_t     base, 
int        indexZone, 
CGNSOO::ZoneType_t zonetype, 
CGNSOO::BC_t       bc, 
std::string     bcname, 
CGNSOO::BCType_t   bctype, 
CGNSOO::PointSetType_t psettype )
{
    // On elimine de facto les types non valides
    bool bc_name_ok = true;
    if ( bctype == CGNSOO::FamilySpecified )
    {
        std::string familyName;
        bc.readFamilyName(familyName);
        map<std::string, CGNSOO::BCType_t>::const_iterator ifam = familymap.find( familyName );
        if ( ifam == familymap.end() )
        {
            Foam::Warning << "Family " << familyName 
                << " was referenced in zone " << indexZone 
                << " but has no BC associated with it" << Foam::endl;
            bc_name_ok = false;
        }
        else
        {
            bctype = ifam->second;
            if ( debug_ )
                Foam::Info << "BC named " << bcname 
                    << " is family specified." << Foam::endl
                    << "It is mapped to new name " << familyName 
                    << " with bctype = " << bctype << Foam::endl;
            bcname = familyName;
        }
    }
    else if ( bctype == UserDefined )
    {
        std::string familyName;
        try 
        {
            bc.readFamilyName(familyName);
        }
        catch( CGNSOO::cgns_exception& e )
        {
            familyName = "undefined";
        }

        // pre-filter Family_Name of type Orphan or other
        bc_name_ok = validateUserDefinedBCName(bcname + "_" + familyName);  
    }
    else if ( bctype == CGNSOO::BCTypeNull )
    {
        Foam::Warning << "The BC " << bcname 
            << " is defined with a type BCTypeNull." 
            << Foam::endl
            << "It will be mapped to Foam's wallFunction type"
            << Foam::endl;
    }

    if ( ! bc_name_ok ) return;

    // map CGNS bc names to names acceptable by Foam 'word' type
    std::string fbcname = FoamReplaceChar(bcname,' ','_');

    // vector of points indices (CGNS style, >=1) defining this patch
    std::vector<int> bc_point_list;

    // used for internal validation
    // = number of faces expected in a patch defined by range.
    int expected_faces = 0;

    // We should make sure the GridLocation is set to Vertex,
    // since faceCenter is not handled here (yet)
    CGNSOO::GridLocation_t gloc = CGNSOO::Vertex;   
    try
    {
        bc.readGridLocation( gloc );
    }
    catch( CGNSOO::cgns_notfound& ) {}

    switch( gloc )
    {
    case CGNSOO::Vertex:
            break; //ok
    case CGNSOO::FaceCenter:
            Foam::FatalError << "BC defined by a GridLocation=FaceCenter"
                << " - not yet supported"
                << exit(Foam::FatalError);
            break;
        case CGNSOO::IFaceCenter:
            Foam::FatalError << "BC defined by a GridLocation=IFaceCenter"
                << " - not yet supported"
                << exit(Foam::FatalError);
    case CGNSOO::JFaceCenter:
            Foam::FatalError << "BC defined by a GridLocation=JFaceCenter"
                << " - not yet supported"
                << exit(Foam::FatalError);
    case CGNSOO::KFaceCenter: 
            Foam::FatalError << "BC defined by a GridLocation=KFaceCenter"
                << " - not yet supported"
                << exit(Foam::FatalError);
    case CGNSOO::GridLocationNull:
    case CGNSOO::GridLocationUserDefined:
    case CGNSOO::CellCenter:
    case CGNSOO::EdgeCenter:
            Foam::FatalError << "BC defined by an invalid GridLocation"
                << " - not yet supported"
                << exit(Foam::FatalError);
    }

    // get the BC point indices and map them to the new numbering
    switch( psettype )
    {
    case CGNSOO::PointRange:
            {
                CGNSOO::range prange;
        bc.readPointRange( prange );
        switch( zonetype )
        {
                    case CGNSOO::Structured:
            {
                            int i_min = prange.low(0);  int i_max = prange.high(0);
                            int j_min = prange.low(1);  int j_max = prange.high(1);
                            int k_min = prange.low(2);  int k_max = prange.high(2);

                            if ( (i_min != i_max && j_min != j_max && k_min != k_max) 
                            || (i_min > i_max) 
                            || (j_min > j_max) 
                            || (k_min > k_max) 
                            )
                            {
                Foam::FatalError << "Invalid BC " << bcname << Foam::endl
                                    << "range min[" << i_min << ":" << j_min << ":" << k_min << "]" << Foam::endl
                                    << "range max[" << i_max << ":" << j_max << ":" << k_max << "]" << Foam::endl
                                    << exit( Foam::FatalError );
                            }
                            else
                            {
                expected_faces = 1;
                if ( i_min != i_max ) expected_faces *= (i_max-i_min);
                if ( j_min != j_max ) expected_faces *= (j_max-j_min);
                if ( k_min != k_max ) expected_faces *= (k_max-k_min);

                int ni=1, nj=1, nk=1;
                try
                {
                                    mapper.getStructuredZoneDimensions(indexZone,ni,nj,nk);
                }
                catch( std::logic_error& e )
                {
                                    FatalErrorIn("CGNSBoundaryConditions::addBoundaryPatch")
                                        << "Cannot define a BC by PointRange in an unstructured zone"
                                            << "This should not happen!!" 
                                            << Foam::endl
                                            << abort( Foam::FatalError );
                }

                // build a point list representing this range
                // dimension of the point indices list
                int nbindices = (i_max-i_min+1) * (j_max-j_min+1) * (k_max-k_min+1);
                bc_point_list.resize( nbindices );

                int count = 0;
                for( int i=i_min ; i<=i_max ; i++ )
                                    for( int j=j_min ; j<=j_max ; j++ )
                                        for( int k=k_min ; k<=k_max ; k++ )
                                        {
                                            int cgnsindex = mapper.getLocalNodeIndexFromStructuredCGNSTriplet( ni,nj,nk, i,j,k );
                                            bc_point_list[count++] = cgnsindex; // CGNS index starts at 1
                                        }
                if ( debug_ )
                                    Foam::Info << "BC " << bcname 
                    << " of type " << getFoamTypeName(bctype,bcname) 
                    << ", defined with a range" << Foam::endl
                    << "\tBC range min: " << i_min << " : " << j_min << " : " << k_min << Foam::endl
                    << "\tBC range max: " << i_max << " : " << j_max << " : " << k_max << Foam::endl
                    << "\texpecting " << expected_faces << " faces" << Foam::endl
                    << "\tnbindices = " << nbindices << ", count = " << count << Foam::endl
                    << "zone dimensions = " << ni << "x" << nj << "x" << nk << Foam::endl;
                            }
            }
            break;
                    case CGNSOO::Unstructured:
            {
                            // Region defined with a point range 
                            // in an unstructured zone
                            // The point range has only two indices
                            int i_min = prange.low(0);  
                            int i_max = prange.high(0);

                            if ( i_min > i_max )
                            {
                Foam::FatalError << "Invalid BC " << bcname 
                                    << "range min=" << i_min 
                                    << ", max= " << i_max 
                                    << exit(Foam::FatalError);
                            }
                            else
                            {
                // build a point list representing this range
                // dimension of the point indices list
                int nbindices = (i_max-i_min+1);
                bc_point_list.resize( nbindices );

                int count = 0;
                for( int cgnsindex=i_min ; cgnsindex<=i_max ; cgnsindex++ )
                {
                                    bc_point_list[count++] = cgnsindex;
                }
                if ( debug_ )
                                    Foam::Info << "BC " << bcname 
                    << " of type " << getFoamTypeName(bctype,bcname) 
                    << ", defined with a range" << Foam::endl
                    << "\tBC range min=" << i_min << ", max=" << i_max << Foam::endl
                    << "\tnbindices = " << nbindices << ", count = " << count << Foam::endl;
                            }
            }
            break;
                    default:
            Foam::FatalError << "Invalid zone type" 
                            << exit(Foam::FatalError);
        } // switch zonetype
            } // case PointRange
            break;

    case CGNSOO::PointList: 
            bc.readPointList( bc_point_list );
            break;

    default:
            Foam::Warning << "BC " << bcname
                << " : patch defined by a unsupported PointSetType_t" 
                << " - skipped."
                << Foam::endl;
            break;
    }   

    // Convert CGNS indices to Foam indices using the zone mapper
    int k = 0;
    std::vector<foamNodeIndex> foam_point_list( bc_point_list.size() );
    for ( std::vector<int>::iterator i  = bc_point_list.begin() ;
          i != bc_point_list.end() ;
          i++ )
    {
        foamNodeIndex index_local = (*i)-1; // CGNS indices start at 1
        foam_point_list[k++] = mapper.n_localToGlobal(indexZone,index_local);
    }

    // Remember the BC name in a std::set, so we can figure out how many DISTINCT one there are
    s_bcname.insert( fbcname );

    // store the patch information (type & pointlist) 
    // in a multimap using fbcname as a key
    bcpatches.insert( mm_bcp::value_type(fbcname, patchinfo(bctype,foam_point_list,expected_faces) ) );
}

// Builds the list of faces defining the boundary of the domain.
// These are the faces that are referenced by only one cell. 
// This list is used in the lookup of boundary condition patches.
void CGNSBoundaryConditions::identifyBoundaryFaces()
{
    const Foam::cellShapeList& cellsAsShapes = mapper.getCells();

    /* the following code was strongly inspired from 
        polyMesh::polyMesh
    (
    const IOobject& io,
    const pointField& points,
    const cellShapeList& cellsAsShapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const wordList& boundaryPatchTypes,
    const word& defaultBoundaryPatchType,
    const wordList& boundaryPatchPhysicalTypes
    )
    see src/OpenFOAM/meshes/polyMesh/polyMeshFromShapeMesh.C
        */

    // Set up a list of face shapes for each cell
    if ( debug_ )
    Foam::Info << "Allocating cellFaceList - size = " << cellsAsShapes.size() << Foam::endl;
    Foam::faceListList cellFaceList(cellsAsShapes.size()); // the list of faces making up each cell
    if ( debug_ )
    Foam::Info << "Allocating cellFaces - size = " << cellsAsShapes.size() << Foam::endl;
    Foam::labelListList cellFaces(cellsAsShapes.size());   // the list of unique face ids making up each cell

    if ( debug_ )
    Foam::Info << "Prefilling cellFaces for each cell" << Foam::endl;
    forAll(cellsAsShapes, cellI)
    {
        cellFaceList[cellI] = cellsAsShapes[cellI].faces();

    // Count maximum possible numer of mesh faces
        //maxFaces += cellsFaceShapes[cellI].size();
    Foam::label nCellFaces = cellFaceList[cellI].size();
        cellFaces[cellI].setSize( nCellFaces );

        // Initialise cell faces ids to -1 to flag undefined faces
        cellFaces[cellI] = -1; // this sets the whole list to -1
    }

    // setup point-to-cell addressing
    if ( debug_ )
    Foam::Info << "Setup point-to-cell adressing" << Foam::endl;
    Foam::List<Foam::DynamicList<Foam::label, Foam::primitiveMesh::cellsPerPoint_> > 
        tmp_pc(mapper.getTotalNodes());
    forAll(cellsAsShapes, i)
    {
        const Foam::labelList& labels = cellsAsShapes[i];

        forAll(labels, j)
        {
            // Enter the cell label in the point's cell list
            Foam::label curPoint = labels[j];
            tmp_pc[curPoint].append(i);
        }
    }
    Foam::labelListList pointToCells(tmp_pc.size());
    forAll (pointToCells, pointI)
    {
        pointToCells[pointI].transfer( tmp_pc[pointI].shrink() );
    }

    // Build the list of boundary faces
    // For each face of each cell, we check if the neighbour cells
    // contains an identical face. We avoid scanning internal faces twice
    // by looking only at neighbours with a bigger index than the current cell.
    if ( debug_ )
    Foam::Info << "Lookup boundary faces" << Foam::endl;
    Foam::label faceIndex = 0;
    forAll(cellsAsShapes, cellI)
    {
        const Foam::faceList& curFaces = cellFaceList[cellI];

        // For all faces of this cell
        forAll(curFaces, faceI)
        {
            // Skip faces that have already been matched
            if (cellFaces[cellI][faceI] >= 0) continue;

            bool found = false;

            const Foam::face& curFace = curFaces[faceI];

            const Foam::labelList& curNeighbours = 
                pointToCells[curFace[0]];

            // For all neighbour cells
            forAll(curNeighbours, neiI)
            {
                Foam::label curNei = curNeighbours[neiI];

                // Reject neighbours with the lower label
                // they have already been scanned
                if (curNei > cellI)
                {
                    // Get the list of search faces
                    const Foam::faceList& neiFaces = cellFaceList[curNei];

                    forAll(neiFaces, neiFaceI)
                    {
                        if (neiFaces[neiFaceI] == curFace)
                        {
                            // Match found!!
                            found = true;
                            cellFaces[curNei][neiFaceI] = faceIndex;
                            cellFaces[cellI][faceI]     = faceIndex;
                            faceIndex++;
                            break;
                        }
                    }
                    if (found) break;
                }
                if (found) break;
            }

        if ( !found )
        {
            // No neighbour was found so this must be a boundary face
        boundary_faces.append( curFace );
        }
        }  // End of current faces
    }  // End of current cells

    if ( debug_ )
    Foam::Info << "Shrinking tables" << Foam::endl;
    boundary_faces.shrink();
}

// local class behaving like a std::set but which also keeps track of
// the number of time an element has been added to the set
template <typename T>
struct hash
{
public:
    hash(){}
    size_t operator()( const T& t ) const 
    {
            register size_t hashVal = 0;
            const Foam::labelList& labels = t;
            forAll( labels, i )
            {
                hashVal += labels[i];
            }
            return hashVal;
    }
    size_t operator() ( const T& t, const size_t tableSize ) const
    {
            return ::abs(operator()(t)) % tableSize;
    }
};

template <typename T>
class CountedSet : public Foam::HashTable<unsigned,T,hash<T> >
{
    CountedSet( CountedSet& ) {} // no copy constructor
public:
    typedef typename Foam::HashTable<unsigned,T,hash<T> >::const_iterator const_iterator;
    typedef typename Foam::HashTable<unsigned,T,hash<T> >::iterator iterator;

    CountedSet() {}
    void operator+=( T& t )
    {
            if ( this->found(t) )
            {
                this->operator[](t) += 1;
            }
            else
            {
                this->insert(t,1);
            }
    }
    unsigned count( T& t ) const
    {
            return ( this->found(t) ) ? 0 : this->operator[](t);
    }
    unsigned count( const_iterator i ) const
    {
            return *i;
    }
};
typedef CountedSet<Foam::face> FaceSet;

// Builds the list of faces defining the boundary of the domain.
// These are the faces that are referenced by only one cell. 
// This list is used in the lookup of boundary condition patches.
// Second version with reduced memory usage (but might take longer to execute)
void CGNSBoundaryConditions::identifyBoundaryFaces2()
{
    const Foam::cellShapeList& cellsAsShapes = mapper.getCells();
    const Foam::label npoints = mapper.getTotalNodes();
        #if 0
    // build a lookup table to identify the cells adjacent to given point
    forAll( cellsAsShapes, icell )
    {
    const Foam::cellShape& cell = cellsAsShapes[icell];
    const Foam::labelList& cellPoints = cell;
    forAll( cellPoints, icp )
    {
            Foam::label point = cellPoints[icp];
            int bank=0;
            while(1)
            {
                for ( int i=0 ; i<7 ; i++ )
                {
                    if ( pointToCell[bank][point][i] < 0 )
                        pointToCell[bank][point][i] = icell;
                }
            }
    }
    }
        #endif
    // We start from the nodes, lookup adjacent faces and keep only 
    // those that have a single cell as neighbour.
    for( Foam::label point=0 ; point!=npoints ; point++ )
    {
        FaceSet faceset;
    // find adjacent cells
    forAll( cellsAsShapes, icell )
    {
            const Foam::cellShape& cell = cellsAsShapes[icell];
            const Foam::labelList& cellPoints = cell;
            forAll( cellPoints, icp )
            {
                if ( cellPoints[icp] == point )
                {
                    // ok we know 'cell' is a neighbour of 'point'
                    // now lookup its faces and append those that 
                    // include 'point' into faceset
                    const Foam::faceList& faces = cell.faces();
                    forAll( faces, iface )
                    {
                        const Foam::face& face = faces[iface];
                        const Foam::labelList& facepoints = face;
                        bool reject = false;
                        bool accept = false;
                        forAll( facepoints, ifp )
                        {
                            if ( facepoints[ifp] < point )
                            {
                                // already looked up
                                reject = true;
                                break;
                            }
                            else if ( facepoints[ifp] == point )
                            {
                                accept = true;
                                // keep on checking the other points
                            } 
                        }
                        if ( accept && !reject )
                        {
                            faceset += const_cast<Foam::face&>(face);
                        }
                    }
                    break;
                }
            }
    }
    // keep only the faces that have a single reference
    //forAll( faceset, face )
    if ( debug_ )
            Foam::Info << "Creating boundary face list" << Foam::endl;
    for ( FaceSet::iterator i  = faceset.begin() ;
              i != faceset.end() ;
              i++ )
    {
            if ( faceset.count(i)==1 )
            {
                boundary_faces.append( i.key() );
            }
    }
    }
    if ( debug_ )
    Foam::Info << "Shrinking tables" << Foam::endl;
    boundary_faces.shrink();
}

// This method uses a set of point indices to
// identify the set of faces that includes those points
void CGNSBoundaryConditions::pointListToFaceList( std::set<foamNodeIndex>& pointlist, 
Foam::faceList&     l_faces )
{
    list<Foam::face> flist; // temporary dynamic storage (we don't know are many faces we will find)
    forAll ( boundary_faces, iface )
    {
        const Foam::face& f = boundary_faces[iface];
        int nfacevertices = f.size();
        if ( nfacevertices < 3 )
        {
            Foam::Warning << "Degenerated face : " << f << Foam::endl;
            continue; // next face
        }

        // Check if all the points of this face are in the pointlist
        // if so, this face belongs to the boundary patch
        bool f_is_included = true;
        for( int ivert=0 ; ivert<nfacevertices ; ivert++ )
        {
            if ( pointlist.find(f[ivert]) == pointlist.end() )
            {
                f_is_included = false;
                break;
            }
        }

        if ( f_is_included ) flist.push_back( f );
    }

    // Convert from a list<face> to a Foam::faceList
    int i=0;
    l_faces.setSize( flist.size() );
    for ( list<Foam::face>::const_iterator  p  = flist.begin() ;
          p != flist.end() ;
          p++ )
    {
        l_faces[i++] = *p;
    }
}

// This method build the data structures defining the patches where
// boundary conditions are applied. These data structures are:
//    patchFaces - list of the faces defining each patch
//    patchNames - name of each patch
//    patchTypes - type of foam patch (always "patch"  in this context)
//    patchPhysicalTypes - foam's boundary condition type associated with aech patch
void CGNSBoundaryConditions::buildPatches( bool separatePatches )
{
    if ( debug_ )
        Foam::Info << "Identification of boundary conditions patches" << Foam::endl;
    // find out how many patches we will have
    int nbBC = 0;
    for ( std::set<std::string>::const_iterator name  = s_bcname.begin() ; 
          name != s_bcname.end() ;
          name++ )
    {
        if ( separatePatches )
        {
            std::string bcname = *name;
            for ( mm_bcp::const_iterator p  = bcpatches.lower_bound( bcname ) ;
                  p != bcpatches.upper_bound( bcname ) ;
                  p++ )
            {
                nbBC++;
            }
        }
        else
        {
            nbBC++;
        }
    }
    if ( debug_ )
        Foam::Info << nbBC << " patches were identified" << Foam::endl;

    // dimension the tables accordingly
    patchFaces.setSize        ( nbBC );
    patchNames.setSize        ( nbBC );
    patchTypes.setSize        ( nbBC );
    patchPhysicalTypes.setSize( nbBC );

    int bc_counter = 0;

    if ( debug_ )
        Foam::Info << "Extracting patches" << Foam::endl;
    // Loop on all BCs, extract the points and identify associated faces
    for ( std::set<std::string>::iterator name  = s_bcname.begin() ; 
          name != s_bcname.end() ;
          name++ )
    {
        std::string bcname = *name;
        if ( debug_ )
            Foam::Info << "\tExtracting faces of patch named " << bcname << Foam::endl;

        // Find out start/end iterators for this bcname
        mm_bcp::const_iterator p    = bcpatches.lower_bound( bcname );
        mm_bcp::const_iterator pend = bcpatches.upper_bound( bcname );

        // Construct a single list containing all the 
        // points referencing this bcname
        // even if they originally belonged to different zones
        int npatches = 0;
        int nptsbc = 0;
        int nfaces = 0;

        CGNSOO::BCType_t bctype = CGNSOO::BCTypeNull;
        std::set<foamNodeIndex> merged_point_list;
        while( p != pend )
        {
            const patchinfo& patch = (*p).second;
            const std::vector<foamNodeIndex>& pl = patch.pointlist;
            merged_point_list.insert( pl.begin(), pl.end() );
            if ( bctype != CGNSOO::BCTypeNull && bctype != patch.bctype )
            {
                Foam::Warning << "The set of patches defining the BC named " << bcname 
                    << " do not have the same BCType" << Foam::endl;
            }
            bctype = patch.bctype;

            // just for statistics...
            int nptsonpatch = pl.size();
            nptsbc += nptsonpatch;
            npatches++;
            nfaces += patch.nfaces;

            p++;

            // if separatePatches is true, we build one patch for each
            // CGNS BC_t. if not, we merge all the BC_t having the 
            // same name (useful when there are multiple zones).
            if ( p==pend || separatePatches )
            {
                Foam::word foam_phystype = getFoamTypeName( bctype, bcname );
                Foam::word foam_bctype = (foam_phystype=="wallFunctions")?"wall":"patch";

                if ( debug_ )
                    Foam::Info << "\tBC " << bcname << " (physical type : " << foam_phystype << ")"
                        << " contains " << (unsigned int)merged_point_list.size() 
                        << " points on " << npatches 
                        << " distinct patches having a total of " << nptsbc << " points."
                        << Foam::endl;

                //if( foam_phystype != "empty" )
                {
                    // Extract faces
                    Foam::faceList l_faces;
                    pointListToFaceList( merged_point_list, l_faces );

                    if ( nfaces != 0 )
                    {
                        if ( debug_ )
                            Foam::Info << "\tNumber of faces found : " << l_faces.size() 
                                << ", expected " << nfaces
                                << Foam::endl;
                        // nfaces==0 means we don't know how many faces to expect
                        // because the patch was not defined by a PointRange
                    }

                    // Memorize
                    Foam::word fbcname = Foam::string::validate<Foam::word>(Foam::word(bcname));
                    patchFaces        [bc_counter] = l_faces; 
                    patchNames        [bc_counter] = fbcname;
                    patchTypes        [bc_counter] = foam_bctype;
                    patchPhysicalTypes[bc_counter] = foam_phystype;
                    bc_counter++;
                }

                // reset
                merged_point_list.erase( merged_point_list.begin(), 
                merged_point_list.end() );
                npatches = 0;
                nptsbc   = 0;
                nfaces   = 0;
                bctype   = CGNSOO::BCTypeNull;
            }
        }
    }
    patchFaces.setSize        ( bc_counter );
    patchNames.setSize        ( bc_counter );
    patchTypes.setSize        ( bc_counter );
    patchPhysicalTypes.setSize( bc_counter );
}

inline Foam::scalar radians( Foam::scalar degree )
{
  return Foam::constant::mathematical::pi*degree/180.0;
}

Foam::tensor CGNSBoundaryConditions::rotationTransform( 
    Foam::vector            rotAxis,
    Foam::scalar            rotAngle // in degrees
)
{
    const Foam::vector xaxis = Foam::vector(1,0,0);
    const Foam::vector zaxis = Foam::vector(0,0,1);

    rotAxis /= mag(rotAxis); // normalize

    // get any reference vector not aligned with rotAxis 
    Foam::vector ref = xaxis; // just a guess for now
    if ( magSqr( ref ^ rotAxis) < 1.0e-6 )
        ref = zaxis; // change for another one

    // compute a unit basis orthogonal to rotAxis
    Foam::vector baseAxis1 = rotAxis ^ ref;
    baseAxis1 /= mag( baseAxis1 );
    Foam::vector baseAxis2 = baseAxis1 ^ rotAxis;

    // compute the vector resulting from the rotation of base2 around rotAxis
    double alpha = radians(rotAngle);
    Foam::vector baseAxis2Rot = ::cos(alpha) * baseAxis1 + ::sin(alpha) * baseAxis2;

    // We need to check for both positive and negative rotations
        #if OPENFOAM_VERSION<140
    const Foam::tensor rot(Foam::transformationTensor(baseAxis1,baseAxis2Rot));
        #else
    const Foam::tensor rot(Foam::rotationTensor(baseAxis1,baseAxis2Rot));
        #endif

    //Foam::Info << "CGNSBoundaryConditions::rotationTransform: posrotT = "  << rot << Foam::endl;

    return rot;
}

// This takes the "Per" list of faces and applies both the positive and 
// negative transform to identify matching faces from the "Ref" list.
// This routine defines a face by its center location.
// The returned faceList is a permutation of 'cyclicFaceListPer'
// such that its order matches 'cyclicFaceRef' i.e. returned[i] is 
// matched with cyclicFaceListRef[i]
Foam::faceList CGNSBoundaryConditions::computeMatchingCyclicFaces( 
    const Foam::faceList&   cyclicFaceListRef,
    const Foam::faceList&   cyclicFaceListPer,
    const Foam::pointField& points,
    Foam::tensor            posTransform,
    Foam::scalar            cyclicToleranceFactor ) // relative factor
{
    // the inverse of an orthonormal matrix is its transpose
    const Foam::tensor negTransform(posTransform.T());

    // Compute face centers... 
    // Will match the reference face centres from cyclicFaceListRef 
    int nbFaces = cyclicFaceListRef.size();
    Foam::vectorField fListCtrsRef( nbFaces );
    Foam::vectorField fListCtrsPer( nbFaces );
    forAll( cyclicFaceListRef, iface )
    {
        fListCtrsRef[iface] = cyclicFaceListRef[iface].centre(points);
        fListCtrsPer[iface] = cyclicFaceListPer[iface].centre(points);
    }

    // Compute the minimum edge length to normalize tolerance
    Foam::scalar minEdgeLength = Foam::GREAT;
    forAll( cyclicFaceListRef, iface )
    {
        Foam::edgeList eList = cyclicFaceListRef[iface].edges();

        forAll( eList, iedge )
        {
            minEdgeLength = Foam::min(
                minEdgeLength,
                eList[iedge].mag(points)
            );

        }
    }
    // We multiply by the relative tolerance factor
    Foam::scalar maximumMatchingDelta = minEdgeLength*cyclicToleranceFactor;  

    // Apply both transforms to the face centers
    Foam::vectorField fListCtrsPosTr( nbFaces );
    Foam::vectorField fListCtrsNegTr( nbFaces );
    forAll( fListCtrsPer, iface )
    {
        fListCtrsPosTr[iface] = Foam::transform( posTransform, fListCtrsPer[iface] );
        fListCtrsNegTr[iface] = Foam::transform( negTransform, fListCtrsPer[iface] );
    }

    Foam::faceList matchingCyclicFaceListPer( nbFaces );
    /*
    Foam::Info    << "\tcomputeMatchingCyclicFaces: " << Foam::endl
        << "\t\tfListCtrsRef            : " << fListCtrsRef[0]       << Foam::endl
        << "\t\tfListCtrsPer            : " << fListCtrsPer[0]       << Foam::endl
        << "\t\tfListCtrsPosTr          : " << fListCtrsPosTr[0]     << Foam::endl
        << "\t\tfListCtrsNegTr          : " << fListCtrsNegTr[0]     << Foam::endl
        << "\t\tminEdgeLength           : " << minEdgeLength         << Foam::endl
        << "\t\tcyclicMatchFaceTolFactor: " << cyclicToleranceFactor << Foam::endl
        << "\t\tmaximumMatchingDelta    : " << maximumMatchingDelta  << Foam::endl;
        */
    // Now, let's find the correct transform and the correct face association
    Foam::scalar smallestDeltaFound = std::numeric_limits<Foam::scalar>::max();
    Foam::scalar curDelta = 0;
    int nbMatch = 0;

    enum { POSITIVE_TRANSFORM, NEGATIVE_TRANSFORM, UNKNOWN } transform_dir = UNKNOWN;
    forAll( fListCtrsRef, iface )
    {
        bool foundMatch = false;

        // Keep note of the minimum delta found.
        // That way, if we have no match for a point, we can at
        // least have an idea about the minimum delta found for this point.
        if ( transform_dir != NEGATIVE_TRANSFORM )
        {
            forAll( fListCtrsPosTr, jface )
            {
                curDelta = mag( fListCtrsRef[iface] - fListCtrsPosTr[jface]);

                if( curDelta < maximumMatchingDelta)
                {
                    matchingCyclicFaceListPer[iface] = cyclicFaceListPer[jface];
                    foundMatch = true;
                    nbMatch++;
                    break;
                }
                else
                {
                    if ( curDelta < smallestDeltaFound )
                        smallestDeltaFound = curDelta;
                }
            }
        }
        if (foundMatch) 
        {
            transform_dir = POSITIVE_TRANSFORM;
            continue; // pick next face
        }

        // check using the other direction
        if ( transform_dir != POSITIVE_TRANSFORM )
        {
            forAll( fListCtrsNegTr, jface )
            {
                curDelta = mag(fListCtrsRef[iface] - fListCtrsNegTr[jface]);

                if( curDelta < maximumMatchingDelta )
                { 
                    matchingCyclicFaceListPer[iface] = cyclicFaceListPer[jface];
                    foundMatch = true;
                    nbMatch++;
                    break;
                }
                else
                {
                    if ( curDelta < smallestDeltaFound )
                        smallestDeltaFound = curDelta;
                }
            }
        }
        if (foundMatch) 
        {
            transform_dir = NEGATIVE_TRANSFORM;
            continue; // pick next face
        }
        else
        {
            if ( debug_ )
                Foam::Info     << "No match found for : " << fListCtrsRef[iface]
                    << " after trying both forward and backward tranformation" << Foam::endl
                    << "Smallest delta found : " << smallestDeltaFound << Foam::endl;
            break;
        }
    }

    if( nbMatch != nbFaces )
    {
        Foam::FatalError << "Cannot find enough cyclic matches." << Foam::endl
            << "\tFound " << nbMatch << ", expecting " << nbFaces << Foam::endl
            << "\tCyclic match tolerance = " << maximumMatchingDelta << Foam::endl
            << exit(Foam::FatalError);
    }

    // Now, lets make sure the first vertex of the matching faces are matching
    // this seems to be a requirement for Foam::globalPoints and Foam::cyclicPolyPatch

    forAll( cyclicFaceListRef, iface )
    {
        const Foam::labelList& vListRef = cyclicFaceListRef[iface];
        Foam::labelList& vListPer = matchingCyclicFaceListPer[iface];

        Foam::point pRef = points[vListRef[0]];
        Foam::label startVert = -1;
        forAll( vListPer, ivert )
        {
            Foam::point pPerTransform = Foam::transform(
                ( transform_dir==POSITIVE_TRANSFORM ) ? 
                posTransform : negTransform,
                points[vListPer[ivert]] );
            if ( mag( pPerTransform - pRef ) < maximumMatchingDelta )
            {
                startVert = ivert;
                break;
            }
        }
        if ( startVert == -1 )
        {
            FatalErrorIn("computeMatchingCyclicFaces")
                << "Cyclic faces do not match - "
                "this is probably an internal error" 
                    << Foam::endl
                    << abort(Foam::FatalError);
        }
        else if ( startVert != 0 )
        {
            Foam::label nvert = vListPer.size();
            Foam::labelList vshift( nvert );
            forAll( vListPer, iv )
                vshift[iv] = vListPer[ (iv+startVert)%nvert ];
            vListPer = vshift;
        }
    }

    return matchingCyclicFaceListPer;
}

// This is the public interface for computing matching boundary conditions
// the criteria used to decide if two patches are matching is the location
// of their constituent nodes.
int CGNSBoundaryConditions::computeCyclicBC( 
    const std::string&       firstPatchName, 
    const std::string&       secondPatchName,
    const Foam::vector& rotAxis,
    Foam::scalar        rotAngle,
    Foam::scalar        cyclicToleranceFactor )
{
    // internally, patches are indexed with FOAM::word
    std::string fpname1 = FoamReplaceChar(firstPatchName ,' ','_');
    std::string fpname2 = FoamReplaceChar(secondPatchName,' ','_');

    // Find patch indices for the two matching BCs
    int ipatch1 = -1;
    int ipatch2 = -1;
    int counter = 0;
    forAll( patchNames, patchNamesI )
    {
        if( patchNames[patchNamesI] == fpname1 )
        {
            ipatch1 = counter;
        }
        else if( patchNames[patchNamesI] == fpname2 )
        {
            ipatch2 = counter;
        }
        counter++;
    }

    // Validation
    if ( ipatch1 == -1 ) return -1;
    if ( ipatch2 == -1 ) return -2;
    if ( patchFaces[ipatch1].size() != patchFaces[ipatch2].size() )
    {
        Foam::FatalError << "cyclic patch size mismatch." << Foam::endl
            << "\tBC #1: " << firstPatchName  
            << " : size=" << patchFaces[ipatch1].size() << Foam::endl
            << "\tBC #2: " << secondPatchName 
            << " : size=" << patchFaces[ipatch2].size() << Foam::endl
            << exit(Foam::FatalError);
    }

    // We have two patch definitions that represent the two sides of
    // the periodic interface. Our aim is to combine these patches into
    // single one

    // Compute correct faces association
    Foam::tensor transformation = rotationTransform( rotAxis, rotAngle );
    const Foam::pointField newPoints = mapper.getPoints();
    Foam::faceList matchingFaces = computeMatchingCyclicFaces( 
        patchFaces[ipatch1],
        patchFaces[ipatch2],
        newPoints,
        transformation,
        cyclicToleranceFactor );

    if ( debug_ )
        Foam::Info << "\tPatches " << firstPatchName << " and " << secondPatchName 
            << " were correctly identified as cyclic."
            << Foam::endl;

    // Copy the faces of theses two patches into a single faceList
    int offset = patchFaces[ipatch1].size();
    Foam::faceList flist_coupled(2*offset);
    forAll( patchFaces[ipatch1], iface )
    {
        flist_coupled[iface]        = patchFaces[ipatch1][iface]; // First part: BC#1 
        flist_coupled[offset+iface] = matchingFaces[iface]; // Complement with BC#2 
    }

    // Update the various tables defining the patches
    // Remove patch entries occupied by ipatch2
    // The new cyclic BC goes into the slot previously used by patch1
    int iLastPatch = patchFaces.size()-1;
    if ( ipatch2 != iLastPatch )
    {
        // move last patch into the ipatch2 slot
        patchFaces[ipatch2]         = patchFaces[iLastPatch];
        patchNames[ipatch2]         = patchNames[iLastPatch];
        patchTypes[ipatch2]         = patchTypes[iLastPatch];
        patchPhysicalTypes[ipatch2] = patchPhysicalTypes[iLastPatch]; 
    }

    // Memorize cyclic patches into slot ipatch1
    patchFaces[ipatch1]         = flist_coupled;
    patchNames[ipatch1]         = fpname1 + "_and_" + fpname2;
    patchTypes[ipatch1]         = "cyclic";
    patchPhysicalTypes[ipatch1] = "cyclic";

    // Resize - removes the last slot which is now useless
    patchFaces        .setSize(iLastPatch);
    patchNames        .setSize(iLastPatch);
    patchTypes        .setSize(iLastPatch);
    patchPhysicalTypes.setSize(iLastPatch);

    return 0;
}

Foam::polyMesh* CGNSBoundaryConditions::buildFoamMesh( const Foam::Time& t ) const
{
    //int n = patchFaces.size();
    if ( debug_ )
    {
        Foam::Info << "Patches:\n";
        forAll( patchNames, i )
        {
            Foam::Info << "\tPatch " << i << " named " << patchNames[i] 
                << ", type = " << patchTypes[i] << "/" << patchPhysicalTypes[i] 
                << Foam::endl;
        }
    }
    Foam::word defaultPatchName = "defaultWall";
    Foam::word defaultPatchType = Foam::wallPolyPatch::typeName;
    auto points = mapper.getPoints();
    return new Foam::polyMesh(
        Foam::IOobject
        (
            Foam::polyMesh::defaultRegion,
            t.constant(),
            t
        ),
        std::move(points),
        mapper.getCells(),
        patchFaces,
        patchNames,
        patchTypes,
        defaultPatchName,
        defaultPatchType,
        patchPhysicalTypes
    );
}

// This method uses the name and the CGNS boundary condition types of a patch
// to return a suitable foam boundary condition type
// Essentially, the conversions follow these lines:
//     *** CGNS ***          *** FOAM ***
//    BCWall           ---> wallFunctions
//    AxisymetricWedge ---> wedge
//    BCInflow         ---> inlet
//    BCOutflow        ---> outlet
//    BCTypeNull 
//    UserDefined      ---> search the patch name for a one of the token { in,out,wall}
//    everything else  ---> empty
/*static*/ Foam::word CGNSBoundaryConditions::getFoamTypeName( CGNSOO::BCType_t bctype, const std::string& bcname )   
{
    Foam::word s_retValue = "empty";  // default value
    Foam::word s_wall     = "wallFunctions";  // "wall" // old id used in version 1.1
    switch( bctype )
    {
    case CGNSOO::BCWall                  :
            s_retValue = s_wall; 
            break;
    case CGNSOO::BCTypeNull              :
    case CGNSOO::BCTypeUserDefined       :
            {
        std::string l_bcname = toLower(bcname);  // search in lowercase

        if( l_bcname.find("wall")    != std::string::npos) s_retValue = s_wall; 
        else if(l_bcname.find("in")  != std::string::npos) s_retValue = "inlet";
        else if(l_bcname.find("out") != std::string::npos) s_retValue = "outlet";
            }
            break;
    case CGNSOO::BCAxisymmetricWedge     :
            s_retValue = "wedge";
            break;
    case CGNSOO::BCInflow                :
            s_retValue = "inlet";
            break;
    case CGNSOO::BCOutflow               :
            s_retValue = "outlet";
            break;
    case CGNSOO::BCSymmetryPlane         :
            s_retValue = "symmetryPlane";
            break;
            /*
                case BCOutflowSubsonic       :
                case BCOutflowSupersonic     :
                case BCSymmetryPolar         :
                case BCTunnelInflow          :
                case BCTunnelOutflow         :
                case BCDegenerateLine        :
                case BCDegeneratePoint       :
                case BCDirichlet             :
                case BCExtrapolate           :
                case BCFarfield              :
                case BCGeneral               :
                case BCInflowSubsonic        :
                case BCInflowSupersonic      :
                case BCNeumann               :
                case BCWallInviscid          :
                case BCWallViscous           :
                case BCWallViscousHeatFlux   :
                case BCWallViscousIsothermal :
                case FamilySpecified         :
                */
    default:
            break;
    }

    return s_retValue;
}

/*static*/ bool CGNSBoundaryConditions::validateUserDefinedBCName( const std::string& bcname )
{
    // convert to lowercase
    std::string l_bcname = toLower(bcname);

    return    l_bcname.find("wall")    != std::string::npos ||
        l_bcname.find("inflow")  != std::string::npos ||
        l_bcname.find("outflow") != std::string::npos ||
        l_bcname.find("inlet")   != std::string::npos ||
        l_bcname.find("outlet")  != std::string::npos
        ;
}
