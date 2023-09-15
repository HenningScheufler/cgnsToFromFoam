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

#ifndef ConnectivityMapper_H
#define ConnectivityMapper_H

namespace Foam
{
class cellModel;
}


#include "point.H"
#include "cellShapeList.H"
#include "pointField.H"
#include "polyMesh.H"
// #include "Xfer.H"

#include "libcgnsoo3/cgnsoo.H"

#include <vector>
#include <string>
#include <list>

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
typedef Foam::label foamNodeIndex;

//------------------------------------------------------------------------------
// This class handles node and cell numbering operations on multizone meshes
// Each zone can be defined separately with its own nodes and connectivity.
// The merge() method builds a unique global numbering and corresponding
// mapping tables from the old numbering to the new one.
//------------------------------------------------------------------------------
class ConnectivityMapper
{
private:
    //!< structure to keep info about a zone of 3D cells
    struct zone_data_t
    {
        bool structured; // is this zone structured
        int  ni, nj, nk; // for a structured zone, number of nodes in each directions
        int  n_cells;    // number of cells in this zone
        int  n_nodes;    // number of nodes in this zone
        int  offset;     // index of the first node of this zone in the global node list
        list<Foam::point> pts; // nodes of this zone
    };
public:
    //!< structure to keep the information about a 2D region
    struct zone2d_data_t
    {
        int           zindex;
        std::string        name;
        std::vector<int>   connec;
        CGNSOO::ElementType_t etype;
		
        zone2d_data_t( int z, const std::string& n, const std::vector<int>& c, CGNSOO::ElementType_t t ) : 
            zindex(z), name(n), connec(c), etype(t) 
            {}
    };
private:
    std::vector<zone_data_t>   zdata;
    int                   g_n_nodes;
    int                   g_n_cells;
    Foam::labelList       oldToNew;
    list<Foam::cellShape> tmpcells;
    Foam::cellShapeList   cellShapes;
    Foam::pointField      newPoints;
    Foam::polyMesh*       pShapeMesh;
    bool                  merged;
    std::vector<zone2d_data_t> uns_z2d;
	
public:
    ConnectivityMapper() : 
        g_n_nodes(0), 
        g_n_cells(0), 
        pShapeMesh(NULL), 
        merged(false) {}
    ~ConnectivityMapper()
	{
            delete pShapeMesh;
	}
	
    // methods for building the datastructures
    void addStructuredZone( int nni, int nnj, int nnk, const list<Foam::point>& p );
    void addUnstructuredZone( int nn, int nc, const list<Foam::point>& p,
    const std::vector<std::string>& sectionnames, 
    const std::vector<CGNSOO::ElementType_t>& etype, 
    const std::vector< std::vector<int> >& connectivity );
	
    // method to compute matching nodes
    void merge( double mergeTolerance = 1.0e-1 );
	
    // methods for inquiries
    int  n_localToGlobal( int z, int localnodeid ) const;
    int  c_localToGlobal( int z, int localcellid ) const;
    int  getTotalNodes() const { return g_n_nodes; }
    int  getTotalCells() const { return g_n_cells; }
    int  getLocalNodes( int zid ) const { return zdata[zid].n_nodes; }
    int  getLocalCells( int zid ) const { return zdata[zid].n_cells; }
    static int getLocalNodeIndexFromStructuredCGNSTriplet( int ni, int nj, int nk, int i, int j, int k );
    void getStructuredZoneDimensions( int zid, int& ni, int& nj, int& nk ) const
	{
            if ( !zdata[zid].structured )
                throw std::logic_error( "zone is not structured" );
            ni = zdata[zid].ni;
            nj = zdata[zid].nj;
            nk = zdata[zid].nk;
	}
    std::vector<zone2d_data_t>      getSectionsInfo() const { return uns_z2d; }
    const Foam::pointField&         getPoints() const { return newPoints; }

    Foam::cellShapeList&       getCells();
    const Foam::cellShapeList& getCells() const;
	
private:	
    void addStructuredConnectivity( int zoneIndex, int nni, int nnj, int nnk );
    void addUnstructuredConnectivity( int zoneIndex, 
    const std::vector<CGNSOO::ElementType_t>& etype, 
    const std::vector< std::vector<int> >& connectivity );
					  
    Foam::cellShape buildOneCellFromCGNSConnectivity( 
        std::vector<int>::const_iterator  connectivity_ptr,
        int                          numberOfNodesPerElement,
        const Foam::cellModel&       cell_model,
        int                          zoneIndex );
		
    void extractCellsFromMixedCgnsElements( 
        const std::vector<int>&      connectivity,
        int                     zoneIndex,
        int                     eSectionIndex,
        list<Foam::cellShape>&  l_cells );
		
    void extractCellsFromUnstructuredCgnsElements( 
        CGNSOO::ElementType_t   elementType,
        const std::vector<int>&      connectivity,
        int                     zoneIndex,
        int                     eSectionIndex,
        list<Foam::cellShape>&  l_cells );
		
    double computeMinEdgeLength( const Foam::pointField& points ) const;
};

#endif
