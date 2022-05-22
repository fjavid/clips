/**
 *  PicoParaview v1.0.4
 *
 *  A small library to write files for visualization in Paraview.
 *
 *  Core functions related to the .vtu, .vts and .pvd file formats are in
 *  'ppv::fmt' namespace. Convenience functions for common use cases are in the
 *  'ppv' namespace. Demo examples are in the 'ppv::demo' namespace.
 *
 *  To use PicoParaview in your project you need to:
 *  1. Include picoparaview.h as a regular header wherever needed.
 *  2. In one instance #define PICOPARAVIEW_IMPLEMENTATION` before the include.
 */

#ifndef INCLUDE_PICOPARAVIEW_H
#define INCLUDE_PICOPARAVIEW_H

#include <algorithm>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace ppv
{
    namespace fmt
    {
        /**
         *  Enum type of VTK elements. We use this enum to distinguish elements rather
         *  than node numbers since elements of different types can have the same number
         *  of nodes (e.g. triangles and quadratic line elements).
         */
        typedef enum e_VtkType {
            UNDEFINED = -1,
            POINT     = 1,
            // linear Lagrange elements
            LINE = 3,
            TRI  = 5,
            QUAD = 9,
            TET  = 10,
            HEX  = 12,
            // quadratic Lagrange elements
            LINE2 = 21,
            TRI2  = 22,
            QUAD2 = 23,
            TET2  = 24,
            HEX2  = 25,
            // cubic Lagrange elements
            LINE3 = 35
        } VtkType;

        /**
         *  Helper function to give number of nodes of given vtk element type
         *
         *  @param vtkType VTK element type
         *  @return        Number of nodes in element
         */
        std::size_t VtkTypeToNumNodes( VtkType vtkType );

        /**
         *  Start writing VTU file, this is the first call
         *
         *  @param  vtu          VTU file output stream
         *  @param  nNumNodes    total number of nodes to write
         *  @param  nNumElements total number of elements (of all types) to write
         */
        void OpenVTU( std::ofstream& vtu, int nNumNodes, int nNumElements );

        /**
         *  Finish writing VTU file, this is the last call
         *
         *  @param  vtu VTU file output stream
         */
        void CloseVTU( std::ofstream& vtu );

        /**
         *  Start writing VTS file, this is the first call.
         *
         *  @param  vts VTU file output stream
         *  @param  vx  number of elements (voxels, *not* points) in x-direction
         *  @param  vy  number of elements (voxels, *not* points) in y-direction
         *  @param  vz  number of elements (voxels, *not* points) in z-direction
         */
        void OpenVTS( std::ofstream& vts, int vx, int vy, int vz );

        /**
         *  Finish writing VTS file, this is the last call
         *
         *  @param  vts VTS file output stream
         */
        void CloseVTS( std::ofstream& vts );

        /**
         *  Start writing VTI file. This is the first call.
         *
         *  @param vti VTI file output stream
         *  @param ex  number of elements in x-direction
         *  @param ey  number of elements in y-direction
         *  @param ez  number of elements in z-direction
         *  @param sx  element size in x-direction
         *  @param sy  element size in y-direction
         *  @param sz  element size in z-direction
         *  @param ox  origin x-component
         *  @param oy  origin y-component
         *  @param oz  origin z-component
         */
        template <typename Real>
        void OpenVTI( std::ofstream& vti, int ex, int ey, int ez,
                      Real sx, Real sy, Real sz,
                      Real ox, Real oy, Real oz) {
          vti << "<?xml version=\"1.0\"?>\n"
              << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
              << "<ImageData WholeExtent=\"0 " << ex << " 0 " << ey << " 0 " << ez << '\"'
              << " Origin=\"" << ox << ' ' << oy << ' ' << oz << '\"'
              << " Spacing=\"" << sx << ' ' << sy << ' ' << sz << "\">\n";
        }

        /**
         *  Finish writing VTI file. This is the last call.
         *
         *  @param  vti VTI file output stream
         */
        void CloseVTI( std::ofstream& vti );

        /**
         *  Open piece of VTI file
         *
         *  @param vti VTI file output stream
         *  @param ex  number of elements in x-direction
         *  @param ey  number of elements in y-direction
         *  @param ez  number of elements in z-direction
         */
        void OpenPiece(std::ostream &vti, int ex, int ey, int ez);

        /**
         *  Close piece of VTI file
         */
        void ClosePiece(std::ostream &out);

        /**
         *  Begin writing mesh node coordinates
         *
         *  @param  vtu VTU file output stream
         */
        void OpenNodeCoordinates( std::ofstream& vtu );

        /**
         *  Finish writing mesh node coordinates
         *
         *  @param  vtu VTU file output stream
         */
        void CloseNodeCoordinates( std::ofstream& vtu );

        /**
         *  Begin writing mesh elements.  Called only once even for mixed element
         *  meshes.
         *
         *  @param  vtu VTU file output stream
         */
        void OpenElementConnectivity( std::ofstream& vtu );

        /**
         *  Finish writing mesh coordinates.  Must also remember to write offset and
         *  element types and offsets prior to making this call.
         *
         *  @param  vtu VTU file output stream
         */
        void CloseElementConnectivity( std::ofstream& vtu );

        /**
         *  Begin writing nodal data arrays
         *
         *  @param  vtu VTU file output stream
         */
        void OpenNodeData( std::ofstream& vtu );

        /**
         *  Finish writing nodal data arrays
         *
         *  @param  vtu VTU file output stream
         */
        void CloseNodeData( std::ofstream& vtu );

        /**
         *  Begin writing element data arrays
         *
         *  @param  vtu VTU file output stream
         */
        void OpenElementData( std::ofstream& vtu );

        /**
         *  Finish writing element data arrays
         *
         *  @param  vtu VTU file output stream
         */
        void CloseElementData( std::ofstream& vtu );

        /**
         *  Begin writing VTU data array. IMPORTANT: if writing VTU element
         *  connectivity must set nNumComp as 1!
         *
         *  @param vtu      VTU file output stream
         *  @param compName "Float64" for floating-point types, "Int32" otherwise
         *  @param nNumComp Number of components in each DataArray record
         *  @param dataName Name of data as it appears for user in Paraview
         */
        void OpenDataArray( std::ofstream& vtu,
                            std::string    compName,
                            int            nNumComp,
                            std::string    dataName );

        /**
         *  Finish writing VTU data array
         *
         *  @param vtu VTU file output stream
         */
        void CloseDataArray( std::ofstream& vtu );

        /**
         *  Append to VTU offset container
         *
         *  @param vctOffsets   Container of element offsets
         *  @param nNumElements Number of components in each DataArray record
         *  @param vtkType      Enum type of vtk element
         */
        void AddToOffsetVector( std::vector<int>& vctOffsets, int nNumElements, VtkType vtkType );

        /**
         *  Append to VTU element type container
         *
         *  @param vctVtkTypes  Container of vtk element types
         *  @param nNumElements Number of components in each DataArray record
         *  @param vtkType      Enum type of vtk element
         */
        void
        AddToVtkTypesVector( std::vector<int>& vctVtkTypes, int nNumElements, VtkType vtkType );

        /**
         *  Open a PVD file.  PVD files are a convenient way to visualise collections
         *  of vtu files written separately
         *
         *  @param pvd PVD file output stream
         */
        void OpenPVD( std::ostream& pvd );

        /**
         *  Add a VTU file to a PVD collection
         *
         *  @param pvd      PVD file output stream
         *  @param filename Name of file to add to collection, usually a VTU file
         *  @param time     All VTU files to be visualised simultaneously should be
         *                  given the same time
         */
        void RecordPVDDataSet( std::ostream& pvd, std::string filename, float time = 0.0f );

        /**
         *  Close a PVD file
         *
         *  @param pvd PVD file output stream
         */
        void ClosePVD( std::ostream& pvd );
    }  // namespace fmt

    /**
     *  Write geometry of mesh.  Not a stand-alone function, must at least
     *  surround with OpenVTU() and CloseVTU() calls
     *
     *  @param  vtu      VTU file output stream
     *  @param  nodes    Nodes to be written
     *  @param  elements Elements to be written
     *  @param  vtkType  Vtk enum type to be written
     */
    template <typename Real>
    void WriteGeometryVTU( std::ofstream &vtu,
                           const std::vector<std::vector<Real>>& nodes,
                           const std::vector<std::vector<int>>& elements,
                           int vtkType )
    {
        fmt::OpenNodeCoordinates( vtu );
        fmt::OpenDataArray( vtu, "Float64", 3, "coordinates" );
        std::for_each(nodes.begin(), nodes.end(),
                      [&vtu](const std::vector<Real> &node) {
                        vtu << node[0] << ' ' << node[1] << ' ' << node[2] << '\n';
                      });
        fmt::CloseDataArray(vtu);
        fmt::CloseNodeCoordinates( vtu );

        auto writeInt = [&vtu]( int i ) { vtu << i << '\n'; };

        fmt::OpenElementConnectivity( vtu );
        // NOTE(kecmanm): nNumComp == 1 is intentional and expected by Paraview
        fmt::OpenDataArray( vtu, "Int32", 1, "connectivity" );
        std::for_each(elements.begin(), elements.end(),
                      [&vtu](const std::vector<int> &e) {
                        for (size_t ei = 0; ei < e.size() - 1; ei++) {
                          vtu << e[ei] << ' ';
                        }
                        vtu << e[e.size() - 1] << '\n';
                      });
        fmt::CloseDataArray(vtu);
        std::vector<int> vctOffsets;
        // TODO(kecmanm): Compiler issue prevented using VtkType
        fmt::AddToOffsetVector( vctOffsets, (int)elements.size(), (fmt::VtkType)vtkType );
        fmt::OpenDataArray( vtu, "Int32", 1, "offsets" );
        std::for_each( vctOffsets.begin(), vctOffsets.end(), writeInt );
        fmt::CloseDataArray( vtu );
        std::vector<int> vctVtkTypes;
        fmt::AddToVtkTypesVector( vctVtkTypes, (int)elements.size(), (fmt::VtkType)vtkType );
        fmt::OpenDataArray( vtu, "Int32", 1, "types" );
        std::for_each( vctVtkTypes.begin(), vctVtkTypes.end(), writeInt );
        fmt::CloseDataArray( vtu );
        fmt::CloseElementConnectivity( vtu );
    }

    namespace demo {
      void WriteExampleVTU();
      void WriteExampleVTS();
      void WriteExampleVTI();
      void WriteExamplePVD();
    }

    namespace experimental {

    }
}  // namespace ppv
#endif // INCLUDE_PICOPARAVIEW_H

#ifdef PICOPARAVIEW_IMPLEMENTATION
#include <algorithm>
#include <cassert>
#include <fstream>
#include <string>

void ppv::fmt::OpenVTU( std::ofstream& vtu, int nNumNodes, int nNumElements )
{
    vtu << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << nNumNodes << "\" NumberOfCells=\"" << nNumElements << "\">\n";
}

void ppv::fmt::CloseVTU( std::ofstream& vtu )
{
    vtu << "</Piece>\n"
        << "</UnstructuredGrid>\n"
        << "</VTKFile>\n";
}

void ppv::fmt::OpenVTS( std::ofstream& vts, int vx, int vy, int vz )
{
    vts << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"StructuredGrid\" version="
        << "\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<StructuredGrid WholeExtent=\"" << "0 " <<  vx << " 0 " << vy << " 0 " << vz << "\">\n"
        << "<Piece Extent=\"" << "0 " <<  vx << " 0 " << vy << " 0 " << vz << "\">\n";
}

void ppv::fmt::CloseVTS( std::ofstream& vts )
{
    vts << "</Piece>\n"
        << "</StructuredGrid>\n"
        << "</VTKFile>\n";
}

void ppv::fmt::CloseVTI( std::ofstream& vti )
{
    vti << "</ImageData>\n"
        << "</VTKFile>\n";
}

void ppv::fmt::OpenPiece(std::ostream &vti, int ex, int ey, int ez) {
  vti << "<Piece Extent=\"0 " << ex << " 0 " << ey << " 0 " << ez << "\">\n";
}

void ppv::fmt::ClosePiece(std::ostream &out) {
  out << "</Piece>\n";
}

void ppv::fmt::OpenNodeCoordinates( std::ofstream& vtu )
{
    vtu << "<Points>\n";
}

void ppv::fmt::CloseNodeCoordinates( std::ofstream& vtu )
{
    vtu << "</Points>\n";
}

void ppv::fmt::OpenElementConnectivity( std::ofstream& vtu )
{
    vtu << "<Cells>\n";
}

void ppv::fmt::CloseElementConnectivity( std::ofstream& vtu )
{
    vtu << "</Cells>\n";
}

void ppv::fmt::OpenNodeData( std::ofstream& vtu )
{
    vtu << "<PointData>\n";
}

void ppv::fmt::CloseNodeData( std::ofstream& vtu )
{
    vtu << "</PointData>\n";
}

void ppv::fmt::OpenElementData( std::ofstream& vtu )
{
    vtu << "<CellData>\n";
}

void ppv::fmt::CloseElementData( std::ofstream& vtu )
{
    vtu << "</CellData>\n";
}

void ppv::fmt::OpenDataArray( std::ofstream& vtu,
                                  std::string    compName,
                                  int            nNumComp,
                                  std::string    dataName )
{
    vtu << "<DataArray type=\"" << compName << "\" NumberOfComponents=\"" << nNumComp
        << "\" Name=\"" << dataName << "\" format=\"ascii\">\n";
}

void ppv::fmt::CloseDataArray( std::ofstream& vtu )
{
    vtu << "</DataArray>\n";
}

std::size_t ppv::fmt::VtkTypeToNumNodes( VtkType vtkType )
{
    std::size_t numNodes = 0;
    switch ( vtkType )
    {
    case POINT:
        numNodes = 1;
        break;
    case LINE:
        numNodes = 2;
        break;
    case TRI:
        numNodes = 3;
        break;
    case QUAD:
        numNodes = 4;
        break;
    case TET:
        numNodes = 4;
        break;
    case HEX:
        numNodes = 8;
        break;
    case LINE2:
        numNodes = 3;
        break;
    case TRI2:
        numNodes = 6;
        break;
    case QUAD2:
        numNodes = 9;
        break;
    case TET2:
        numNodes = 10;
        break;
    case HEX2:
        numNodes = 27;
        break;
    case LINE3:
        numNodes = 4;
        break;
    case UNDEFINED:
        numNodes = 0;
        break;
    }
    return numNodes;
}

void ppv::fmt::AddToOffsetVector( std::vector<int>& vctOffsets,
                                      int               nNumElements,
                                      VtkType           vtkType )
{
    int currOffset = vctOffsets.empty() ? 0 : vctOffsets[vctOffsets.size()-1];
    size_t oldSize = vctOffsets.size();
    size_t newSize = oldSize + nNumElements;
    vctOffsets.resize( newSize );
    size_t numNodes = VtkTypeToNumNodes( vtkType );
    for ( size_t e = oldSize; e < newSize; e++ )
    {
        vctOffsets[e] = currOffset + (int)((e - oldSize + 1) * numNodes);
    }
}

void ppv::fmt::AddToVtkTypesVector( std::vector<int>& vctVtkTypes,
                                        int               nNumElements,
                                        VtkType           vtkType )
{
    size_t oldSize = vctVtkTypes.size();
    size_t newSize = oldSize + nNumElements;
    vctVtkTypes.resize( newSize );
    for ( size_t e = oldSize; e < newSize; e++ )
    {
        vctVtkTypes[e] = vtkType;
    }
}

void ppv::fmt::OpenPVD( std::ostream& pvd )
{
    pvd << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<Collection>\n";
};

void ppv::fmt::RecordPVDDataSet( std::ostream& pvd, std::string filename, float time )
{
    pvd << "<DataSet timestep=\"" << std::to_string( time ) << "\" file=\"" << filename << "\"/>\n";
};

void ppv::fmt::ClosePVD( std::ostream& pvd )
{
    pvd << "</Collection>\n" << "</VTKFile>\n";
};

void ppv::demo::WriteExampleVTU()
{
    const int nn = 8; // number of nodes
    const int ne = 6; // number of elements
    const float nodes[nn][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    const int elements[ne][4] = {{0,1,2,4},{5,2,1,4},{2,5,6,4},{1,3,2,5},{2,3,6,5},{7,6,3,5}};

    std::ofstream vtu("picoparaview-example.vtu");
    using namespace ppv::fmt;

    auto writeInt = [&vtu]( int i ) { vtu << i << '\n'; };

    // open vtu file (required)
    OpenVTU(vtu, nn, ne);

    { // write coordinates (required)
        OpenNodeCoordinates(vtu);
        OpenDataArray(vtu, "Float64", 3, "coordinates");
        for (int n = 0; n < nn; n++) {
            vtu << nodes[n][0] << ' ' << nodes[n][1] << ' ' << nodes[n][2] << '\n';
        }
        CloseDataArray(vtu);
        CloseNodeCoordinates(vtu);
    }

    { // write connectivity (required)
        OpenElementConnectivity( vtu );
        OpenDataArray( vtu, "Int32", 1, "connectivity" ); // CARE! must be 1!
        for (int e = 0; e < ne; e++) {
            vtu << elements[e][0] << ' ' << elements[e][1] << ' ' << elements[e][2] << ' ' << elements[e][3] << '\n';
        }
        CloseDataArray( vtu );
        // vtu file format requires additional 'offsets' and 'types' data (these allow for different element types in same mesh)
        OpenDataArray( vtu, "Int32", 1, "offsets" );
        std::vector<int> vctOffsets;
        AddToOffsetVector( vctOffsets, ne, TET );
        std::for_each( vctOffsets.begin(), vctOffsets.end(), writeInt );
        CloseDataArray( vtu );
        OpenDataArray( vtu, "Int32", 1, "types" );
        std::vector<int> vctVtkTypes;
        AddToVtkTypesVector( vctVtkTypes, ne, TET );
        std::for_each( vctVtkTypes.begin(), vctVtkTypes.end(), writeInt );
        CloseDataArray( vtu );
        CloseElementConnectivity( vtu );
    }

    { // write node data (optional)
        OpenNodeData(vtu);
        OpenDataArray(vtu, "Int32", 1, "node integer scalar");
        for (int n = 0; n < nn; n++) { vtu << n << '\n'; }
        CloseDataArray(vtu);
        OpenDataArray(vtu, "Int32", 3, "node bool vector");
        for (int n = 0; n < nn; n++) { vtu << (bool)(n%2) << ' ' << (bool)(n%3) << ' ' << (bool)(n%5) << '\n'; }
        CloseDataArray(vtu);
        CloseNodeData(vtu);
    }

    { // write element data (optional)
        OpenElementData(vtu);
        OpenDataArray(vtu, "Float64", 6, "element float tensor");
        for (int e = 0; e < ne; e++) {
            vtu << 0*e << ' ' << 1*e << ' ' << 2*e << ' ' << 3*e << ' ' << 4*e << ' ' << 5*e << '\n';
        }
        CloseDataArray(vtu);
        CloseElementData(vtu);
    }

    // close vtu file (required)
    CloseVTU(vtu);
    vtu.close();
    return;
}

void ppv::demo::WriteExampleVTS()
{
   std::ofstream vts("picoparaview-example.vts");
    using namespace ppv::fmt;

    int vpd = 2; // elements (cells/voxels) per coordinate direction
    int npd = vpd+1; // nodes (vertices/coordinates) per coordinate direction

    // open vts file (required)
    OpenVTS(vts, vpd, vpd, vpd);

    { // write coordinates of cartesian grid, connectivity is implied (required)
        OpenNodeCoordinates(vts);
        OpenDataArray(vts, "Float64", 3, "coordinates");
        float h = 0.5;
        for (int k = 0; k < npd; k++) {
            for (int j = 0; j < npd; j++) {
                for (int i = 0; i < npd; i++) {
                    vts << i*h << ' ' << j*h << ' ' << k*h << '\n';
                }
            }
        }
        CloseDataArray(vts);
        CloseNodeCoordinates(vts);
    }

    { // write node data (optional)
        OpenNodeData(vts);
        OpenDataArray(vts, "Int32", 1, "node integer scalar");
        for (int k = 0; k < npd; k++) {
            for (int j = 0; j < npd; j++) {
                for (int i = 0; i < npd; i++) {
                    vts << i+j+k << '\n';
                }
            }
        }
        CloseDataArray(vts);
        OpenDataArray(vts, "Float64", 3, "node float vector");
        for (int k = 0; k < npd; k++) {
            for (int j = 0; j < npd; j++) {
                for (int i = 0; i < npd; i++) {
                    vts << (float)(i*i) << ' ';
                    vts << (float)(i*j) << ' ';
                    vts << (float)(i*k) << '\n';
                }
            }
        }
        CloseDataArray(vts);
        CloseNodeData(vts);
    }

    { // write element data (optional)
        OpenElementData(vts);
        OpenDataArray(vts, "Int32", 6, "element integer tensor");
        for (int k = 0; k < vpd; k++) {
            for (int j = 0; j < vpd; j++) {
                for (int i = 0; i < vpd; i++) {
                    int e = i + j*vpd + k*vpd*vpd;
                    vts << 0*e << ' ' << 1*e << ' ' << 2*e << ' ' << 3*e << ' ' << 4*e << ' ' << 5*e << '\n';
                }
            }
        }
        CloseDataArray(vts);
        CloseElementData(vts);
    }

    // close vts file (required)
    CloseVTS(vts);
    vts.close();
    return;
}

void ppv::demo::WriteExampleVTI()
{
  std::ofstream vti("picoparaview-example.vti");
  using namespace ppv::fmt;

  int ex = 1;
  int ey = 2;
  int ez = 3;
  float hx = 1.0f;
  float hy = 2.0f;
  float hz = 3.0f;
  float ox = 0.0f;
  float oy = 0.0f;
  float oz = 0.0f;

  // open vti file (required)
  OpenVTI(vti, ex, ey, ez, hx, hy, hz, ox, oy, oz);
  OpenPiece(vti, ex, ey, ez);

  { // write node data (optional)
    OpenNodeData(vti);
    OpenDataArray(vti, "Int32", 1, "node integer scalar");
    for (int k = 0; k < ez+1; k++) {
      for (int j = 0; j < ey+1; j++) {
        for (int i = 0; i < ex+1; i++) {
          vti << i+j+k << '\n';
        }
      }
    }
    CloseDataArray(vti);
    CloseNodeData(vti);
  }

  { // write element data (optional)
    OpenElementData(vti);
    OpenDataArray(vti, "Int32", 6, "element integer tensor");
    for (int k = 0; k < ez; k++) {
      for (int j = 0; j < ey; j++) {
        for (int i = 0; i < ex; i++) {
          int e = i + j*ex + k*ex*ey;
          vti << 0*e << ' ' << 1*e << ' ' << 2*e << ' ' << 3*e << ' ' << 4*e << ' ' << 5*e << '\n';
        }
      }
    }
    CloseDataArray(vti);
    CloseElementData(vti);
  }

  // close vti file (required)
  ClosePiece(vti);
  CloseVTI(vti);
  vti.close();
  return;
}

void ppv::demo::WriteExamplePVD()
{
    using namespace ppv::fmt;
    typedef std::vector<float> Node;
    typedef std::vector<int> Element;

    const int nv = 8;
    const int ne = 12;
    const int nf = 6;
    float cubeVerts[nv][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
    int cubeEdges[ne][2] = {{0,1},{0,2},{1,3},{2,3},{0,4},{1,5},{2,6},{3,7},{4,5},{4,6},{5,7},{6,7}};
    int cubeFaces[nf][4] = {{0,2,3,1},{0,1,5,4},{1,3,7,5},{2,6,7,3},{0,4,6,2},{4,5,7,6}};

    std::ofstream pvd("picoparaview-example.pvd");
    OpenPVD( pvd );
    {
        int time = 0;
        std::vector<Node> nodes;
        nodes.reserve( nv );
        for (int n = 0; n < nv; n++) {
            nodes.push_back(Node(cubeVerts[n], cubeVerts[n]+3));
        }

        // add edges
        for (int e = 0; e < ne; e++) {
            std::vector<Element> elements;
            elements.push_back(Element(cubeEdges[e], cubeEdges[e]+2));

            std::string filename = "picoparaview-example-T0-Edge" + std::to_string(e) + ".vtu";
            std::ofstream vtu(filename);
            OpenVTU(vtu, nv, 1);
            WriteGeometryVTU(vtu, nodes, elements, LINE);
            CloseVTU(vtu);
            vtu.close();

            RecordPVDDataSet(pvd, filename, (float)(time++));
        }

        // add faces
        for (int e = 0; e < nf; e++) {
            std::vector<Element> elements;
            elements.push_back(Element(cubeFaces[e], cubeFaces[e]+4));

            std::string filename = "picoparaview-example-T0-Face" + std::to_string(e) + ".vtu";
            std::ofstream vtu(filename);
            OpenVTU(vtu, nv, 1);
            WriteGeometryVTU(vtu, nodes, elements, QUAD);
            CloseVTU(vtu);
            vtu.close();

            RecordPVDDataSet(pvd, filename, (float)(time++));
        }
    }
    ClosePVD( pvd );
}

#endif // PICOPARAVIEW_IMPLEMENTATION
