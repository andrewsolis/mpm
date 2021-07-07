#ifndef GXY_WRITER_H_
#define GXY_WRITER_H_

#ifdef USE_VTK
#ifdef USE_GALAXY

#include <Eigen/Dense>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTensorGlyph.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkFloatArray.h>

#include "Skt.h"
#include "bufhdr.h"

class GXY_Data
{
  public:
    size_t dsz = -1;
    char *data = NULL;
    float xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1, dmin = 0, dmax = 1;
    int step; //sender_id
    int nPts = -1;

    GXY_Data( vtkPoints* pts, const std::string& fieldname, unsigned timestep );
};

//! Galaxy Writer class
//! \brief Galaxy writer class
class GxyWriter {
  public:
    GxyWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

    void setup( const GXY_Data *data );

    void write( const std::string& data_field, unsigned step );

    void set_destination(string hst, int prt );
    
    //! Vector of nodal coordinates
    vtkSmartPointer<vtkPoints> points_;

    bool first_run = false;

  private:
    gxy::ClientSkt *master_socket;

    string master_host = "localhost";
    int    base_port = 1900;
    
    string host = "";
    int    port = -1;
    
    int    n_senders = 1;

};

#endif
#endif
#endif  // GXY_WRITER_H_