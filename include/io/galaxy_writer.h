#ifndef GXY_WRITER_H_
#define GXY_WRITER_H_

#ifdef USE_VTK

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

#include "Skt.hpp"

//! VTK Writer class
//! \brief VTK writer class
class GxyWriter {
  public:
      GxyWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

      void create_data( const std::string& data_field, unsigned step );
    
  private:
    //! Vector of nodal coordinates
    vtkSmartPointer<vtkPoints> points_;

      int sender_id;
      size_t dsz = -1;
      char *data = NULL;
      gxy::ClientSkt *master_socket;
      float xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1, dmin = 0, dmax = 1;
      int n_senders = 1; 
      int nPts = -1;

      std::string dest_host;
      int dest_port;
};

#endif
#endif  // GXY_WRITER_H_