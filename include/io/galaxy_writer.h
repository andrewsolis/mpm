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

#include "Skt.hpp"
#include "bufhdr.h"

//! VTK Writer class
//! \brief VTK writer class
class GxyWriter {
  public:
      GxyWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

      void write( const std::string& data_field, unsigned step );
    
  private:
    //! Vector of nodal coordinates
    vtkSmartPointer<vtkPoints> points_;

    bool first_run = false;

};

#endif
#endif
#endif  // GXY_WRITER_H_