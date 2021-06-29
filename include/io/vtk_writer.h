#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#ifdef USE_VTK
#include <fstream>

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

#include <fstream>
#include <string>
#include <vector>

#include "data_types.h"
#include "Skt.hpp"
#include "bufhdr.h"

//! VTK Writer class
//! \brief VTK writer class
class VtkWriter {
 public:
  // Constructor with coordinates
  VtkWriter(const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates);

  //! Write coordinates
  void write_geometry(const std::string& filename);

  //! \brief Write scalar data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Scalar field data
  //! \param[in] data_field Field name ("pdstrain")
  void write_scalar_point_data(const std::string& filename,
                               const std::vector<double>& data,
                               const std::string& data_field);
  //! Write vector data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Vector data
  //! \param[in] data_field Field name ("Displacement", "Forces")
  void write_vector_point_data(const std::string& filename,
                               const std::vector<Eigen::Vector3d>& data,
                               const std::string& data_fields);

  //! Write tensor data
  //! \param[in] filename Output file to write geometry
  //! \param[in] data Tensor data
  //! \param[in] data_field Field name ("Displacement", "Forces")
  void write_tensor_point_data(
      const std::string& filename,
      const std::vector<Eigen::Matrix<double, 6, 1>>& data,
      const std::string& data_fields);

  //! Write mesh
  //! \param[in] filename Mesh VTP file
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] node_pairs Node pairs ids
  void write_mesh(const std::string& filename,
                  const std::vector<Eigen::Vector3d>& coordinates,
                  const std::vector<std::array<mpm::Index, 2>>& node_pairs);

  //! Write Parallel VTK file
  //! \param[in] filename Mesh PVTP file name
  //! \param[in] attribute VTK data attribute to be written
  //! \param[in] mpi_size Number of MPI tasks
  //! \param[in] step Current time step
  //! \param[in] max_steps Maximum number of steps in the simulation
  //! \param[in] ncomponents Number of components to write
  void write_parallel_vtk(const std::string& filename,
                          const std::string& attribute, int mpi_size,
                          unsigned step, unsigned max_steps,
                          unsigned ncomponents = 3);

int setup_galaxy(int mpi_rank = 0, int mpi_size = 1);

void create_data(const std::string& data_field, unsigned step);

bool send_data();


 private:
  //! Vector of nodal coordinates
  vtkSmartPointer<vtkPoints> points_;

  //! galaxy variables
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
#endif  // VTK_WRITER_H_
