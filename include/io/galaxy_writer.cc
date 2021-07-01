#include "galaxy_writer.h"

#ifdef USE_GALAXY

class gxy_data
{
  size_t dsz = -1;
  char *data = NULL;
  float xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1, dmin = 0, dmax = 1;
  int step; //sender_id
  int port = -1;
  int nPts = -1;
};

//! Gxy Writer class Constructor with coordniates
//! \param[in] coordinate Point coordinates
//! \param[in] node_pairs Node ID pairs to form elements
GxyWriter::GxyWriter(
    const std::vector<Eigen::Matrix<double, 3, 1>>& coordinates) {
  // Assign points
  points_ = vtkSmartPointer<vtkPoints>::New();
  unsigned long long id = 0;
  for (const auto& coordinate : coordinates) {
    const double* point = coordinate.data();
    points_->InsertPoint(id, point);
    ++id;
  }
}

void GxyWriter::create_data( const std::string& data_field )
{
  
}

void GxyWriter::write( const std::string& data_field )
{
    // create data to send

    if( !first_run )
    {
        // setup master_socket
    }
}


#endif