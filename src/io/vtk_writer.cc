#include "vtk_writer.h"

#ifdef USE_VTK

//! VTK Writer class Constructor with coordniates
//! \param[in] coordinate Point coordinates
//! \param[in] node_pairs Node ID pairs to form elements
VtkWriter::VtkWriter(
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

//! Write coordinates
//! \param[in] filename Output file to write geometry
void VtkWriter::write_geometry(const std::string& filename) {

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points_);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

  writer->SetDataModeToBinary();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif

  writer->Write();
}

//! Write Tensor data
void VtkWriter::write_tensor_point_data(
    const std::string& filename,
    const std::vector<Eigen::Matrix<double, 6, 1>>& data,
    const std::string& data_field) {

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points_);

  // Create an array to hold distance information
  auto tensordata = vtkSmartPointer<vtkDoubleArray>::New();
  tensordata->SetNumberOfComponents(9);
  tensordata->SetNumberOfTuples(pdata->GetNumberOfPoints());
  tensordata->SetName(data_field.c_str());

  // Evaluate the signed distance function at all of the grid points
  for (vtkIdType id = 0; id < pdata->GetNumberOfPoints(); ++id) {
    const double* vdata = data.at(id).data();
    tensordata->InsertTuple9(id, vdata[0], vdata[3], vdata[5], vdata[3],
                             vdata[1], vdata[4], vdata[5], vdata[4], vdata[2]);
  }

  // Add the SignedDistances to the grid
  pdata->GetPointData()->SetTensors(tensordata);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

  writer->SetDataModeToBinary();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif

  writer->SetCompressor(vtkZLibDataCompressor::New());
  writer->Write();
}

//! Write vector data
void VtkWriter::write_vector_point_data(
    const std::string& filename, const std::vector<Eigen::Vector3d>& data,
    const std::string& data_field) {

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points_);

  // Create an array to hold distance information
  auto vectordata = vtkSmartPointer<vtkDoubleArray>::New();
  vectordata->SetNumberOfComponents(3);
  vectordata->SetNumberOfTuples(pdata->GetNumberOfPoints());
  vectordata->SetName(data_field.c_str());

  // Evaluate the signed distance function at all of the grid points
  for (vtkIdType id = 0; id < pdata->GetNumberOfPoints(); ++id) {
    const double* vdata = data.at(id).data();
    vectordata->SetTuple(id, vdata);
  }

  // Add the SignedDistances to the grid
  pdata->GetPointData()->SetVectors(vectordata);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

  writer->SetDataModeToBinary();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif

  writer->SetCompressor(vtkZLibDataCompressor::New());
  writer->Write();
}

//! Write scalar data
void VtkWriter::write_scalar_point_data(const std::string& filename,
                                        const std::vector<double>& data,
                                        const std::string& data_field) {

  // Create an array to hold distance information
  auto scalardata = vtkSmartPointer<vtkDoubleArray>::New();
  scalardata->SetNumberOfComponents(1);
  scalardata->SetName(data_field.c_str());

  // Add value to scalarfield
  for (auto value : data) scalardata->InsertNextValue(value);

  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();
  // Add the points to the dataset
  pdata->SetPoints(points_);
  // Add the scalar data to the points
  pdata->GetPointData()->SetScalars(scalardata);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

  writer->SetDataModeToBinary();

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif

  writer->SetCompressor(vtkZLibDataCompressor::New());
  writer->Write();
}

//! VTK Write mesh
void VtkWriter::write_mesh(
    const std::string& filename,
    const std::vector<Eigen::Vector3d>& coordinates,
    const std::vector<std::array<mpm::Index, 2>>& node_pairs) {
  // Assing points
  auto points = vtkSmartPointer<vtkPoints>::New();
  unsigned id = 0;
  for (const auto& coordinate : coordinates) {
    const double* point = coordinate.data();
    points->InsertPoint(id, point);
    ++id;
  }
  // Assign elements
  // Create a cell array to store the lines
  auto lines = vtkSmartPointer<vtkCellArray>::New();

  for (const auto& node_pair : node_pairs) {
    // Create a line between the two points
    auto line = vtkSmartPointer<vtkLine>::New();
    // the SetId(A,B) call is the following
    line->GetPointIds()->SetId(0, node_pair[0]);
    // A = the id of the point relative to the line - this can only be 0 or 1.
    // B = the index into the vtkPoints object of the point that you would
    // like to set the Ath point to.
    line->GetPointIds()->SetId(1, node_pair[1]);

    lines->InsertNextCell(line);
  }
  // Create a polydata to store everything in it
  auto pdata = vtkSmartPointer<vtkPolyData>::New();

  // Add the points to the dataset
  pdata->SetPoints(points);

  // Add the lines to the dataset
  pdata->SetLines(lines);

  // Write file
  auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  writer->SetFileName(filename.c_str());

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(pdata);
#else
  writer->SetInputData(pdata);
#endif
  writer->Write();
}

//! Write Parallel VTK file
void VtkWriter::write_parallel_vtk(const std::string& filename,
                                   const std::string& attribute, int mpi_size,
                                   unsigned step, unsigned max_steps,
                                   unsigned ncomponents) {

  // If the number of components is 1, set as scalar or vector / tensor
  std::string data_type;
  switch (ncomponents) {
    case 1:
      data_type = "Scalars";
      break;
    case 9:
      data_type = "Tensors";
      break;
    default:
      data_type = "Vectors";
      break;
  }

  std::string ppolydata =
      "<?xml version=\"1.0\"?>\n<VTKFile type=\"PPolyData\" version=\"0.1\" "
      "byte_order=\"LittleEndian\" "
      "compressor=\"vtkZLibDataCompressor\">\n<PPolyData "
      "GhostLevel=\"0\">\n\t<PPointData " +
      data_type + "=\"" + attribute +
      "\">\n\t\t<PDataArray "
      "type=\"Float64\" Name=\"" +
      attribute + "\"";
  if (ncomponents != 1)
    ppolydata += " NumberOfComponents=\"" + std::to_string(ncomponents) + "\"";
  ppolydata +=
      "/>\n\t</"
      "PPointData>\n\n\t<PPoints>\n\t\t<PDataArray "
      "type=\"Float32\" Name=\"Points\" "
      "NumberOfComponents=\"3\"/>\n\t</PPoints>\n";

  for (unsigned i = 0; i < mpi_size; ++i) {
    std::stringstream file_name;
    file_name.str(std::string());
    file_name << attribute;
    const std::string rank_size =
        "-" + std::to_string(i) + "_" + std::to_string(mpi_size) + "-";
    file_name << rank_size;

    file_name.fill('0');
    int digits = log10(max_steps) + 1;
    file_name.width(digits);
    file_name << step;
    file_name << ".vtp";
    ppolydata += "\n\t<Piece Source=\"" + file_name.str() + "\"/>";
  }
  ppolydata += "\n</PPolyData>\n\n</VTKFile>";

  // Write parallel VTK file
  std::ofstream pvtk;
  pvtk.open(filename);
  pvtk << ppolydata;
  pvtk.close();
}

void VtkWriter::setup_galaxy(int mpi_size, int mpi_rank)
{

  // master_socket = (mpi_rank == 0) ? new gxy::ClientSkt("localhost", 1900) : NULL;

  // int status = 0;
  // int sz;
  // char *hoststring = NULL;

  // if ( master_socket->Connect() )
  // {
  //   status = 1;

  //   if (status && !master_socket->Send(std::string("box")))
  //     status = 0;

  //   float box[] = {gxmin, gymin, gzmin, gxmax, gymax, gzmax};
  //   if (status && !master_socket->Send((void *)box, 6*sizeof(float)))
  //     status = 0;

  //   std::stringstream ss;
  //   ss << "nsenders " << n_senders;
  //   if (status && !master_socket->Send(ss.str()))
  //     status = 0;

  //   if (status && !master_socket->Send(std::string("sendhosts")))
  //     status = 0;
  
  //   if (status && ((hoststring = master_socket->Receive()) == NULL))
  //     status = 0;

  //   if (status && !master_socket->Send(std::string("go")))
  //     status = 0;

  //   master_socket->Disconnect();

  // }
}

void VtkWriter::create_data(const std::string& data_field, unsigned step)
{
  sender_id = step;
  bufhdr *hdr;

  vtkPointSet *ps;
  ps->SetPoints(points_);

  vtkPoints *vpts = ps->GetPoints();
  nPts = ps->GetNumberOfPoints();

  dsz = sizeof(int) + sizeof(bufhdr)+ nPts*4*sizeof(float);
  data = (char *)malloc(dsz);
  *(int *)data = dsz;

  hdr = (bufhdr *)(data + sizeof(int));
  float *pdst = (float *)(hdr + 1);
  float *ddst = pdst + 3*nPts;

  vtkFloatArray *pArray = vtkFloatArray::SafeDownCast(vpts->GetData());
  if (pArray)
    memcpy(pdst, pArray->GetVoidPointer(0), nPts * 3 * sizeof(float));
  else
  {
    vtkDoubleArray *pArray = vtkDoubleArray::SafeDownCast(vpts->GetData());
    if (! pArray)
    {
      std::cerr << "points have to be float or double\n";
      exit(1);
    }

    double *psrc = (double *)pArray->GetVoidPointer(0);
    for (int i = 0; i < 3*nPts; i++)
      *pdst++ = (float)*psrc++;
  }

  vtkFloatArray *fArray = vtkFloatArray::SafeDownCast(ps->GetPointData()->GetArray(data_field.c_str()));
  if (fArray)
  {
    float *fsrc = (float *)fArray->GetVoidPointer(0);

    if (fArray->GetNumberOfComponents() == 1)
    {
      memcpy(ddst, fArray->GetVoidPointer(0), nPts * sizeof(float));
    }
    else if (fArray->GetNumberOfComponents() == 3)
    {
      for (int i = 0; i < nPts; i++)
      {
        double a = fsrc[0]*fsrc[0] + fsrc[1]*fsrc[1] + fsrc[2]*fsrc[2];
        *ddst++ = sqrt(a);
        fsrc += 3;
      }
    }
    else
    {
      std::cerr << "data have to be scalar or 3-vector\n";
      exit(1);
    }
  }
  else
  {
    vtkDoubleArray *dArray = vtkDoubleArray::SafeDownCast(ps->GetPointData()->GetArray("displacements"));
    if (! dArray)
    {
      std::cerr << "data have to be float or double\n";
      exit(1);
    }

    double *dsrc = (double *)dArray->GetVoidPointer(0);

    if (dArray->GetNumberOfComponents() == 1)
    {
      for (int i = 0; i < nPts; i++)
        *ddst++ = *dsrc++;
    }
    else if (dArray->GetNumberOfComponents() == 3)
    {
      for (int i = 0; i < nPts; i++)
      {
        double a = dsrc[0]*dsrc[0] + dsrc[1]*dsrc[1] + dsrc[2]*dsrc[2];
        *ddst++ = sqrt(a);
        dsrc += 3;
      }
    }
    else
    {
      std::cerr << "data have to be scalar or 3-vector\n";
      exit(1);
    }
  }

  hdr->type          = bufhdr::Particles; 
  hdr->origin        = sender_id;
  hdr->has_data      = true;
  hdr->npoints       = nPts;
  hdr->nconnectivity = 0;

  float *pptr = (float *)(data + sizeof(int) + sizeof(bufhdr));
  float *dptr = pptr + 3*nPts;

  xmin = xmax = pptr[0];
  ymin = ymax = pptr[1];
  zmin = zmax = pptr[2];
  dmin = dmax = *dptr;

  for (int i = 0; i < nPts; i++)
  {
    if (xmin > pptr[0]) xmin = pptr[0];
    if (xmax < pptr[0]) xmax = pptr[0];
    if (ymin > pptr[1]) ymin = pptr[1];
    if (ymax < pptr[1]) ymax = pptr[1];
    if (zmin > pptr[2]) zmin = pptr[2];
    if (zmax < pptr[2]) zmax = pptr[2];
    if (dmin > *dptr) dmin = *dptr;
    if (dmax < *dptr) dmax = *dptr;
    pptr += 3;
    dptr ++;
  }
}

#endif
