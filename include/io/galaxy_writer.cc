#include "galaxy_writer.h"

#ifdef USE_GALAXY

class GXY_Data
{
  size_t dsz = -1;
  char *data = NULL;
  float xmin = -1, xmax = 1, ymin = -1, ymax = 1, zmin = -1, zmax = 1, dmin = 0, dmax = 1;
  int step; //sender_id
  int port = -1;
  int nPts = -1;

  GXY_Data( vtkPoints* pts, const std::string& fieldname, int timestep )
  {
    bufhdr *hdr;

    vtkPointSet *ps;
    ps->SetPoints(pts);

    dsz = sizeof(int) + sizeof(bufhdr)+ nPts*4*sizeof(float);
    data = (char *)malloc(dsz);
    *(int *)data = dsz;

    hdr = (bufhdr *)(data + sizeof(int));
    float *pdst = (float *)(hdr + 1);
    float *ddst = pdst + 3*nPts;

    vtkFloatArray *pArray = vtkFloatArray::SafeDownCast(pts->GetData());
    if (pArray)
    {
      memcpy(pdst, pArray->GetVoidPointer(0), nPts * 3 * sizeof(float));
    }
    else
    {
      vtkDoubleArray *pArray = vtkDoubleArray::SafeDownCast(pts->GetData());
      if (! pArray)
      {
        std::cerr << "points have to be float or double\n";
        exit(1);
      }

      double *psrc = (double *)pArray->GetVoidPointer(0);
      for (int i = 0; i < 3*nPts; i++)
        *pdst++ = (float)*psrc++;
    }
    vtkFloatArray *fArray = vtkFloatArray::SafeDownCast(ps->GetPointData()->GetArray(fieldname.c_str()));
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
    hdr->origin        = timestep;
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
  bufhdr *hdr;

  vtkPointSet *ps;
  ps->SetPoints(points_);


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