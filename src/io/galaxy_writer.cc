#include "galaxy_writer.h"

#ifdef USE_GALAXY


GXY_Data::GXY_Data( vtkPoints* pts, const std::string& fieldname, unsigned timestep )
{
  bufhdr *hdr;

  step = timestep;

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

void GxyWriter::set_destination( string hst, int prt ) 
{
  host = hst;
  port = prt;
}


void GxyWriter::setup(const GXY_Data* data)
{
  master_socket = new gxy::ClientSkt( master_host, base_port );

  int status = 0;
  int sz;
  char *hoststring = NULL;

  if( master_socket && master_socket->Connect() )
  {
    if (status && !master_socket->Send(std::string("box")))
      status = 0;

    float box[] = {data->xmin, data->ymin, data->zmin, data->xmax, data->ymax, data->zmax};
    if (status && !master_socket->Send((void *)box, 6*sizeof(float)))
      status = 0;

    std::stringstream ss;
    ss << "nsenders " << n_senders;
    if (status && !master_socket->Send(ss.str()))
      status = 0;

    if (status && !master_socket->Send(std::string("sendhosts")))
      status = 0;

    if (status && ((hoststring = master_socket->Receive()) == NULL))
      status = 0;

    if (status && !master_socket->Send(std::string("go")))
      status = 0;

    master_socket->Disconnect();
  
  }

  sz = status ? strlen(hoststring) : -1;
  
  
  if(status)
  {

    std::vector<char *> hosts;
    hosts.push_back(hoststring);

    for (char *c = hoststring; *c; c++)
    {
      if (*c == ':') *c = ' ';
      else if (*c == ';')
      {
        *c = '\0';
        if (*(c+1) != '\0')
          hosts.push_back(c+1);
      }
    }

    char host[256]; int port;

    stringstream ss(hosts[0]);

    ss >> host >> port;

    set_destination( host, port );
  }

}

void GxyWriter::write( const std::string& data_field, unsigned step )
{
    GXY_Data *gxy_data = new GXY_Data::GXY_Data(points_, data_field, step );

    if( !first_run )
    {
      setup( gxy_data ); 
    }

    if( port == -1 )
    {
      std::cerr << "Can't send data without knowing where!\n";
    }
    else
    {
      gxy::ClientSkt c(host, port);
      c.Connect() && c.Send(gxy_data->data, gxy_data->dsz);
    }

}


#endif