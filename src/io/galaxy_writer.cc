#include "galaxy_writer.h"

#ifdef USE_GALAXY


GXY_Data::GXY_Data( vtkSmartPointer<vtkPoints> pts, const std::string& fieldname, unsigned timestep )
{
  std::cout << "inside gxy_data constructor" << std::endl;
  bufhdr *hdr;

  step = timestep;
  std::cout << "set step" << std::endl;

  std::cout << "classname : " << pts->GetClassName() << std::endl;
  std::cout << "# points : " << pts->GetNumberOfPoints() << std::endl;


 vtkPointSet *ps;
  ps->SetPoints( (vtkPoints *) pts);

  std::cout << "set ps" << std::endl;

  dsz = sizeof(int) + sizeof(bufhdr)+ nPts*4*sizeof(float);
  data = (char *)malloc(dsz);
  *(int *)data = dsz;

  hdr = (bufhdr *)(data + sizeof(int));
  float *pdst = (float *)(hdr + 1);
  float *ddst = pdst + 3*nPts;

  std::cout << "creating parray  from pts" << std::endl;

  vtkFloatArray *pArray = vtkFloatArray::SafeDownCast(pts->GetData());
  if (pArray)
  {
    std::cout << "Inside pArray if" << std::endl;
    memcpy(pdst, pArray->GetVoidPointer(0), nPts * 3 * sizeof(float));
  }
  else
  {
    std::cout << "inside pArray else" << std::endl;
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
  std::cout << "creating fArray" << std::endl;

  vtkFloatArray *fArray = vtkFloatArray::SafeDownCast(ps->GetPointData()->GetArray(fieldname.c_str()));
  if (fArray)
  {
    std::cout << "inide fArray if" << std::endl;
    float *fsrc = (float *)fArray->GetVoidPointer(0);

    if (fArray->GetNumberOfComponents() == 1)
    {
      std::cout << "inside fArray if if" << std::endl;

      memcpy(ddst, fArray->GetVoidPointer(0), nPts * sizeof(float));
    }
    else if (fArray->GetNumberOfComponents() == 3)
    {
      std::cout << "inside else if" << std::endl;

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
    std::cout << "inside else for fArray" << std::endl;

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

  std::cout << "creating hdr" << std::endl;

  hdr->type          = bufhdr::Particles;
  hdr->origin        = timestep;
  hdr->has_data      = true;
  hdr->npoints       = nPts;
  hdr->nconnectivity = 0;

  std::cout << "creating pointers" << std::endl;

  float *pptr = (float *)(data + sizeof(int) + sizeof(bufhdr));
  float *dptr = pptr + 3*nPts;

  std::cout << "creating minmax" << std::endl;

  xmin = xmax = pptr[0];
  ymin = ymax = pptr[1];
  zmin = zmax = pptr[2];
  dmin = dmax = *dptr;


  std::cout << "min max for loop" << std::endl;
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
  std::cout << "assinging host and port" << std::endl;
  host = hst;
  port = prt;
}


void GxyWriter::setup(const GXY_Data* data)
{
  std::cout << "creating new master socket" << std::endl;

  master_socket = new gxy::ClientSkt( master_host, base_port );

  int status = 0;
  int sz;
  char *hoststring = NULL;

  std::cout << "attempting to connect..." << std::endl;

  if( master_socket && master_socket->Connect() )
  {
    std::cout << "inside if block of setup. sending string" << std::endl;
    if (status && !master_socket->Send(std::string("box")))
      status = 0;

    std::cout << "sending min/max data" << std::endl;

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

    std::cout << "made it to end of if block. disconnect...." << std::endl;

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
    
    std::cout << "calling set destination for host/port" << std::endl;
    set_destination( host, port );
  }

}

void GxyWriter::write( const std::string& data_field, unsigned step )
{
    std::cout << "creating new galaxy data object" << std::endl;

    GXY_Data *gxy_data = new GXY_Data::GXY_Data(points_, data_field, step );

    if( !first_run )
    {
      std::cout << "calling setup to setup galaxy server" << std::endl;
      setup( gxy_data ); 
    }

    if( port == -1 )
    {
      std::cerr << "Can't send data without knowing where!\n";
    }
    else
    {
      std::cout << "creating temp socket to send data" << std::endl;
      gxy::ClientSkt c(host, port);
      std::cout << "attempting to send data to galaxy" << std::endl;
      c.Connect() && c.Send(gxy_data->data, gxy_data->dsz);
    }

}


#endif
