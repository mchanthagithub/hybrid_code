#include "outputWriter.h"
//
// Created by maytee on 9/18/19.
//

void writeRegionVTU(std::string file_name, std::vector<double>& grid_min, std::vector<double>& grid_max)
{
  std::cout<<"Writing VTU"<<std::endl;
  std::ofstream fp;
  fp.open(file_name);
  //for(int ii = 0; ii < m_boundary_samples.nPoints;ii++)
  //  fp << m_boundary_samples.x.col(ii)(0)<<","<<m_boundary_samples.x.col(ii)(1)<<std::endl;
  int width = 10;
  int num_cells = 1;
  int num_pts = 4;

  fp <<"<?xml version=\"1.0\"?>\n";
  fp <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  fp <<" <UnstructuredGrid>\n";
  fp <<"  <Piece NumberOfPoints=\""<<num_pts<< "\" NumberOfCells=\""<<num_cells<<"\">\n";
  fp <<"   <CellData>\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"DEM_m_deficit\" format=\"ascii\">\n";
  int ctr = 0;

  // DEM Deficit ===================================================
  bool spaceFlag = false;
  for(int cell_num = 0; cell_num < num_cells; cell_num++)
  {
    spaceFlag = false;
    fp<<10.0<<" ";
    ctr++;
    if(ctr % width == 0)
    {
      fp<<"\n";
      spaceFlag = true;
    }
  }
  if(spaceFlag == false)
    fp<<"\n";
  fp <<"    </DataArray>\n";
  fp <<"   </CellData>\n";
  fp <<"   <Points>\n";
  fp <<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  {
    fp << grid_min[0] <<" "<< grid_min[1] <<" 0.0\n";
    fp << grid_max[0] <<" "<< grid_min[1] <<" 0.0\n";
    fp << grid_max[0] <<" "<< grid_max[1] <<" 0.0\n";
    fp << grid_min[0] <<" "<< grid_max[1] <<" 0.0\n";
  }
  fp <<"    </DataArray>\n";
  fp <<"   </Points>\n";
  fp <<"   <Cells>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  {
    fp <<0<<" "<<1<<" "<<2<<" "<<3<<std::endl;
  }
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  ctr = 0;
  for(int cell_num = 0; cell_num < num_cells; cell_num++)
  {
    spaceFlag = false;
    fp<<(cell_num+1)*4<<" ";
    ctr++;
    if(ctr % width == 0)
    {
      fp<<"\n";
      spaceFlag = true;
    }
  }
  if(spaceFlag == false)
    fp<<"\n";
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  ctr = 0;
  for(int cell_num = 0; cell_num < num_cells; cell_num++)
  {
    spaceFlag = false;
    fp<<9<<" ";
    ctr++;
    if(ctr % width == 0)
    {
      fp<<"\n";
      spaceFlag = true;
    }
  }
  if(spaceFlag == false)
    fp<<"\n";
  fp <<"    </DataArray>\n";
  fp <<"   </Cells>\n";
  fp <<"  </Piece>\n";
  fp <<" </UnstructuredGrid>\n";
  fp <<"</VTKFile>\n";

  fp.close();
}

void writeGrainsVTU(std::string file_name, std::vector<double>& q, std::vector<double>& v, std::vector<double>& r,
                    std::vector<int>& unique_id,std::vector<double>& f)
{
  //std::cout<<"Writing grains VTU"<<std::endl;
  std::ofstream fp;
  fp.open(file_name);
  //for(int ii = 0; ii < m_boundary_samples.nPoints;ii++)
  //  fp << m_boundary_samples.x.col(ii)(0)<<","<<m_boundary_samples.x.col(ii)(1)<<std::endl;

  int num_pts = q.size()/2;
  fp <<"<?xml version=\"1.0\"?>\n";
  fp <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  fp <<" <UnstructuredGrid>\n";
  fp <<"  <Piece NumberOfPoints=\""<<num_pts << "\" NumberOfCells=\"0\">\n";
  fp <<"   <Points>\n";
  fp <<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int kk = 0; kk < num_pts; kk ++)
     fp << q[kk*2]<<" "<<q[kk*2+1]<<" 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"   </Points>\n";
  fp <<"   <PointData Scalars=\"scalars\">\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"Radius\" format=\"ascii\">\n";
  for (int kk = 0; kk < num_pts; kk ++)
    fp << r[kk]<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"UniqueID\" format=\"ascii\">\n";
  for (int kk = 0; kk < num_pts; kk ++)
    fp << unique_id[kk]<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int kk = 0; kk < num_pts; kk ++)
    fp <<v[kk*2]<<" "<<v[kk*2+1]<<" 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"NetForce\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int kk = 0; kk < num_pts; kk ++)
    fp <<f[kk*2]<<" "<<f[kk*2+1]<<" 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"   </PointData>\n";

  /*
  // Prints out sigmaxx, sigmayy, sigmaxy, in that order
  fp <<"    <DataArray type=\"Float32\" Name=\"Cauchy Stress\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int p = 0; p < points.m_nPoints; p++)
    fp <<"     "<<points.m_pointStressNew[p](0,0)<<" "<<points.m_pointStressNew[p](0,1)<<" "<<points.m_pointStressNew[p](1,1)<<"\n";
  fp <<"    </DataArray>\n";

  // Prints out velocity components
  fp <<"    <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int p = 0; p < points.m_nPoints; p++)
    fp <<"     "<<points.m_vel(p*2)<<" "<<points.m_vel(p*2+1)<<" 0.0 \n";
  fp <<"    </DataArray>\n";

  // Prints out coordinates so don't need to use calculator in paraview
  fp <<"    <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int p = 0; p < points.m_nPoints; p++)
    fp <<"     "<<points.m_pos(p*2)<<" "<<points.m_pos(p*2+1)<<" 0.0 \n";
  fp <<"    </DataArray>\n";
  */
  fp <<"   <Cells>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  fp <<"    </DataArray>\n";
  fp <<"   </Cells>\n";
  fp <<"  </Piece>\n";
  fp <<" </UnstructuredGrid>\n";
  fp <<"</VTKFile>\n";

  fp.close();

}

void writeForcesVTU(std::string file_name, std::map<std::pair<int,int>,ContactGrainGrain> grain_grain_contacts,
                    std::map<std::pair<int,int>,ContactGrainWall> grain_wall_contacts)
{
  //std::cout<<"Writing grains VTU"<<std::endl;
  std::ofstream fp;
  fp.open(file_name);
  //for(int ii = 0; ii < m_boundary_samples.nPoints;ii++)
  //  fp << m_boundary_samples.x.col(ii)(0)<<","<<m_boundary_samples.x.col(ii)(1)<<std::endl;

  int num_grain_grain_contacts = grain_grain_contacts.size();
  int num_grain_wall_contacts = grain_wall_contacts.size();
  int num_total_contacts = num_grain_grain_contacts + num_grain_wall_contacts;

  // Add a dummy force to have some output (paraview doesn't like 0 elements)
  num_total_contacts++;
  fp <<"<?xml version=\"1.0\"?>\n";
  fp <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

  fp <<" <UnstructuredGrid>\n";
  fp <<"  <Piece NumberOfPoints=\""<<num_total_contacts << "\" NumberOfCells=\""<<num_total_contacts<<"\">\n";
  fp <<"   <Points>\n";
  fp <<"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(auto const& contact_grain_grain_member : grain_grain_contacts ) {
    ContactGrainGrain contact = contact_grain_grain_member.second;
    fp << contact.contact_pt[0]<<" "<<contact.contact_pt[1]<<" 0.0"<<std::endl;
  }
  for(auto const& contact_grain_wall_member : grain_wall_contacts ) {
    ContactGrainWall contact = contact_grain_wall_member.second;
    fp << contact.contact_pt[0]<<" "<<contact.contact_pt[1]<<" 0.0"<<std::endl;
  }
  fp<<" 0.0 0.0 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"   </Points>\n";
  fp <<"   <PointData Scalars=\"scalars\">\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"Overlap\" format=\"ascii\">\n";
  for(auto const& contact_grain_grain_member : grain_grain_contacts ) {
    ContactGrainGrain contact = contact_grain_grain_member.second;
    fp << contact.overlap<<std::endl;
  }
  for(auto const& contact_grain_wall_member : grain_wall_contacts ) {
    ContactGrainWall contact = contact_grain_wall_member.second;
    fp << contact.overlap<<std::endl;
  }
  fp<<" 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"Normal\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(auto const& contact_grain_grain_member : grain_grain_contacts ) {
    ContactGrainGrain contact = contact_grain_grain_member.second;
    fp << contact.n[0]<<" "<<contact.n[1]<<" 0.0"<<std::endl;
  }
  for(auto const& contact_grain_wall_member : grain_wall_contacts ) {
    ContactGrainWall contact = contact_grain_wall_member.second;
    fp << contact.n[0]<<" "<<contact.n[1]<<" 0.0"<<std::endl;
  }
  fp<<" 0.0 0.0 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Float32\" Name=\"Force\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(auto const& contact_grain_grain_member : grain_grain_contacts ) {
    ContactGrainGrain contact = contact_grain_grain_member.second;
    fp << contact.force[0]<<" "<<contact.force[1]<<" 0.0"<<std::endl;
  }
  for(auto const& contact_grain_wall_member : grain_wall_contacts ) {
    ContactGrainWall contact = contact_grain_wall_member.second;
    fp << contact.force[0]<<" "<<contact.force[1]<<" 0.0"<<std::endl;
  }
  fp<<" 0.0 0.0 0.0"<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"   </PointData>\n";
  fp <<"   <Cells>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(int ii = 0; ii < num_total_contacts; ii++)
    fp<<" "<<ii<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(int ii = 0; ii < num_total_contacts; ii++)
    fp<<" "<<ii+1<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
  for(int ii = 0; ii < num_total_contacts; ii++)
    fp<<" "<<1<<std::endl;
  fp <<"    </DataArray>\n";
  fp <<"   </Cells>\n";
  fp <<"  </Piece>\n";
  fp <<" </UnstructuredGrid>\n";
  fp <<"</VTKFile>\n";

  fp.close();

}

std::string generateOutputFileStr(std::string prefix, int freq, double t_f, double dt, int& output_num)
{
  int num_steps_total = static_cast<int>(t_f / dt);
  int num_digits = std::to_string(num_steps_total).length();
  std::string output_str = prefix;
  int num_digits_output_num = std::to_string(output_num).length();
  for(int ii = num_digits - num_digits_output_num; ii > 0; ii--){
    output_str.append("0");
  }
  output_str.append(std::to_string(output_num));
  output_str.append(".vtu");
  return output_str;
}
