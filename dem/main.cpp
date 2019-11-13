#include <iostream>
#include <vector>
#include "outputWriter.h"
#include "Region.h"
#include <random>
#include "chrono"
#include "omp.h"

int main(int argc, char *argv[]) {

  double r_mean = 0.0005;
  double bin_size = r_mean*2.0*3.0;

  if(argc < 2) {
    std::cout<<"Using 3 grain d as bin size"<<std::endl;
  }
  else {
    bin_size = r_mean*2.0*std::atof(argv[1]);
    std::cout<<"Bin size of "<<std::atof(argv[1])<<" grain d"<<std::endl;
  }

  srand(100);
  std::cout << "Hello, World!" << std::endl;

  // Generate initial region
  std::vector<double> global_region_min{0.0,0.0};
  std::vector<double> global_region_max{0.40,0.24};
  int num_total_add = 64000;
  std::vector<bool> is_edge{1,1,1,1};
  Region init_region(global_region_min,global_region_max,0, is_edge,bin_size);
  init_region.generateRandomInitialPacking(r_mean,num_total_add);
  int num_total_init = init_region.m_num_grains;

  // Get number of threads for domain decomposition
  int num_threads = 1;
  #pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }
  std::cout<<"Num threads: "<<num_threads<<std::endl;
  
  // Domain decomposition
  std::vector<std::vector<double> > region_mins;
  std::vector<std::vector<double> > region_maxes;
  std::vector<Region> region_list;
  std::vector<int> domain_division;
  init_region.domainDecomposition(num_threads,region_list, domain_division);

  writeGrainsVTU("init_grains.vtu",init_region.m_q,init_region.m_v,
                   init_region.m_r,init_region.m_unique_id,init_region.m_f);

  for(int ii = 0; ii < region_list.size();ii++){
    std::cout<<"Region "<<region_list[ii].m_region_id<<" ngrains: "<<region_list[ii].m_num_grains<<std::endl;
  }

  double t = 0.0;
  double t_f = 6.0;
  double dt = 0.000005;
  int num_regions = region_list.size();
  double delta = r_mean*1.21;
  int step_num = 0;
  int freq = 200;
  //freq = static_cast<int>(1.0/dt);
  int skip = static_cast<int>(1.0/dt)/freq;

  int output_num = 0;

  // Write initial configuration
  for(int region_num = 0; region_num < num_regions; region_num++) {
    std::string prefix = "grain";
    prefix.append(std::to_string(region_num));
    prefix.append("_");
    std::string output_region_grain_str = generateOutputFileStr(prefix,freq,t_f,dt,output_num);
    writeGrainsVTU(output_region_grain_str,region_list[region_num].m_q,region_list[region_num].m_v,
                   region_list[region_num].m_r,region_list[region_num].m_unique_id,region_list[region_num].m_f);

    std::string region_prefix = "region";
    region_prefix.append(std::to_string(region_num));
    region_prefix.append(".vtu");
    writeRegionVTU(region_prefix,region_list[region_num].m_region_min,region_list[region_num].m_region_max);

  }

  output_num++;
  while(t < t_f){
    //std::cout<<"t: "<<t<<std::endl;
    step_num++;

    int right_dir = 1, left_dir = -1;
    double start_comms = omp_get_wtime();

    int num_total_across_regions = 0;
    for(int ii = 0; ii < region_list.size(); ii++) {
     region_list[ii].clearCollectionData(region_list[ii].m_surrounding_collection);
     num_total_across_regions += region_list[ii].m_num_grains;
    }

    if(num_total_across_regions != num_total_init)
      std::cout << "LOST A GRAIN: "<<num_total_across_regions<<" vs "<<num_total_init<<std::endl;
    assert(num_total_across_regions == num_total_init);

    // Left-Right communication
    for(int x_dr = -1; x_dr < 2; x_dr += 2) {

      for (int y_div = 0; y_div < domain_division[1]; y_div++) {
        for (int x_div = 0; x_div < domain_division[0]; x_div++) {
          int region_num = x_div + y_div * domain_division[0];

          // If check left boundary and is left element, or check right boundary and check right element, skip
          if( (x_div == 0 && x_dr == -1) || x_div == domain_division[0]-1 && x_dr == 1) {
            continue;
          }
          //std::cout<<"LR Region "<<region_num<<" send message to "<<region_num+x_dr<<" "<<x_dr<<std::endl;
          GrainCollection message;
          region_list[region_num].buildCommunicationMessage(message,x_dr,delta);
          region_list[region_num+x_dr].receiveCommunicationMessage(message);
        }
      }
    }

    // Up-Down communication
    for(int y_dr = -1; y_dr < 2; y_dr += 2) {

      for (int y_div = 0; y_div < domain_division[1]; y_div++) {
        for (int x_div = 0; x_div < domain_division[0]; x_div++) {
          int region_num = x_div + y_div * domain_division[0];

          // If check bottom boundary and is bottom element, or check top boundary and check top element, skip
          if( (y_div == 0 && y_dr == -1) || y_div == domain_division[1]-1 && y_dr == 1) {
            continue;
          }
          std::cout<<"UD Region "<<region_num<<" send message to "<<region_num+y_dr*domain_division[0]<<" "<<y_dr<<std::endl;
          GrainCollection message;
          region_list[region_num].buildCommunicationMessage(message,y_dr,delta);
          region_list[region_num+y_dr*domain_division[0]].receiveCommunicationMessage(message);
        }
      }
    }
    double end_comms = omp_get_wtime();

    double start_map = omp_get_wtime();
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      region_list[region].buildGrainIDMap();
      region_list[region].buildGrainCollectionIDMap(region_list[region].m_surrounding_collection);
    }
    double end_map = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].rasterizeGrainsToBins(delta);
    }
    double end_raster = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].findNeighborsFromBinsVect(delta);
    }
    double end_bins = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].buildContactList();
    }
    double end_contact = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].calculateContactForces();
    }
    double end_forces = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].applyBodyForces();
    }
    double end_body = omp_get_wtime();

    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++) {
      region_list[region].forwardEuler(dt);
    }
    double end_forward = omp_get_wtime();

    std::cout<<end_comms-start_comms<<","<<end_map-start_map<<","<<end_raster-end_map<<","<<end_bins-end_raster;
    std::cout<<","<<end_contact-end_bins<<","<<end_forces-end_contact;
    std::cout<<","<<end_body-end_forces<<","<<end_forward-end_body<<","<<end_forward-start_map<<std::endl;

    /*
    double start_sim_loop = omp_get_wtime();
    #pragma omp parallel for
    for(int region = 0; region < num_regions; region++){
      region_list[region].buildGrainIDMap();
      region_list[region].buildGrainCollectionIDMap(region_list[region].m_surrounding_collection);
      region_list[region].rasterizeGrainsToBins(delta);
      region_list[region].findNeighborsFromBinsVect(delta);
      region_list[region].buildContactList();
      region_list[region].calculateContactForces();
      region_list[region].applyBodyForces();
      region_list[region].forwardEuler(dt);
    }
    double end_sim_loop = omp_get_wtime();
    std::cout<<end_comms-start_comms<<","<<end_sim_loop-start_sim_loop<<std::endl;
    */

    bool verbose = false;

    if (verbose) {
      std::cout<<"Time: "<<t<<std::endl;
      for(int region = 0; region < num_regions; region++){
        for(int grain_num = 0; grain_num < region_list[region].m_num_grains; grain_num++){
          if(region_list[region].m_contact_list_unique[grain_num].size() > 0) {
            std::cout<<"Contact: "<<region_list[region].m_unique_id[grain_num]<<" with: ";
            for(int ii = 0; ii < region_list[region].m_contact_list[grain_num].size();ii++)
            {
              std::cout<<region_list[region].m_contact_list_unique[grain_num][ii]<<" ";
            }
            std::cout<<std::endl;

            std::cout<<"  Contact: "<<grain_num<<" with: ";
            for(int ii = 0; ii < region_list[region].m_contact_list[grain_num].size();ii++)
            {
              std::cout<<region_list[region].m_contact_list[grain_num][ii]<<" ";
            }
            std::cout<<std::endl;
          }
        }
      }
    }



    // Output to file
    if(step_num % skip == 0){
      //std::cout<<"output: "<<output_num<<" "<<t<<std::endl;
      for(int region_num = 0; region_num < num_regions; region_num++) {
        std::string prefix = "grain";
        prefix.append(std::to_string(region_num));
        prefix.append("_");
        std::string output_region_grain_str = generateOutputFileStr(prefix,freq,t_f,dt,output_num);
        writeGrainsVTU(output_region_grain_str,region_list[region_num].m_q,region_list[region_num].m_v,
                       region_list[region_num].m_r,region_list[region_num].m_unique_id,region_list[region_num].m_f);

      }
      output_num++;
    }
    t += dt;
  }


  return 0;
}