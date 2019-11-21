//
// Created by maytee on 10/23/19.
//

#include <set>
#include <cfloat>
#include <algorithm>
#include <iomanip>
#include "Region.h"

int Region::m_total_num_grains = 0;
int Region::m_id_tracker = 0;
bool verbose = false;

GrainCollection::GrainCollection() {
  num_grains_in_collection = 0;
}

// Should only be used for input setting purposes as it does not update unique id
void GrainCollection::addGrainToCollection(double x_in, double y_in, double theta_in, double r_in, double vx_in,
                                           double vy_in, double omega_in)
{
  q.push_back(x_in); q.push_back(y_in);
  theta.push_back(theta_in);
  r.push_back(r_in);
  v.push_back(vx_in); v.push_back(vy_in);
  omega.push_back(omega_in);
  f.push_back(0.0); f.push_back(0.0);
  std::vector<int> dummy;
  neighbor_list.push_back(dummy);
  contact_list.push_back(dummy);
  contact_list_unique.push_back(dummy);
  unique_id.push_back(0); // WATCH FOR THI
  num_grains_in_collection++;
}

Region::Region(InputSettings& settings) {
  // Generate regions
  std::vector<Region> region_list;
  int num_create_regions = settings.list_of_region_bounds.size();
  std::vector<double> global_region_min = {DBL_MAX, DBL_MAX};
  std::vector<double> global_region_max = {DBL_MIN, DBL_MIN};

  // Create all sub-regions
  for(int ii = 0; ii < num_create_regions; ii++) {
    Region new_region({settings.list_of_region_bounds[ii][0],settings.list_of_region_bounds[ii][1]},
                      {settings.list_of_region_bounds[ii][2],settings.list_of_region_bounds[ii][3]},ii,
                      {0,0,0,0},settings.bin_size);
    double r_mean = settings.list_radii_dist_params[ii][0];
    double r_min = r_mean*(1.0 - settings.list_radii_dist_params[ii][1]);
    double r_max = r_mean*(1.0 + settings.list_radii_dist_params[ii][2]);
    int num_add = static_cast<int>(settings.list_radii_dist_params[ii][3]);
    new_region.generateRandomInitialPacking(r_mean,r_min,r_max,num_add);
    region_list.push_back(new_region);

    if(settings.list_of_region_bounds[ii][0] < global_region_min[0])
      global_region_min[0] = settings.list_of_region_bounds[ii][0];
    if(settings.list_of_region_bounds[ii][1] < global_region_min[1])
      global_region_min[1] = settings.list_of_region_bounds[ii][1];

    if(settings.list_of_region_bounds[ii][2] > global_region_max[0])
      global_region_max[0] = settings.list_of_region_bounds[ii][2];
    if(settings.list_of_region_bounds[ii][3] > global_region_max[1])
      global_region_max[1] = settings.list_of_region_bounds[ii][3];
  }

  // Initialize this region
  m_region_min = global_region_min;
  m_region_max = global_region_max;
  m_num_grains = 0;
  m_bin_size = settings.bin_size;
  m_num_bins_x = static_cast<int>((m_region_max[0] - m_region_min[0])/m_bin_size);
  m_num_bins_y = static_cast<int>((m_region_max[1] - m_region_min[1])/m_bin_size);
  m_bin_list.resize(m_num_bins_x*m_num_bins_y);
  m_region_id = 0;
  m_is_edge = {true,true,true,true};


  for(const Region & reg : region_list) {
    for(int grain_num = 0; grain_num < reg.m_num_grains; grain_num++) {
      double x = reg.m_q[grain_num*2]; double y = reg.m_q[grain_num*2+1];
      double theta = reg.m_theta[grain_num]; double r = reg.m_r[grain_num];
      double vx = reg.m_v[grain_num*2]; double vy = reg.m_v[grain_num*2+1];
      double omega = reg.m_omega[grain_num];
      addGrainToRegion(x,y,theta,r,vx,vy,omega);
    }
  }

  for(const std::vector<double> & grain : settings.list_of_individual_grains) {
    double x = grain[0]; double y = grain[1];
    double theta = grain[2]; double r = grain[3];
    double vx = grain[4]; double vy = grain[5];
    double omega = grain[6];
    addGrainToRegion(x,y,theta,r,vx,vy,omega);
  }

}

Region::Region(std::vector<double> in_min, std::vector<double> in_max, int id, std::vector<bool> is_edge, double bin_size)
{
  m_region_min = in_min;
  m_region_max = in_max;
  m_num_grains = 0;
  double r_mean = 0.0005;
  m_bin_size = bin_size;
  m_num_bins_x = static_cast<int>((m_region_max[0] - m_region_min[0])/m_bin_size);
  m_num_bins_y = static_cast<int>((m_region_max[1] - m_region_min[1])/m_bin_size);
  m_bin_list.resize(m_num_bins_x*m_num_bins_y);
  m_region_id = id;
  m_is_edge = is_edge;
  std::cout<<"binsx: "<<m_num_bins_x<<" binsy: "<<m_num_bins_y<<std::endl;

  std::cout<<"Region: "<<m_region_id<<" edge: "<<m_is_edge[0]<<" "<<m_is_edge[1]<<" "<<m_is_edge[2]<<" "<<m_is_edge[3]<<std::endl;
}

Region::Region()
{
  m_region_min.push_back(0.0); m_region_min.push_back(0.0);
  m_region_max.push_back(2.0); m_region_max.push_back(2.0);
  m_num_grains = 0;
  double r_mean = 0.0005;
  m_bin_size = r_mean*2.0*3.0;
  m_num_bins_x = static_cast<int>((m_region_max[0] - m_region_min[0])/m_bin_size);
  m_num_bins_y = static_cast<int>((m_region_max[1] - m_region_min[1])/m_bin_size);
  m_bin_list.resize(m_num_bins_x*m_num_bins_y);
  m_region_id = 0;
}

void Region::generateRandomInitialPacking(double r_mean, double r_min, double r_max, int num_total_add)
{
  std::cout<<"r range: "<<r_min<<" "<<r_mean<<" "<<r_max<<std::endl;
  double x_in = r_max + m_region_min[0];
  double y_in = r_max + m_region_min[1];
  double v_max = r_min*50.0, v_min = r_min*(-50.0);
  for(int grain = 0; grain < num_total_add; grain++) {
    double f = (double)rand() / RAND_MAX;
    double r_in = r_min + f * (r_max - r_min);
    double f_vx = (double)rand() / RAND_MAX;
    double f_vy = (double)rand() / RAND_MAX;
    double v_x = v_min + f_vx * (v_max - v_min);
    double v_y = v_min + f_vy * (v_max - v_min);

    if(x_in + r_max > m_region_max[0]) {
      x_in = r_max + m_region_min[0];
      y_in += r_max*2.0;
    }
    addGrainToRegion(x_in,y_in,0.0,r_in,v_x,v_y,0.0);
    x_in += r_max*2.0;
  }
}

void Region::addGrainToRegion(double x_in, double y_in, double theta_in, double r_in, double vx_in, double vy_in, double omega_in)
{
  m_q.push_back(x_in); m_q.push_back(y_in);
  m_theta.push_back(theta_in);
  m_r.push_back(r_in);
  m_v.push_back(vx_in); m_v.push_back(vy_in);
  m_omega.push_back(omega_in);
  m_f.push_back(0.0); m_f.push_back(0.0);
  std::vector<int> dummy;
  m_neighbor_list.push_back(dummy);
  m_contact_list.push_back(dummy);
  m_contact_list_unique.push_back(dummy);
  m_unique_id.push_back(m_id_tracker);
  m_total_num_grains++;
  m_id_tracker++;
  m_num_grains++;
}

void Region::addGrainToRegion(double x_in, double y_in, double theta_in, double r_in, double vx_in, double vy_in, double omega_in,
        double fx_in, double fy_in, std::vector<int> neighbor_list_in, std::vector<int> contact_list_in, int unique_id)
{
  m_q.push_back(x_in); m_q.push_back(y_in);
  m_theta.push_back(theta_in);
  m_r.push_back(r_in);
  m_v.push_back(vx_in); m_v.push_back(vy_in);
  m_omega.push_back(omega_in);
  m_f.push_back(fx_in); m_f.push_back(fy_in);
  m_neighbor_list.push_back(neighbor_list_in);
  m_contact_list.push_back(contact_list_in);
  m_contact_list_unique.push_back(contact_list_in);
  m_unique_id.push_back(unique_id);
  //m_total_num_grains++;
  // ID tracker doesn't change as it is keeping track across all Regions
  //m_id_tracker++;
  m_num_grains++;
}

void Region::removeGrain(int idx)
{
  //std::cout<<"Remove "<<idx<<" "<<m_unique_id[idx]<<" from region "<<m_region_id<<std::endl;
  m_q.erase(m_q.begin()+idx*2,m_q.begin()+idx*2+1+1);
  m_theta.erase(m_theta.begin()+idx);
  m_r.erase(m_r.begin()+idx);
  m_v.erase(m_v.begin()+idx*2,m_v.begin()+idx*2+1+1);
  m_omega.erase(m_omega.begin()+idx);
  m_f.erase(m_f.begin()+idx*2,m_f.begin()+idx*2+1+1);
  m_neighbor_list.erase(m_neighbor_list.begin()+idx);
  m_contact_list.erase(m_contact_list.begin()+idx);
  m_contact_list_unique.erase(m_contact_list_unique.begin()+idx);
  m_unique_id.erase(m_unique_id.begin()+idx);
  // WARNING, BECAUSE M_TOTAL_NUM_GRAINS IS STATIC ACROSS ALL MEMBERS THIS MAY NOT BE THREAD SAFE
  m_total_num_grains--;
  m_num_grains--;
}

void Region::addGrainToCollection(GrainCollection& collection, int idx)
{
  collection.q.push_back(m_q[idx*2]); collection.q.push_back(m_q[idx*2+1]);
  collection.theta.push_back(m_theta[idx]);
  collection.r.push_back(m_r[idx]);
  collection.v.push_back(m_v[idx*2]); collection.v.push_back(m_v[idx*2+1]);
  collection.omega.push_back(m_omega[idx]);
  collection.f.push_back(m_f[idx*2]); collection.f.push_back(m_f[idx*2+1]);
  collection.neighbor_list.push_back(m_neighbor_list[idx]);
  collection.contact_list.push_back(m_contact_list[idx]);
  collection.contact_list_unique.push_back(m_contact_list_unique[idx]);
  collection.unique_id.push_back(m_unique_id[idx]);
  collection.region_id.push_back(m_region_id);
  collection.num_grains_in_collection++;
}

void Region::addGrainFromCollectionToCollection(GrainCollection &collection_in, GrainCollection &collection_receive, int grain_idx)
{
  collection_receive.q.push_back(collection_in.q[grain_idx*2]); collection_receive.q.push_back(collection_in.q[grain_idx*2+1]);
  collection_receive.theta.push_back(collection_in.theta[grain_idx]);
  collection_receive.r.push_back(collection_in.r[grain_idx]);
  collection_receive.v.push_back(collection_in.v[grain_idx*2]); collection_receive.v.push_back(collection_in.v[grain_idx*2+1]);
  collection_receive.omega.push_back(collection_in.omega[grain_idx]);
  collection_receive.f.push_back(collection_in.f[grain_idx*2]); collection_receive.f.push_back(collection_in.f[grain_idx*2+1]);
  collection_receive.neighbor_list.push_back(collection_in.neighbor_list[grain_idx]);
  collection_receive.contact_list.push_back(collection_in.contact_list[grain_idx]);
  collection_receive.contact_list_unique.push_back(collection_in.contact_list_unique[grain_idx]);
  collection_receive.unique_id.push_back(collection_in.unique_id[grain_idx]);
  collection_receive.region_id.push_back(collection_in.region_id[grain_idx]);
  collection_receive.num_grains_in_collection++;
}

void Region::buildGrainIDMap()
{
  m_unique_to_local_map.clear();
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    m_unique_to_local_map[m_unique_id[grain_num]] = grain_num;
  }
}

void Region::buildGrainCollectionIDMap(GrainCollection& collection)
{
  collection.unique_to_local_map.clear();
  for(int grain_num = 0; grain_num < collection.num_grains_in_collection; grain_num++) {
    collection.unique_to_local_map[collection.unique_id[grain_num]] = grain_num;
  }
}

void Region::domainDecomposition(int num_regions, std::vector<Region>& region_list, std::vector<int>& domain_division)
{
  // For now break up domain into rectangles with only single face-face connections, so must have
  // even number of regions
  assert(num_regions == 1 || num_regions%2 == 0);

  // Find avg x coord, avg y coord, and max and min x and y coords to determine how to split up domain
  double avg_x = 0.0, avg_y = 0.0;
  double max_r = DBL_MIN;
  double min_x = DBL_MAX, min_y = DBL_MAX, max_x = DBL_MIN, max_y = DBL_MIN;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    double x_coord = m_q[grain_num*2];
    double y_coord = m_q[grain_num*2+1];
    avg_x += x_coord;
    avg_y += y_coord;

    if(x_coord < min_x)
      min_x = x_coord;
    if(x_coord > max_x)
      max_x = x_coord;
    if(y_coord < min_y)
      min_y = y_coord;
    if(y_coord > max_y)
      max_y = y_coord;

    if(m_r[grain_num] > max_r)
      max_r = m_r[grain_num];
  }
  avg_x /= m_num_grains;
  avg_y /= m_num_grains;

  double margin = 3.0;
  min_x -= max_r*margin;
  min_y -= max_r*margin;
  max_x += max_r*margin;
  max_y += max_r*margin;
  //min_x = m_region_min[0];
  //min_y = m_region_min[1];
  //max_x = m_region_max[0];
  max_y = m_region_max[1];
  double x_length = max_x - min_x;
  double y_length = max_y - min_y;

  std::vector<std::pair<int,int>> divisors;
  for(int ii = 1; ii <= num_regions; ii++) {
    if(num_regions%ii == 0) {
      std::cout<<"ii: "<<ii<<","<<num_regions/ii<<std::endl;
      std::pair<int,int> in_pair = std::make_pair(ii,num_regions/ii);
      divisors.push_back(in_pair);
    }
  }

  std::vector<int> count_in_region(num_regions,0);
  double min_ratio = DBL_MAX;
  int min_idx = divisors.size();
  for(int ii = 0; ii < divisors.size(); ii++) {
    int x_div = divisors[ii].first;
    int y_div = divisors[ii].second;
    double dx = x_length/x_div;
    double dy = y_length/y_div;
    for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
      for(int y_idx = 0; y_idx < y_div; y_idx++) {
        for(int x_idx = 0; x_idx < x_div; x_idx++) {
          if(m_q[grain_num*2] >= min_x + x_idx*dx && m_q[grain_num*2] < min_x + (x_idx+1)*dx &&
              m_q[grain_num*2+1] >= min_y + y_idx*dy && m_q[grain_num*2+1] < min_y + (y_idx+1)*dy) {
            count_in_region[x_idx+y_idx*x_div] += 1;
          }
        }
      }
    }
    int max_num = *std::max_element(count_in_region.begin(),count_in_region.end());
    int min_num = *std::min_element(count_in_region.begin(),count_in_region.end());
    if(min_num == 0)
      continue;
    double ratio = static_cast<double>(max_num)/min_num;
    std::cout<<"ratio: "<<ratio<<std::endl;
    if(ratio < min_ratio) {
      min_ratio = ratio;
      min_idx = ii;
    }
  }

  std::cout<<"minidx: "<<min_idx<<" Divisions: "<<divisors[min_idx].first<<","<<divisors[min_idx].second<<std::endl;

  int x_div_final = divisors[min_idx].first;
  int y_div_final = divisors[min_idx].second;
  double dx_final = x_length/x_div_final;
  double dy_final = y_length/y_div_final;
  std::vector<std::vector<int> >grain_indexes_in_regions;
  grain_indexes_in_regions.resize(num_regions);
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    for(int y_idx = 0; y_idx < y_div_final; y_idx++) {
      for(int x_idx = 0; x_idx < x_div_final; x_idx++) {
        if(m_q[grain_num*2] >= min_x + x_idx*dx_final && m_q[grain_num*2] < min_x + (x_idx+1)*dx_final &&
            m_q[grain_num*2+1] >= min_y + y_idx*dy_final && m_q[grain_num*2+1] < min_y + (y_idx+1)*dy_final) {
          grain_indexes_in_regions[x_idx+y_idx*x_div_final].push_back(grain_num);
        }
      }
    }
  }
  domain_division.push_back(divisors[min_idx].first);
  domain_division.push_back(divisors[min_idx].second);

  for(int y_regions = 0 ; y_regions < divisors[min_idx].second; y_regions++) {
    for(int x_regions = 0 ; x_regions < divisors[min_idx].first; x_regions++) {
      std::vector<double> region_min = {min_x + x_regions*dx_final,min_y+y_regions*dy_final};
      std::vector<double> region_max = {min_x + (x_regions+1)*dx_final,min_y+(y_regions+1)*dy_final};
      int region_idx = x_regions + y_regions*divisors[min_idx].first;
      std::cout<<"region: "<<region_idx<<" has "<<grain_indexes_in_regions[region_idx].size()<<" grains "<<std::endl;
      std::cout<<"      minx: "<<region_min[0]<<" miny: "<<region_min[1]<<std::endl;
      std::cout<<"      maxx: "<<region_max[0]<<" maxy: "<<region_max[1]<<std::endl;
      std::vector<bool> is_edge{false,false,false,false};
      if(x_regions == 0)
        is_edge[0] = true;
      if(x_regions == divisors[min_idx].first-1)
        is_edge[1] = true;
      if(y_regions == 0)
        is_edge[2] = true;
      if(y_regions == divisors[min_idx].second-1)
        is_edge[3] = true;

      Region region(region_min,region_max,region_idx,is_edge,m_bin_size);

      for(int ii = 0; ii < grain_indexes_in_regions[region_idx].size(); ii++) {
        int idx = grain_indexes_in_regions[region_idx][ii];
        region.addGrainToRegion(m_q[idx*2],m_q[idx*2+1],m_theta[idx], m_r[idx],m_v[idx*2],m_v[idx*2+1],m_omega[idx], m_f[idx*2],m_f[idx*2+1],
            m_neighbor_list[idx],m_contact_list[idx],m_unique_id[idx]);
      }
      region_list.push_back(region);
    }
  }
}

void Region::findNeighborsBruteForce(double delta)
{
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    m_neighbor_list[grain_num].clear();
    std::vector<int> list;
    for(int check_num = grain_num+1; check_num < m_num_grains; check_num++){
      double dist_squared = (m_q[grain_num*2]-m_q[check_num*2])*(m_q[grain_num*2]-m_q[check_num*2])
                          + (m_q[grain_num*2+1]-m_q[check_num*2+1])*(m_q[grain_num*2+1]-m_q[check_num*2+1]);
      double r_squared = (m_r[grain_num] + m_r[check_num]+delta) * (m_r[grain_num] + m_r[check_num]+delta);
      if(dist_squared <= r_squared){
        //std::cout<<"grain: "<<grain_num<<" check: "<<check_num<<std::endl;
        list.push_back(check_num);
      }
    }
    m_neighbor_list[grain_num] = list;
  }
}

void Region::rasterizeGrainsToBins(double delta)
{
  for(int ii = 0; ii < m_bin_list.size();ii++)
    m_bin_list[ii].clear();
  //m_bin_list.resize(m_num_bins_x*m_num_bins_y);

  // By looping over grains and then inserting them into bins, this ensures that
  // the smallest grain ID will be first in the bin. Note that this can change if
  // using the unique ID and not just the grain number
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){

    // Cache quantities
    double q_x = m_q[grain_num*2];
    double q_y = m_q[grain_num*2+1];
    double r = m_r[grain_num];
    // Inflate r by a little bit to avoid floating point issues
    r *= 1.001;

    // Generate AABB for this grain
    std::vector<double> aabb_min{q_x-r,q_y-r};
    std::vector<double> aabb_max{q_x+r,q_y+r};

    // Find IDs of bins that the AABB intersects
    std::vector<int> bin_id_min{std::max( std::min(static_cast<int>((aabb_min[0]-m_region_min[0])/m_bin_size),m_num_bins_x-1),0),
                                std::max( std::min(static_cast<int>((aabb_min[1]-m_region_min[1])/m_bin_size),m_num_bins_y-1),0)};
    std::vector<int> bin_id_max{std::min( std::max(static_cast<int>((aabb_max[0]-m_region_min[0])/m_bin_size),0),m_num_bins_x-1),
                                std::min( std::max(static_cast<int>((aabb_max[1]-m_region_min[1])/m_bin_size),0),m_num_bins_y-1)};

    /*
    std::vector<int> bin_id_min{std::max(static_cast<int>((aabb_min[0]-m_region_min[0])/m_bin_size),0),
                                std::max(static_cast<int>((aabb_min[1]-m_region_min[1])/m_bin_size),0)};
    std::vector<int> bin_id_max{std::min(static_cast<int>((aabb_max[0]-m_region_min[0])/m_bin_size),m_num_bins_x-1),
                                std::min(static_cast<int>((aabb_max[1]-m_region_min[1])/m_bin_size),m_num_bins_y-1)};
    */

    assert(bin_id_max[0] >= bin_id_min[0]);
    assert(bin_id_max[1] >= bin_id_min[1]);

    if(verbose) {
      std::cout<<"grain: "<<grain_num<<std::endl;
      std::cout<<"  minx "<<bin_id_min[0]<<" miny "<<bin_id_min[1]<<std::endl;
      std::cout<<"  maxx "<<bin_id_max[0]<<" maxy "<<bin_id_max[1]<<std::endl;
    }


    for(int kk = bin_id_min[1]; kk <= bin_id_max[1]; kk++){
      for(int ii = bin_id_min[0]; ii <= bin_id_max[0]; ii++) {
        int bin_idx = ii + kk*m_num_bins_x;
        m_bin_list[bin_idx].push_back(grain_num);
      }
    }
  }

  // Loop over surrounding grains
  for(int surrounding_grain_num = 0; surrounding_grain_num < m_surrounding_collection.num_grains_in_collection; surrounding_grain_num++){

    // Cache quantities
    double q_x = m_surrounding_collection.q[surrounding_grain_num*2];
    double q_y = m_surrounding_collection.q[surrounding_grain_num*2+1];
    double r = m_surrounding_collection.r[surrounding_grain_num];
    // Inflate r by a little bit to avoid floating point issues
    r *= 1.001;
    // Generate AABB for this grain
    std::vector<double> aabb_min{q_x-r,q_y-r};
    std::vector<double> aabb_max{q_x+r,q_y+r};

    // Find IDs of bins that the AABB intersects
    // Need additional min and max checks because neighboring grains will be above and below the region
    std::vector<int> bin_id_min{std::max( std::min(static_cast<int>((aabb_min[0]-m_region_min[0])/m_bin_size),m_num_bins_x-1),0),
                                std::max( std::min(static_cast<int>((aabb_min[1]-m_region_min[1])/m_bin_size),m_num_bins_y-1),0)};
    std::vector<int> bin_id_max{std::min( std::max(static_cast<int>((aabb_max[0]-m_region_min[0])/m_bin_size),0),m_num_bins_x-1),
                                std::min( std::max(static_cast<int>((aabb_max[1]-m_region_min[1])/m_bin_size),0),m_num_bins_y-1)};

    assert(bin_id_max[0] >= bin_id_min[0]);
    assert(bin_id_max[1] >= bin_id_min[1]);

    if(verbose) {
      std::cout<<"neighbor grain: "<<surrounding_grain_num<<std::endl;
      std::cout<<" calcminx "<<static_cast<int>((aabb_min[0]-m_region_min[0])/m_bin_size);
      std::cout<<" calcminy "<<static_cast<int>((aabb_min[1]-m_region_min[1])/m_bin_size)<<std::endl;
      std::cout<<" calcmaxx "<<static_cast<int>((aabb_max[0]-m_region_min[0])/m_bin_size);
      std::cout<<" calcmaxy "<<static_cast<int>((aabb_max[1]-m_region_min[1])/m_bin_size)<<std::endl;
      std::cout<<"  minx "<<bin_id_min[0]<<" miny "<<bin_id_min[1]<<std::endl;
      std::cout<<"  maxx "<<bin_id_max[0]<<" maxy "<<bin_id_max[1]<<std::endl;
    }


    for(int kk = bin_id_min[1]; kk <= bin_id_max[1]; kk++){
      for(int ii = bin_id_min[0]; ii <= bin_id_max[0]; ii++) {
        int bin_idx = ii + kk*m_num_bins_x;
        m_bin_list[bin_idx].push_back(surrounding_grain_num+m_num_grains);
      }
    }
  }


  if(verbose) {
    for(int bin_idx = 0; bin_idx < m_num_bins_x*m_num_bins_y; bin_idx++){
      if(m_bin_list[bin_idx].size() == 0)
        continue;
      std::cout<<"bin: "<<bin_idx<<std::endl;
      for(int ii = 0; ii < m_bin_list[bin_idx].size();ii++) {
        std::cout<<m_bin_list[bin_idx][ii]<<" ";
      }
      std::cout<<std::endl;

    }
  }

}

void Region::findNeighborsFromBins(double delta)
{
  std::vector<std::set<int>> check_set;
  //check_set.resize(m_num_grains+m_surrounding_collection.num_grains_in_collection);
  check_set.resize(m_num_grains);
  // Loop through all bins to see what grains intersect them
  for(int bin_idx = 0; bin_idx < m_num_bins_x*m_num_bins_y; bin_idx++){
    // Take the first grain in the bin as the one that collects the contacts
    // This assumes the grains are in ascending ID order in the bins
    for(int grain_num_1 = 0; grain_num_1 < m_bin_list[bin_idx].size(); grain_num_1++) {
      if(m_bin_list[bin_idx][grain_num_1] >= m_num_grains)
        continue;
      for(int grain_num_2 = grain_num_1+1; grain_num_2 < m_bin_list[bin_idx].size(); grain_num_2++){
        check_set[m_bin_list[bin_idx][grain_num_1]].insert(m_bin_list[bin_idx][grain_num_2]);
      //m_neighbor_list[m_bin_list[bin_idx][0]].push_back(m_bin_list[bin_idx][grain_num]);
      }
    }
  }
  for(int ii = 0; ii < m_num_grains; ii++) {
    m_neighbor_list[ii].clear();
    if(check_set[ii].size() > 0) {
      //std::vector<int> list(check_set[ii].begin(),check_set[ii].end());
      //m_neighbor_list[ii] = list;
      m_neighbor_list[ii].assign(check_set[ii].begin(),check_set[ii].end());
    }
  }
}

void Region::findNeighborsFromBinsVect(double delta)
{
  std::vector<std::vector<int>> check_set;
  //check_set.resize(m_num_grains+m_surrounding_collection.num_grains_in_collection);
  check_set.resize(m_num_grains);
  // Loop through all bins to see what grains intersect them
  for(int bin_idx = 0; bin_idx < m_num_bins_x*m_num_bins_y; bin_idx++){
    // Take the first grain in the bin as the one that collects the contacts
    // This assumes the grains are in ascending ID order in the bins
    for(int grain_num_1 = 0; grain_num_1 < m_bin_list[bin_idx].size(); grain_num_1++) {
      if(m_bin_list[bin_idx][grain_num_1] >= m_num_grains)
        continue;
      for(int grain_num_2 = grain_num_1+1; grain_num_2 < m_bin_list[bin_idx].size(); grain_num_2++){
        check_set[m_bin_list[bin_idx][grain_num_1]].push_back(m_bin_list[bin_idx][grain_num_2]);
        //m_neighbor_list[m_bin_list[bin_idx][0]].push_back(m_bin_list[bin_idx][grain_num]);
      }
    }
  }
  for(int ii = 0; ii < m_num_grains; ii++) {
    m_neighbor_list[ii].clear();
    if(check_set[ii].size() > 0) {
      std::sort(check_set[ii].begin(),check_set[ii].end());
      check_set[ii].erase( std::unique( check_set[ii].begin(),check_set[ii].end() ), check_set[ii].end() );
      m_neighbor_list[ii] = check_set[ii];
    }
  }
}

// Because surrounding grains will always have larger indexes (not unique IDs) than the grains in the region,
// they will always show up in the neighbor lists of the grains in the region, and thus do not need to loop over
// surrounding grains
void Region::buildContactList()
{
  int count = 0;
  int count_extra = 0;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    count += m_neighbor_list[grain_num].size();
    //std::cout<<"grain: "<<grain_num<<" pos: ";
    //for(int check = 0; check < m_neighbor_list[grain_num].size(); check++)
    //  std::cout<<m_neighbor_list[grain_num][check]<<" ";
    //std::cout<<std::endl;
  }
  //std::cout<<"# possible contacts: "<<count<<std::endl;

  m_contacts.clear();
  int actual_contacts = 0;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    std::vector<int> contact_list;
    std::vector<int> contact_list_unique;
    m_contact_list[grain_num].clear();
    m_contact_list_unique[grain_num].clear();
    for(int check_num = 0; check_num < m_neighbor_list[grain_num].size(); check_num++){
      int check_id = m_neighbor_list[grain_num][check_num];
      int surround_id = check_id - m_num_grains;
      int surround_unique_id = -1;

      double dist_squared;
      double r_squared;
      if(check_id < m_num_grains) {
        dist_squared = (m_q[grain_num*2]-m_q[check_id*2])*(m_q[grain_num*2]-m_q[check_id*2])
                          + (m_q[grain_num*2+1]-m_q[check_id*2+1])*(m_q[grain_num*2+1]-m_q[check_id*2+1]);
        r_squared = (m_r[grain_num] + m_r[check_id]) * (m_r[grain_num] + m_r[check_id]);
      } else {
        dist_squared = (m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])*(m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])
                          + (m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1])*(m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1]);
        r_squared = (m_r[grain_num] + m_surrounding_collection.r[surround_id]) * (m_r[grain_num] + m_surrounding_collection.r[surround_id]);
        surround_unique_id = m_surrounding_collection.unique_id[surround_id];
      }

      if(dist_squared <= r_squared){
        contact_list.push_back(check_id);
        ContactGrainGrain newContact;
        newContact.grain_idx_1 = m_unique_id[grain_num];
        newContact.overlap = sqrt(r_squared)-sqrt(dist_squared);

        if(check_id >= m_num_grains) {
          contact_list_unique.push_back(surround_unique_id);
          newContact.grain_idx_2 = surround_unique_id;
        }
        else {
          contact_list_unique.push_back(m_unique_id[check_id]);
          newContact.grain_idx_2 = m_unique_id[check_id];
        }
        actual_contacts++;



        if(verbose) {
          if(check_id < m_num_grains) {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" x1 "<<m_q[grain_num*2]<<" x2 "<<m_q[check_id*2];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_r[check_id]<<std::endl;
          }
          else {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" y1 "<<m_q[grain_num*2+1]<<" y2 "<<m_surrounding_collection.q[surround_id*2+1];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_surrounding_collection.r[surround_id]<<std::endl;
          }
        }
      }
    }
    m_contact_list[grain_num] = contact_list;
    m_contact_list_unique[grain_num] = contact_list_unique;
  }
  if(verbose) {
    for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
      if(m_contact_list[grain_num].size() > 0) {
        std::cout<<"Contact: "<<grain_num<<" with: ";
        for(int ii = 0; ii < m_contact_list[grain_num].size();ii++)
        {
          std::cout<<m_contact_list[grain_num][ii]<<" ";
        }
        std::cout<<std::endl;
      }
    }
  }

  //std::cout<<"#actual contacts: "<<actual_contacts<<std::endl;
}

void Region::buildContactListWithContactObjects()
{
  int count = 0;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    count += m_neighbor_list[grain_num].size();
    //std::cout<<"grain: "<<grain_num<<" pos: ";
    //for(int check = 0; check < m_neighbor_list[grain_num].size(); check++)
    //  std::cout<<m_neighbor_list[grain_num][check]<<" ";
    //std::cout<<std::endl;
  }
  //std::cout<<"# possible contacts: "<<count<<std::endl;

  m_contacts.clear();
  int actual_contacts = 0;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    std::vector<int> contact_list;
    std::vector<int> contact_list_unique;
    m_contact_list[grain_num].clear();
    m_contact_list_unique[grain_num].clear();
    int unique_id_1 = m_unique_id[grain_num];
    double x_coord_1 = m_q[grain_num*2];
    double y_coord_1 = m_q[grain_num*2+1];
    double r_1 = m_r[grain_num];
    for(int check_num = 0; check_num < m_neighbor_list[grain_num].size(); check_num++){
      int check_id = m_neighbor_list[grain_num][check_num];
      int surround_id = check_id - m_num_grains;
      int surround_unique_id = -1;

      double x_coord_2, y_coord_2, r_2;
      int unique_id_2;
      double dist_squared;
      double r_squared;
      if(check_id < m_num_grains) {
        unique_id_2 = m_unique_id[check_id];
        x_coord_2 = m_q[check_id*2];
        y_coord_2 = m_q[check_id*2+1];
        r_2 = m_r[check_id];
      } else {
        unique_id_2 = m_surrounding_collection.unique_id[surround_id];
        x_coord_2 = m_surrounding_collection.q[surround_id*2];
        y_coord_2 = m_surrounding_collection.q[surround_id*2+1];
        r_2 = m_surrounding_collection.r[surround_id];
        surround_unique_id = m_surrounding_collection.unique_id[surround_id];
      }

      dist_squared = (x_coord_1 - x_coord_2)*(x_coord_1 - x_coord_2) + (y_coord_1 - y_coord_2)*(y_coord_1 - y_coord_2);
      r_squared = (r_1 + r_2)*(r_1 + r_2);

      double dist = sqrt(dist_squared);

      // Contact found
      if(dist_squared <= r_squared){

        // Always make sure smallest ID goes first as convention, and to stay consistent
        if(unique_id_1 >= unique_id_2) {
          std::swap(x_coord_1,x_coord_2);
          std::swap(y_coord_1,y_coord_2);
          std::swap(r_1,r_2);
          std::swap(unique_id_1,unique_id_2);
        }

        ContactGrainGrain newContact;
        newContact.grain_idx_1 = unique_id_1;
        newContact.grain_idx_2 = unique_id_2;
        newContact.overlap = sqrt(r_squared)-dist;
        contact_list_unique.push_back(surround_unique_id);

        std::vector<double> n{0.0,0.0};
        std::vector<double> contact_pt{0.0,0.0};
        n[0] = (x_coord_2 - x_coord_1)/dist;
        n[1] = (y_coord_2 - y_coord_1)/dist;
        contact_pt[0] = x_coord_1 + n[0]*dist*0.5;
        contact_pt[1] = y_coord_1 + n[1]*dist*0.5;
        newContact.n = n;
        newContact.contact_pt = contact_pt;

        std::pair<int,int> id_pair = std::make_pair(unique_id_1,unique_id_2);
        // Check to see if this contact existed at the last timestep
        if(m_contacts_cache.count(id_pair) == 1) {
          // Found an old contact, take its integrated history
          newContact.s = m_contacts_cache.at(id_pair).s;
        } else {
          newContact.s = 0.0;
        }
        m_contacts.insert(std::make_pair(id_pair,newContact));

        contact_list.push_back(check_id);
        actual_contacts++;

        if(verbose) {
          if(check_id < m_num_grains) {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" x1 "<<m_q[grain_num*2]<<" x2 "<<m_q[check_id*2];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_r[check_id]<<std::endl;
          }
          else {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" y1 "<<m_q[grain_num*2+1]<<" y2 "<<m_surrounding_collection.q[surround_id*2+1];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_surrounding_collection.r[surround_id]<<std::endl;
          }
        }
      }
    }
    m_contact_list[grain_num] = contact_list;
    m_contact_list_unique[grain_num] = contact_list_unique;
  }
  m_contacts_cache = m_contacts;

  if(verbose) {
    for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
      if(m_contact_list[grain_num].size() > 0) {
        std::cout<<"Contact: "<<grain_num<<" with: ";
        for(int ii = 0; ii < m_contact_list[grain_num].size();ii++)
        {
          std::cout<<m_contact_list[grain_num][ii]<<" ";
        }
        std::cout<<std::endl;
      }
    }
  }

  //std::cout<<"#actual contacts: "<<actual_contacts<<std::endl;
}

double dotProduct(std::vector<double> a, std::vector<double> b) {
  double result = 0.0;
  for(int ii = 0; ii < a.size(); ii++) {
    result += a[ii]*b[ii];
  }
  return result;
}

void Region::calculateContactForces()
{
  // Zero out force
  std::fill(m_f.begin(), m_f.end(), 0.0);

  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    for(int contact_num = 0; contact_num < m_contact_list[grain_num].size(); contact_num++){
      int check_id = m_contact_list[grain_num][contact_num];

      if(check_id < m_num_grains) {
        double dist_squared = (m_q[grain_num*2]-m_q[check_id*2])*(m_q[grain_num*2]-m_q[check_id*2])
                          + (m_q[grain_num*2+1]-m_q[check_id*2+1])*(m_q[grain_num*2+1]-m_q[check_id*2+1]);
        double r_squared = (m_r[grain_num] + m_r[check_id]) * (m_r[grain_num] + m_r[check_id]);

        double dist = sqrt(dist_squared);
        double r_tot = m_r[grain_num]+m_r[check_id];

        double overlap = r_tot-dist;
        std::vector<double> n{0.0,0.0};
        n[0] = (m_q[check_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_q[check_id*2+1] - m_q[grain_num*2+1])/dist;

        std::vector<double> v_rel{0.0,0.0};
        v_rel[0] = m_v[check_id*2] - m_v[grain_num*2];
        v_rel[1] = m_v[check_id*2+1] - m_v[grain_num*2+1];
        std::vector<double> v_n{0.0,0.0};
        std::vector<double> v_t{0.0,0.0};
        v_t[0] = v_rel[0] - dotProduct(n,v_rel)*n[0];
        v_t[1] = v_rel[1] - dotProduct(n,v_rel)*n[1];
        v_n[0] = v_rel[0] - v_t[0];
        v_n[1] = v_rel[1] - v_t[1];

        m_f[grain_num*2] -= n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[grain_num*2+1] -= n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        m_f[check_id*2] += n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[check_id*2+1] += n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        double fx = n[0]*overlap*k_n;
        double fy = n[1]*overlap*k_n;
        if(verbose) {
          std::cout<<std::fixed;
          std::cout<<std::setprecision(10);
          std::cout<<"F btwn "<<m_unique_id[grain_num]<<","<<m_unique_id[check_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" fx "<<fx<<" fy "<<fy<<std::endl;

          //std::cout<<"F btwn "<<m_unique_id[grain_num]<<","<<m_unique_id[check_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" f1x "<<m_f[grain_num*2]<<" f1y "<<m_f[grain_num*2+1]<<std::endl;
          //std::cout<<" f2x "<<m_f[check_id*2]<<" f2y "<<m_f[check_id*2+1]<<std::endl;
        }


      } else {
        int surround_id = check_id - m_num_grains;
        double dist_squared = (m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])*(m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])
                          + (m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1])*(m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1]);
        double r_squared = (m_r[grain_num] + m_surrounding_collection.r[surround_id]) * (m_r[grain_num] + m_surrounding_collection.r[surround_id]);

        double r_tot = m_r[grain_num]+m_surrounding_collection.r[surround_id];
        double dist = sqrt(dist_squared);
        double overlap = r_tot-dist;

        //double overlap = sqrt(r_squared - dist_squared);
        std::vector<double> n{0.0,0.0};
        n[0] = (m_surrounding_collection.q[surround_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_surrounding_collection.q[surround_id*2+1] - m_q[grain_num*2+1])/dist;

        std::vector<double> v_rel{0.0,0.0};
        v_rel[0] = m_surrounding_collection.v[surround_id*2] - m_v[grain_num*2];
        v_rel[1] = m_surrounding_collection.v[surround_id*2+1] - m_v[grain_num*2+1];
        std::vector<double> v_n{0.0,0.0};
        std::vector<double> v_t{0.0,0.0};
        v_t[0] = v_rel[0] - dotProduct(n,v_rel)*n[0];
        v_t[1] = v_rel[1] - dotProduct(n,v_rel)*n[1];
        v_n[0] = v_rel[0] - v_t[0];
        v_n[1] = v_rel[1] - v_t[1];

        // Only apply forces on grain in the region
        m_f[grain_num*2] -= n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[grain_num*2+1] -= n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        double fx = n[0]*overlap*k_n;
        double fy = n[1]*overlap*k_n;

        if(verbose) {
          std::cout<<std::fixed;
          std::cout<<std::setprecision(10);
          std::cout<<"F surround btwn "<<m_unique_id[grain_num]<<","<<m_surrounding_collection.unique_id[surround_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" fx "<<fx<<" fy "<<fy<<std::endl;

          //std::cout<<"F surround btwn "<<m_unique_id[grain_num]<<","<<m_surrounding_collection.unique_id[surround_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<std::endl;
          //std::cout<<" fx "<<m_f[grain_num*2]<<" fy "<<m_f[grain_num*2+1]<<std::endl;
        }
      }

    }
    // Check against floor
    if(m_q[grain_num*2+1]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2+1] += (m_r[grain_num]-m_q[grain_num*2+1])*k_n;
    }

    // Check against ceiling
    if(m_q[grain_num*2+1]+m_r[grain_num] >= 0.24)
    {
      m_f[grain_num*2+1] += (0.24-m_q[grain_num*2+1]-m_r[grain_num])*k_n;
    }

    // Check against left side
    if(m_q[grain_num*2]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2] += (m_r[grain_num]-m_q[grain_num*2])*k_n;
    }
    // Check against right side
    if(m_q[grain_num*2]+m_r[grain_num] >= 0.40)
    {
      m_f[grain_num*2] += (0.40-m_q[grain_num*2]-m_r[grain_num])*k_n;
    }
  }
}

void Region::calculateContactForcesWithContactObjects()
{
  // Zero out force
  std::fill(m_f.begin(), m_f.end(), 0.0);

  
  
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    for(int contact_num = 0; contact_num < m_contact_list[grain_num].size(); contact_num++){
      int check_id = m_contact_list[grain_num][contact_num];

      if(check_id < m_num_grains) {
        double dist_squared = (m_q[grain_num*2]-m_q[check_id*2])*(m_q[grain_num*2]-m_q[check_id*2])
                              + (m_q[grain_num*2+1]-m_q[check_id*2+1])*(m_q[grain_num*2+1]-m_q[check_id*2+1]);
        double r_squared = (m_r[grain_num] + m_r[check_id]) * (m_r[grain_num] + m_r[check_id]);

        double dist = sqrt(dist_squared);
        double r_tot = m_r[grain_num]+m_r[check_id];

        double overlap = r_tot-dist;
        std::vector<double> n{0.0,0.0};
        n[0] = (m_q[check_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_q[check_id*2+1] - m_q[grain_num*2+1])/dist;

        std::vector<double> v_rel{0.0,0.0};
        v_rel[0] = m_v[check_id*2] - m_v[grain_num*2];
        v_rel[1] = m_v[check_id*2+1] - m_v[grain_num*2+1];
        std::vector<double> v_n{0.0,0.0};
        std::vector<double> v_t{0.0,0.0};
        v_t[0] = v_rel[0] - dotProduct(n,v_rel)*n[0];
        v_t[1] = v_rel[1] - dotProduct(n,v_rel)*n[1];
        v_n[0] = v_rel[0] - v_t[0];
        v_n[1] = v_rel[1] - v_t[1];

        m_f[grain_num*2] -= n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[grain_num*2+1] -= n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        m_f[check_id*2] += n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[check_id*2+1] += n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        double fx = n[0]*overlap*k_n;
        double fy = n[1]*overlap*k_n;
        if(verbose) {
          std::cout<<std::fixed;
          std::cout<<std::setprecision(10);
          std::cout<<"F btwn "<<m_unique_id[grain_num]<<","<<m_unique_id[check_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" fx "<<fx<<" fy "<<fy<<std::endl;

          //std::cout<<"F btwn "<<m_unique_id[grain_num]<<","<<m_unique_id[check_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" f1x "<<m_f[grain_num*2]<<" f1y "<<m_f[grain_num*2+1]<<std::endl;
          //std::cout<<" f2x "<<m_f[check_id*2]<<" f2y "<<m_f[check_id*2+1]<<std::endl;
        }


      } else {
        int surround_id = check_id - m_num_grains;
        double dist_squared = (m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])*(m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])
                              + (m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1])*(m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1]);
        double r_squared = (m_r[grain_num] + m_surrounding_collection.r[surround_id]) * (m_r[grain_num] + m_surrounding_collection.r[surround_id]);

        double r_tot = m_r[grain_num]+m_surrounding_collection.r[surround_id];
        double dist = sqrt(dist_squared);
        double overlap = r_tot-dist;

        //double overlap = sqrt(r_squared - dist_squared);
        std::vector<double> n{0.0,0.0};
        n[0] = (m_surrounding_collection.q[surround_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_surrounding_collection.q[surround_id*2+1] - m_q[grain_num*2+1])/dist;

        std::vector<double> v_rel{0.0,0.0};
        v_rel[0] = m_surrounding_collection.v[surround_id*2] - m_v[grain_num*2];
        v_rel[1] = m_surrounding_collection.v[surround_id*2+1] - m_v[grain_num*2+1];
        std::vector<double> v_n{0.0,0.0};
        std::vector<double> v_t{0.0,0.0};
        v_t[0] = v_rel[0] - dotProduct(n,v_rel)*n[0];
        v_t[1] = v_rel[1] - dotProduct(n,v_rel)*n[1];
        v_n[0] = v_rel[0] - v_t[0];
        v_n[1] = v_rel[1] - v_t[1];

        // Only apply forces on grain in the region
        m_f[grain_num*2] -= n[0]*overlap*k_n - 0.5*eta_n*v_n[0];
        m_f[grain_num*2+1] -= n[1]*overlap*k_n - 0.5*eta_n*v_n[1];

        double fx = n[0]*overlap*k_n;
        double fy = n[1]*overlap*k_n;

        if(verbose) {
          std::cout<<std::fixed;
          std::cout<<std::setprecision(10);
          std::cout<<"F surround btwn "<<m_unique_id[grain_num]<<","<<m_surrounding_collection.unique_id[surround_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" fx "<<fx<<" fy "<<fy<<std::endl;

          //std::cout<<"F surround btwn "<<m_unique_id[grain_num]<<","<<m_surrounding_collection.unique_id[surround_id]<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<std::endl;
          //std::cout<<" fx "<<m_f[grain_num*2]<<" fy "<<m_f[grain_num*2+1]<<std::endl;
        }
      }

    }
    // Check against floor
    if(m_q[grain_num*2+1]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2+1] += (m_r[grain_num]-m_q[grain_num*2+1])*k_n;
    }

    // Check against ceiling
    if(m_q[grain_num*2+1]+m_r[grain_num] >= 0.1)
    {
      m_f[grain_num*2+1] += (0.1-m_q[grain_num*2+1]-m_r[grain_num])*k_n;
    }

    // Check against left side
    if(m_q[grain_num*2]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2] += (m_r[grain_num]-m_q[grain_num*2])*k_n;
    }
    // Check against right side
    if(m_q[grain_num*2]+m_r[grain_num] >= 0.05)
    {
      m_f[grain_num*2] += (0.05-m_q[grain_num*2]-m_r[grain_num])*k_n;
    }

  }

}

void Region::applyBodyForces()
{
  double g = -9.81;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    m_f[grain_num*2+1] += g*M_PI*m_r[grain_num]*m_r[grain_num]*m_rho;
  }
}

void Region::forwardEuler(double delta_t)
{
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++) {
    double mass = M_PI*m_r[grain_num]*m_r[grain_num]*m_rho;
    m_v[grain_num*2] += delta_t*m_f[grain_num*2]/mass;
    m_v[grain_num*2+1] += delta_t*m_f[grain_num*2+1]/mass;
    m_q[grain_num*2] += delta_t*m_v[grain_num*2];
    m_q[grain_num*2+1] += delta_t*m_v[grain_num*2+1];
    //m_f[grain_num*2] = 0.0;
    //m_f[grain_num*2+1] = 0.0;
  }
}

void Region::clearCollectionData(GrainCollection& collection)
{
  collection.q.clear();
  collection.r.clear();
  collection.v.clear();
  collection.f.clear();
  collection.unique_id.clear();
  collection.neighbor_list.clear();
  collection.contact_list.clear();
  collection.contact_list_unique.clear();
  collection.region_id.clear();
  collection.num_grains_in_collection = 0;
}

// dir = 1 -> x, 2 -> y, 3 -> z. negative value: check west/south boundary, positive check east/north boundary
void Region::buildCommunicationMessage(GrainCollection& message, int& dir, double& cutoff)
{
  // Clear old messages; this may not be necessary if the object is being created and destroyed due to variable scope
  // Keep for safety
  clearCollectionData(message);
  int sign = (dir < 0)?-1:1;
  bool eastwest = (abs(dir) == 1);
  bool northsouth = (abs(dir) == 2);
  int dir_idx = abs(dir)-1;
  if(abs(dir) != 1 && abs(dir) != 2)
    std::cerr<<"Error, direction must be 1 or 2"<<std::endl;

  // Do not delete grains until message construction is over, to avoid shifting index issues
  std::vector<int> idx_to_del;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++)
  {
    // Cache quantities
    std::vector<double> q = {m_q[grain_num*2],m_q[grain_num*2+1]};
    double r = m_r[grain_num];

    if(sign == -1) {
      // First check if grain has left the boundaries
      if(q[dir_idx] < m_region_min[dir_idx]) {
        idx_to_del.push_back(grain_num);
        addGrainToCollection(message, grain_num);
      } else if (q[dir_idx] - r < m_region_min[dir_idx] + cutoff) {
        addGrainToCollection(message, grain_num);
      }
    } else if (sign == 1 ) {
      // First check if grain has left the boundaries
      if(q[dir_idx] >= m_region_max[dir_idx]) {
        idx_to_del.push_back(grain_num);
        addGrainToCollection(message, grain_num);
      } else if (q[dir_idx] + r > m_region_max[dir_idx] - cutoff) {
        addGrainToCollection(message, grain_num);
      }
    }
  }


  if(dir_idx > 0) {
    // Also go through surroundings for y and z-dirs
    for(int grain_num = 0; grain_num < m_surrounding_collection.num_grains_in_collection; grain_num++)
    {
      // Cache quantities
      std::vector<double> q = {m_surrounding_collection.q[grain_num*2],m_surrounding_collection.q[grain_num*2+1]};
      double r = m_surrounding_collection.r[grain_num];

      if(sign == -1) {
        if (q[dir_idx] - r < m_region_min[dir_idx] + cutoff) {
          addGrainFromCollectionToCollection(m_surrounding_collection,message, grain_num);
        }
      } else if (sign == 1 ) {
        if (q[dir_idx] + r > m_region_max[dir_idx] - cutoff) {
          addGrainFromCollectionToCollection(m_surrounding_collection,message, grain_num);
        }
      }
    }
  }



  // This only works if grains are in order
  int offset = 0;
  for(int grain_num = 0; grain_num < idx_to_del.size(); grain_num++) {
    removeGrain(idx_to_del[grain_num]-offset);
    offset++;
  }
}

void Region::receiveCommunicationMessage(GrainCollection& in_message)
{
  for(int grain_num = 0; grain_num < in_message.num_grains_in_collection; grain_num++) {
    double x_coord = in_message.q[grain_num*2];
    double y_coord = in_message.q[grain_num*2+1];
    double theta = in_message.theta[grain_num];
    double r = in_message.r[grain_num];
    double vx = in_message.v[grain_num*2];
    double vy = in_message.v[grain_num*2+1];
    double omega = in_message.omega[grain_num];
    double fx = in_message.f[grain_num*2];
    double fy = in_message.f[grain_num*2+1];
    std::vector<int> neighbor_list = in_message.neighbor_list[grain_num];
    std::vector<int> contact_list = in_message.contact_list[grain_num];
    int unique_id = in_message.unique_id[grain_num];

    // Check if grain is inside the current region
    if( (x_coord >= m_region_min[0] || m_is_edge[0]) && (x_coord < m_region_max[0] || m_is_edge[1]) &&
        (y_coord >= m_region_min[1] || m_is_edge[2]) && (y_coord < m_region_max[1] || m_is_edge[3]) ) {
      addGrainToRegion(x_coord,y_coord,theta,r,vx,vy,omega,fx,fy,neighbor_list,contact_list,unique_id);
    } else { // Add to surrounding grains
      addGrainFromCollectionToCollection(in_message,m_surrounding_collection,grain_num);
    }
  }
}
