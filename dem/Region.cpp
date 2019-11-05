//
// Created by maytee on 10/23/19.
//

#include <unordered_set>
#include <set>
#include <cfloat>
#include <algorithm>
#include "Region.h"

int Region::m_total_num_grains = 0;
int Region::m_id_tracker = 0;
bool verbose = false;

Region::Region(std::vector<double> in_min, std::vector<double> in_max, int id)
{
  m_region_min = in_min;
  m_region_max = in_max;
  m_num_grains = 0;
  double r_mean = 0.05;
  m_bin_size = r_mean*2.0*2.0;
  m_num_bins_x = static_cast<int>((m_region_max[0] - m_region_min[0])/m_bin_size);
  m_num_bins_y = static_cast<int>((m_region_max[1] - m_region_min[1])/m_bin_size);
  m_region_id = id;
  std::cout<<"binsx: "<<m_num_bins_x<<" binsy: "<<m_num_bins_y<<std::endl;
}

Region::Region()
{
  m_region_min.push_back(0.0); m_region_min.push_back(0.0);
  m_region_max.push_back(2.0); m_region_max.push_back(2.0);
  m_num_grains = 0;
  double r_mean = 0.05;
  m_bin_size = r_mean*2.0*3.0;
  m_num_bins_x = static_cast<int>((m_region_max[0] - m_region_min[0])/m_bin_size);
  m_num_bins_y = static_cast<int>((m_region_max[1] - m_region_min[1])/m_bin_size);
  m_region_id = 0;
}

void Region::generateRandomInitialPacking(double r_mean, int num_total_add)
{
  double r_min = r_mean*0.999, r_max = r_mean*1.001;
  std::cout<<"r range: "<<r_min<<" "<<r_mean<<" "<<r_max<<std::endl;
  double x_in = r_max + m_region_min[0];
  double y_in = r_max + m_region_min[1];
  double v_max = r_min*30.0, v_min = r_min*(-30.0);
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
    addGrainToRegion(x_in,y_in,r_in,v_x,v_y);
    x_in += r_max*2.0;
  }
}

void Region::addGrainToRegion(double x_in, double y_in, double r_in, double vx_in, double vy_in)
{
  m_q.push_back(x_in); m_q.push_back(y_in);
  m_r.push_back(r_in);
  m_v.push_back(vx_in); m_v.push_back(vy_in);
  m_f.push_back(0.0); m_f.push_back(0.0);
  std::vector<int> dummy;
  m_neighbor_list.push_back(dummy);
  m_contact_list.push_back(dummy);
  m_unique_id.push_back(m_id_tracker);
  m_total_num_grains++;
  m_id_tracker++;
  m_num_grains++;
}

void Region::addGrainToRegion(double x_in, double y_in, double r_in, double vx_in, double vy_in, double fx_in, double fy_in,
                              std::vector<int> neighbor_list_in, std::vector<int> contact_list_in, int unique_id)
{
  m_q.push_back(x_in); m_q.push_back(y_in);
  m_r.push_back(r_in);
  m_v.push_back(vx_in); m_v.push_back(vy_in);
  m_f.push_back(fx_in); m_f.push_back(fy_in);
  m_neighbor_list.push_back(neighbor_list_in);
  m_contact_list.push_back(contact_list_in);
  m_unique_id.push_back(unique_id);
  //m_total_num_grains++;
  // ID tracker doesn't change as it is keeping track across all Regions
  //m_id_tracker++;
  m_num_grains++;
}

void Region::removeGrain(int idx)
{
  std::cout<<"Remove "<<idx<<" "<<m_unique_id[idx]<<" from region "<<m_region_id<<std::endl;
  m_q.erase(m_q.begin()+idx*2,m_q.begin()+idx*2+1+1);
  m_r.erase(m_r.begin()+idx);
  m_v.erase(m_v.begin()+idx*2,m_v.begin()+idx*2+1+1);
  m_f.erase(m_f.begin()+idx*2,m_f.begin()+idx*2+1+1);
  m_neighbor_list.erase(m_neighbor_list.begin()+idx);
  m_contact_list.erase(m_contact_list.begin()+idx);
  m_unique_id.erase(m_unique_id.begin()+idx);
  // WARNING, BECAUSE M_TOTAL_NUM_GRAINS IS STATIC ACROSS ALL MEMBERS THIS MAY NOT BE THREAD SAFE
  m_total_num_grains--;
  m_num_grains--;
}

void Region::addGrainToCollection(GrainCollection& collection, int idx)
{
  collection.q.push_back(m_q[idx*2]); collection.q.push_back(m_q[idx*2+1]);
  collection.r.push_back(m_r[idx]);
  collection.v.push_back(m_v[idx*2]); collection.v.push_back(m_v[idx*2+1]);
  collection.f.push_back(m_f[idx*2]); collection.f.push_back(m_f[idx*2+1]);
  collection.neighbor_list.push_back(m_neighbor_list[idx]);
  collection.contact_list.push_back(m_contact_list[idx]);
  collection.unique_id.push_back(m_unique_id[idx]);
  collection.region_id.push_back(m_region_id);
  collection.num_grains_in_collection++;
}

void Region::addGrainFromCollectionToCollection(GrainCollection &collection_in, GrainCollection &collection_receive, int grain_idx)
{
  collection_receive.q.push_back(collection_in.q[grain_idx*2]); collection_receive.q.push_back(collection_in.q[grain_idx*2+1]);
  collection_receive.r.push_back(collection_in.r[grain_idx]);
  collection_receive.v.push_back(collection_in.v[grain_idx*2]); collection_receive.v.push_back(collection_in.v[grain_idx*2+1]);
  collection_receive.f.push_back(collection_in.f[grain_idx*2]); collection_receive.f.push_back(collection_in.f[grain_idx*2+1]);
  collection_receive.neighbor_list.push_back(collection_in.neighbor_list[grain_idx]);
  collection_receive.contact_list.push_back(collection_in.contact_list[grain_idx]);
  collection_receive.unique_id.push_back(collection_in.unique_id[grain_idx]);
  collection_receive.region_id.push_back(collection_in.region_id[grain_idx]);
  collection_receive.num_grains_in_collection++;
}


void Region::domainDecomposition(int num_regions, std::vector<Region>& region_list)
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
  //max_y = m_region_max[1];
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

  for(int y_regions = 0 ; y_regions < divisors[min_idx].second; y_regions++) {
    for(int x_regions = 0 ; x_regions < divisors[min_idx].first; x_regions++) {
      std::vector<double> region_min = {min_x + x_regions*dx_final,min_y+y_regions*dy_final};
      std::vector<double> region_max = {min_x + (x_regions+1)*dx_final,min_y+(y_regions+1)*dy_final};
      int region_idx = x_regions + y_regions*divisors[min_idx].first;
      std::cout<<"region: "<<region_idx<<" has "<<grain_indexes_in_regions[region_idx].size()<<" grains "<<std::endl;
      std::cout<<"      minx: "<<region_min[0]<<" miny: "<<region_min[1]<<std::endl;
      std::cout<<"      maxx: "<<region_max[0]<<" maxy: "<<region_max[1]<<std::endl;
      Region region(region_min,region_max,region_idx);

      for(int ii = 0; ii < grain_indexes_in_regions[region_idx].size(); ii++) {
        int idx = grain_indexes_in_regions[region_idx][ii];
        region.addGrainToRegion(m_q[idx*2],m_q[idx*2+1],m_r[idx],m_v[idx*2],m_v[idx*2+1],m_f[idx*2],m_f[idx*2+1],
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
  m_bin_list.clear();
  m_bin_list.resize(m_num_bins_x*m_num_bins_y);

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
  check_set.resize(m_num_grains+m_surrounding_collection.num_grains_in_collection);
  // Loop through all bins to see what grains intersect them
  for(int bin_idx = 0; bin_idx < m_num_bins_x*m_num_bins_y; bin_idx++){
    // Take the first grain in the bin as the one that collects the contacts
    // This assumes the grains are in ascending ID order in the bins
    for(int grain_num_1 = 0; grain_num_1 < m_bin_list[bin_idx].size(); grain_num_1++) {
      for(int grain_num_2 = grain_num_1+1; grain_num_2 < m_bin_list[bin_idx].size(); grain_num_2++){
        check_set[m_bin_list[bin_idx][grain_num_1]].insert(m_bin_list[bin_idx][grain_num_2]);
      //m_neighbor_list[m_bin_list[bin_idx][0]].push_back(m_bin_list[bin_idx][grain_num]);
      }
    }

  }
  for(int ii = 0; ii < m_num_grains; ii++) {
    m_neighbor_list[ii].clear();
    if(check_set[ii].size() > 0) {
      std::vector<int> list(check_set[ii].begin(),check_set[ii].end());
      m_neighbor_list[ii] = list;
    }
  }

}

// Because surrounding grains will always have larger indexes (not unique IDs) than the grains in the region,
// they will always show up in the neighbor lists of the grains in the region, and thus do not need to loop over
// surrounding grains
void Region::buildContactList()
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

  int actual_contacts = 0;
  for(int grain_num = 0; grain_num < m_num_grains; grain_num++){
    std::vector<int> contact_list;
    m_contact_list[grain_num].clear();
    for(int check_num = 0; check_num < m_neighbor_list[grain_num].size(); check_num++){
      int check_id = m_neighbor_list[grain_num][check_num];
      int surround_id = check_id - m_num_grains;

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
      }

      if(dist_squared <= r_squared){
        contact_list.push_back(check_id);
        actual_contacts++;
        if(verbose) {
          if(check_id < m_num_grains) {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" x1 "<<m_q[grain_num*2]<<" x2 "<<m_q[check_id*2];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_surrounding_collection.r[check_id]<<std::endl;
          }
          else {
            std::cout<<"grain: "<<grain_num<<"/"<<m_num_grains<<" check neighbor: "<<check_id<<" y1 "<<m_q[grain_num*2+1]<<" y2 "<<m_surrounding_collection.q[surround_id*2+1];
            std::cout<<" r1 "<<m_r[grain_num]<<" r2 "<<m_surrounding_collection.r[surround_id]<<std::endl;
          }
        }
      }
    }
    m_contact_list[grain_num] = contact_list;
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

  std::cout<<"#actual contacts: "<<actual_contacts<<std::endl;
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

        double overlap = sqrt(r_squared - dist_squared);
        std::vector<double> n{0.0,0.0};
        n[0] = (m_q[check_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_q[check_id*2+1] - m_q[grain_num*2+1])/dist;

        m_f[grain_num*2] -= n[0]*overlap*k;
        m_f[grain_num*2+1] -= n[1]*overlap*k;

        m_f[check_id*2] += n[0]*overlap*k;
        m_f[check_id*2+1] += n[1]*overlap*k;
        if(verbose) {
          std::cout<<"F btwn "<<grain_num<<","<<check_id<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<" f1x "<<m_f[grain_num*2]<<" f1y "<<m_f[grain_num*2+1]<<std::endl;
          std::cout<<" f2x "<<m_f[check_id*2]<<" f2y "<<m_f[check_id*2+1]<<std::endl;
        }


      } else {
        int surround_id = check_id - m_num_grains;
        double dist_squared = (m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])*(m_q[grain_num*2]-m_surrounding_collection.q[surround_id*2])
                          + (m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1])*(m_q[grain_num*2+1]-m_surrounding_collection.q[surround_id*2+1]);
        double r_squared = (m_r[grain_num] + m_surrounding_collection.r[surround_id]) * (m_r[grain_num] + m_surrounding_collection.r[surround_id]);

        double dist = sqrt(dist_squared);

        double overlap = sqrt(r_squared - dist_squared);
        std::vector<double> n{0.0,0.0};
        n[0] = (m_surrounding_collection.q[surround_id*2] - m_q[grain_num*2])/dist;
        n[1] = (m_surrounding_collection.q[surround_id*2+1] - m_q[grain_num*2+1])/dist;

        if(verbose)
          std::cout<<"F surround btwn "<<grain_num<<","<<surround_id<<" overlap "<<overlap<<" n0 "<<n[0]<<" n1 "<<n[1]<<std::endl;

        // Only apply forces on grain in the region
        m_f[grain_num*2] -= n[0]*overlap*k;
        m_f[grain_num*2+1] -= n[1]*overlap*k;
      }

    }
    // Check against floor
    if(m_q[grain_num*2+1]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2+1] += (m_r[grain_num]-m_q[grain_num*2+1])*k;
    }
    // Check against left side
    if(m_q[grain_num*2]-m_r[grain_num] <= 0.0)
    {
      m_f[grain_num*2] += (m_r[grain_num]-m_q[grain_num*2])*k;
    }
    // Check against right side
    if(m_q[grain_num*2]+m_r[grain_num] >= 4.0)
    {
      m_f[grain_num*2] += (4.0-m_q[grain_num*2]-m_r[grain_num])*k;
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
    double r = in_message.r[grain_num];
    double vx = in_message.v[grain_num*2];
    double vy = in_message.v[grain_num*2+1];
    double fx = in_message.f[grain_num*2];
    double fy = in_message.f[grain_num*2+1];
    std::vector<int> neighbor_list = in_message.neighbor_list[grain_num];
    std::vector<int> contact_list = in_message.contact_list[grain_num];
    int unique_id = in_message.unique_id[grain_num];

    // Check if grain is inside the current region
    if(x_coord >= m_region_min[0] && x_coord < m_region_max[0] && y_coord >= m_region_min[1] && y_coord < m_region_max[1]) {
      addGrainToRegion(x_coord,y_coord,r,vx,vy,fx,fy,neighbor_list,contact_list,unique_id);
    } else { // Add to surrounding grains
      addGrainFromCollectionToCollection(in_message,m_surrounding_collection,grain_num);
    }
  }
}
