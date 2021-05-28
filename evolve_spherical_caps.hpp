//
//  evolve_spherical_caps.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef evolve_spherical_caps_hpp
#define evolve_spherical_caps_hpp

#include <stdio.h>
#include "parallel_process.hpp"
#include "Stripes.hpp"

enum CapType{
    CENTRE_GOOD,
    GOOD,
    BAD
};

struct spherical_cap_output
{
    std::vector<int> clstr_centre_location_modifier_global_arrays;
    std::vector<double> radii_global_arrays;
    std::vector<double> Number_particles;
    std::vector<double> Volume_array;
    std::vector<double> SA_array;
    std::vector<double> projected_SA_array;
}


class EvolveSphericalCap
{
public:
    EvolveSphericalCap();
    ~EvolveSphericalCap();
    EvolveSphericalCap(ParallelProcess, Stripes);
    
    void init_cap_identifier();
    
    void evolve (spherical_cap_output* output_vars);
    
protected:
    /* Member classes */
    ParallelProcess parallel_process;
    Stripes stripes;
    
    
    /* Quantity variables */
    double Volume;
    double SA;
    std::vector<double> projected_SA;
    double N;   //Number of particles per point (i.e. per cluster)
    double Radius;
    double clstr_centre_location_modifier;
    std::vector<double> centre_left_cap(3,0.0);
    std::vector<double> centre_right_cap (3,0.0);
    std::vector<double> mc_volume_SA(2,0.0);
    
    /* Cap-related variables */
    CapType cap_type;
    
    /* Surface-related variables */
    std::vector<int> stripes_bounds;
    
    
    //An array of 0,1 where the current growing cap index and all indices for previous (existing) caps is 1.
    // Ex: If the second pair of good caps is being grown on a surface of 7 patches i.e. B2|G1|B1|G0|B1|G1|B2
    //then the array will be [1, 1, 1, 0] ===> [central_good, next pair of bad, next pair of good ....] where
    //the latest 1 represents the growing caps.
    std::vector<int> current_growing_cap_identifier;
    
    /* Current and previous caps */
    //centres of all caps corresponding to the current configuration
    //of the composite cluster. [centre good, bad (left), good (left), ...]
    std::vector<std::vector<double> > cap_centres;
    
    //radii of all caps corresponding the current configuration of the composite cluster.
    //[centre good, bad (left), good (left), ...]
    std::vector<double> cap_projected_radii;
    
    
};


#endif /* evolve_spherical_caps_hpp */
