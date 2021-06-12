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
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "parallel_process.hpp"
#include "Stripes.hpp"
#include "MC_parallel.hpp"
#include "surface_setup.hpp"
#include "DynamicBox.hpp"
#include "CheckBoundary.hpp"
#include "FreeEnergy.hpp"


enum CapType{
    CENTRE_GOOD,
    GOOD,
    BAD
};

struct SphericalCapOutput
{
    std::vector<int> clstr_centre_location_modifier_global_array;
    std::vector<double> radii_global_array;
    std::vector<double> Number_particles_global_array;
    std::vector<double> Volume_global_array;
    std::vector<double> SA_global_array;
    std::vector<double> projected_SA_global_array;
};

struct SphericalCapInput
{
    double Rg_max, Rb_max, d_Rg, d_Rb;
    double theta_good;
    double theta_bad;
    double delta;
    double Rho;
};

class EvolveSphericalCap
{
public:
    EvolveSphericalCap();
    ~EvolveSphericalCap();
    EvolveSphericalCap(ParallelProcess*, Stripes, SphericalCapOutput* , SphericalCapInput*, CheckBoundary*, MC* , int);
    
    void init_cap_identifier();
    
    int evolve ();
    int dummy_evolve(int);
    void dummy_set_next_evolve();
    
    //Calculating Volume and SA
    void calc_volume_SA_central_good_cap(Spherical_cap&);
    void calc_volume_SA(Composite_cluster&);
    
    //Bounds checking
    bool patch_bounds_crossed();
    bool check_breaking_condition();
    
    //printing
    void print_initial_variables();
    void print_growing_cap_identifier();
    
    
    bool has_one_cap();
    bool has_no_caps();
    bool has_more_than_one_caps();
    

    //Setting MC_parallel communicator
    void set_MC_comm ();
    
    //printing surface points
    
protected:
    int myRank, nProcs;
    /* Member classes */
    ParallelProcess* parallel_process;
    Stripes stripes;
    
    Shape* Cluster_shape_ptr;
    CheckBoundary* check_boundary;
    MC* mc_engine;
    
    int parallelisation_lvl; // [Rg: 0, dB:1, Rb:2]
    
    /* Quantity variables */
    int num_patches;
    double Radius;
    double Rg_max, Rb_max, d_Rg, d_Rb;
    double theta_good;
    double theta_bad;
    double delta;
    double Rho;
    std::vector<double> centre_left_cap;
    std::vector<double> centre_right_cap;
    std::vector<double> mc_volume_SA;
    
    
    /* Cap-related variables */
    CapType cap_type;
    
    //An array of 0,1 where the current growing cap index and all indices for previous (existing) caps is 1.
    // Ex: If the second pair of good caps is being grown on a surface of 7 patches i.e. B2|G1|B1|G0|B1|G1|B2
    //then the array will be [1, 1, 1, 0] ===> [central_good, next pair of bad, next pair of good ....] where
    //the latest 1 represents the growing caps.
    std::vector<int> current_growing_cap_identifier;
    int growing_cap_index;
    
    /* Current and previous caps */
    //centres of all caps corresponding to the current configuration
    //of the composite cluster. [centre good, bad (left), good (left), ...]
    //radii of all caps corresponding the current configuration of the composite cluster.
    //[centre good, bad (left), good (left), ...]
    std::vector<std::vector<double> > current_cap_centres;
    std::vector<double> current_centre_to_left_patch_boundary_distance;
    std::vector<double> current_projected_radii_array;
    std::vector<double> current_radii_array;
    std::vector<double> current_clstr_centre_modifier_array;
    std::vector<Spherical_cap> growing_capsList;
    
    /* calculate centre to left patch boundary */
    void calc_and_set_centre_to_left_patch_boundary_distance();
    
    /* cap identifier related variables and functions */
    void identify_growing_cap();
    void identify_cap_type();
    bool is_valid_growing_cap_identifier ();
    
    /* current cap growth conditions */
    double previous_cap_radius, previous_cap_projected_radius,  previous_cap_theta;
    double previous_centre_to_left_patch_boundary_distance;
    double current_radius_min, current_radius_max, d_current_radius;
    double current_cap_radius, current_projected_radius;
    double current_theta;
    double current_centre_to_previous_patch_boundary_distance;
    double previous_patch_symmetric_boundary;
    int length_current_radius;
    int clstr_centre_location_modifier;
    
    
    std::vector<int> Radius_loop_start_end;
    
    int identify_current_cap_growth_conditions ();
    
    /* Output variables struct */
    SphericalCapOutput* output_variables;
    
    /* cap evolution related functions */
    void update_growing_cap_identifier();
    void reset_growing_cap_identifier();
    void reset_growing_capsList();
    void reset_current_cap_growth_conditions();
    /* Setting up next cluster evolution*/
    void setup_next_evolve ();
    void reset_evolve ();

    /* Adding quantities to the global arrays*/
    void add_radii_and_centre_modifier();
    
    //Dummy evolve
    int dummy_decider;
};


#endif /* evolve_spherical_caps_hpp */
