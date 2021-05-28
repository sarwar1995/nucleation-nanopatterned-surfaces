//
//  manager_spherical_caps.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/22/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef manager_spherical_caps_hpp
#define manager_spherical_caps_hpp

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <cmath>
//#include "Composite_cluster.hpp"
//#include "Surface/Surface.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC_parallel.hpp"
#include "VolumeSA_calculations.hpp"
#include "FreeEnergy.hpp"
#include "surface_setup.hpp"
#include "DynamicBox.hpp"
#include "CheckBoundary.hpp"
#include "parallel_process.hpp"
#include "evolve_spherical_caps.hpp"
#include <chrono>


class ManagerSphericalCaps {

public:
    ManagerSphericalCaps();
    ~ManagerSphericalCaps();
    ManagerSphericalCaps(int, int , int); //int, MPI_Comm
    
    //Setup
    void setup(char* argv[], int start_index);
    
    
    //cap evolution
    //Recursively evolve spherical caps based on whether boundaries have been crossed.
    void evolve_cap ();
    void evolve_three_caps();
    void evolve_bad_cap(double);
    
    
    //Bounds checking
    bool check_bounds();
    bool check_breaking_condition();
    
    //Calculating Volume and SA
    void calc_volume_SA();
    
    //MPI gathering
    void gather ();
    void print_nelements();
    
protected:
    
    ParallelProcess parallel_process;
    EvolveSphericalCap evolve_spherical_cap;
    /*Member classes*/
    Shape* Cluster_shape_ptr;
    CheckBoundary check_boundary;
    MC mc_engine;
    DynamicBox dynamic_box;
    Stripes stripes;
    Surface* surface_ptr;
    
    
    /* Input variables */
    double Rg_max, Rb_max, d_Rg, d_Rb;
    int num_patches;
    double z_surface;
    double good_patch_x_width, good_patch_y_width, bad_patch_x_width, bad_patch_y_width;
    double theta_good;
    double theta_bad;
    double starting_box_dim;
    double extension_length;
    /* MC variables */
    int n_points;
    double point_density;
    int mc_seed[3];
    double delta;
    /* physical variables */
    double Rho ;
    /*IO variables */
    FILE* V_SA_DataFile;
    std::string tag;
    
    /* Reading input related functions */
    
    //The integer is the starting index,
    //as some inputs will be read before (start_index != 0)
    void read_from_command_line(char* argv[], int start_index);
    void broadcast_data_to_workers();

    /* Setup functions */
    void init_cap_identifier();
    void setup_surface();
    void setup_box();
    void setup_mc_and_boundary();
    void setup_output_variables();
    
    
    /*Output variables: These vectors consist of vectors of location_modifiers and radii_arrays
     for each pair of symmetric patches surrounding the good patch and the good patch.*/
    spherical_cap_output output_variables;


    //Gathering data from all processes
    int* counts;

};


#endif /* manager_spherical_caps_hpp */
