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
#include "FreeEnergy.hpp"
#include "surface_setup.hpp"
#include "DynamicBox.hpp"
#include "CheckBoundary.hpp"
#include "parallel_process.hpp"
#include "evolve_spherical_caps.hpp"
#include "VolumeSA_calculations.hpp"


class ManagerSphericalCaps {

public:
    ManagerSphericalCaps();
    ~ManagerSphericalCaps();
    ManagerSphericalCaps(int , int); //int, MPI_Comm
    
    //Setup
    void setup(char* argv[], int start_index);
    
    //evolve
    int evolve();
    void dummy_evolve(int);
    
    //MPI gathering
    void gather ();
    
    //Related to output variables
    bool is_empty_output_variables();
    
    //output file
    void open_output_file();
    void print_to_file();
    
    //printing
    void print_quants(int);
    void print_surface_ptr();
    void print_box();
    void print_mc_and_check_boundary();
    void print_output_variables();
    
    //free comms
    void free_MPI_comms();
    
protected:
    int myRank, nProcs;
    ParallelProcess parallel_process;
    EvolveSphericalCap evolve_spherical_cap;
    /*Member classes*/
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
    FILE* SurfacePointsFile;
    std::string tag;
    
    /* Reading input related functions */
    
    //The integer is the starting index,
    //as some inputs will be read before (start_index != 0)
    void read_from_command_line(char* argv[], int start_index);
    void broadcast_data_to_workers();

    /* Setup functions */
    void setup_surface();
    void setup_box();
    void setup_mc_and_boundary();
    void setup_input_variables();
    void setup_evolve_spherical_caps();
    
    
    /*Output variables: These vectors consist of vectors of location_modifiers and radii_arrays
     for each pair of symmetric patches surrounding the good patch and the good patch.*/
    SphericalCapOutput output_variables;
    SphericalCapInput input_variables;
    SphericalCapOutput gathered_output_variables;

    /* Size variables related to gathering*/
    int clstr_centre_modifier_size;
    int V_SA_size;
    int proj_SA_size;
    int radii_size;
};


#endif /* manager_spherical_caps_hpp */
