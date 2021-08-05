//
//  periodic_io.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 6/23/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef periodic_io_hpp
#define periodic_io_hpp

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <cmath>
#include "parallel_process.hpp"
#include <mpi.h>
#include "VolumeSA_calculations.hpp"

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
    int n_max_points;
};

class PeriodicIO
{
public:
    PeriodicIO();
    ~PeriodicIO();
    PeriodicIO(ParallelProcess*, SphericalCapOutput*, std::string, int); //Add profiling vectors
    
    int gather_and_write();
    int write_individually();
    int write_profiling_individually(std::vector<std::vector<double>>&, std::vector<int>&);
    
protected:
    int myRank, nProcs;
    ParallelProcess* parallel_process;
    
    SphericalCapOutput* output_variables;
    SphericalCapOutput gathered_output_variables;
    int num_patches;
    
    FILE* V_SA_DataFile;
    std::string V_SA_DataFileName;
    
    FILE* V_SA_DataFile_per_root_rank;
    std::string V_SA_DataFileName_per_root_rank;
    
    FILE* profilingFile_per_root_rank;
    std::string profiling_FileName_per_root_rank;
    
    //MPI gathering
    int gather ();
    bool gathered_output_variables_isempty();
    
    //output file
    std::string tag;
    void get_fileName();
    void get_fileName_per_root_rank();
    void get_profiling_fileName_per_root_rank();
    int open_output_file_for_appending();
    int print_to_file();
    
    int conditinoal_open_per_root_rank_output_file_for_appending();
    int print_per_root_rank_to_file();
    int conditinoal_open_per_root_rank_profiling_file_for_appending();
    int print_profiling_per_root(std::vector<std::vector<double>>&, std::vector<int>&);
    
    /* Size variables related to gathering*/
    int clstr_centre_modifier_size;
    int V_SA_size;
    int proj_SA_size;
    int radii_size;
};


#endif /* periodic_io_hpp */
