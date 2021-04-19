//
//  VolumeSA_calculations.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/17/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef VolumeSA_calculations_hpp
#define VolumeSA_calculations_hpp

#include <stdio.h>
#include <vector>
#include "FreeEnergy.hpp"
#include "miscformulas.hpp"

void add_Volume_SA (std::vector<double> Radii, std::vector<double> mc_volume_SA, std::vector<double>compcluster_projected_SA, double Rho, int db, int dg_secondary, FILE* VolumeSAouput);

void add_Volume_SA_parallel(std::vector<double> Nparticles_global_array, std::vector<double> Radii_global_array, std::vector<double> Volume_global_array, std::vector<double>SA_global_array, std::vector<double>proj_SA_global_array, std::vector<int> db_global_array, std::vector<int> dg_secondary_global_array , FILE* VolumeSAouput);

void add_Volume_SA_spherocylinder_parallel(std::vector<int> dB_global_array, std::vector<double> Number_particles_global, std::vector<double> Volume_global_array, std::vector<double> SA_global_array, std::vector<double> projected_global_SA_array, std::vector<double> Rg_global_array, std::vector<double> Rb_global_array, std::vector<double> cyl_length_global_array, std::vector<double> chord_length_global_array, FILE* outputfile);

#endif /* VolumeSA_calculations_hpp */
