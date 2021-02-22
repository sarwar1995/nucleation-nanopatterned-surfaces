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

#endif /* VolumeSA_calculations_hpp */
