//
//  evolve_spherical_caps.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "evolve_spherical_caps.hpp"

EvolveSphericalCap::EvolveSphericalCap()
{}

EvolveSphericalCap::~EvolveSphericalCap()
{
    
}

EvolveSphericalCap::EvolveSphericalCap(ParallelProcess process, Stripes stripes):
parallel_process(process), stripes(stripes)
{
    
}

void EvolveSphericalCap::init_cap_identifier()
{
    int n_unique_patches = stripes.get_n_unique_patches();
    current_growing_cap_identifier.resize(n_unique_patches);
    for(int i=0; i<n_unique_patches; i++)
    {
        current_growing_cap_identifier[i] = 0;
    }
    
}


void EvolveSphericalCap::evolve (std::vector<int>& clstr_centre_location_modifier_global_arrays, std::vector<double>& radii_global_arrays, std::vector<double>& Number_particles, std::vector<double>&Volume_array, std::vector<double>&SA_array, std::vector<double>& projected_SA_array)
