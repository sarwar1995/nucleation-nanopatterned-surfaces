//
//  VolumeSA_calculations.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/17/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "VolumeSA_calculations.hpp"

void add_Volume_SA (std::vector<double> Radii, std::vector<double> mc_volume_SA, std::vector<double>compcluster_projected_SA, double Rho, int db, int dg_secondary, FILE* VolumeSAouput)
{
    /* Radii = [R_good, R_bad, R_good_secondary]
     VolumeSAouput = [N(not rounded), Volume_MC, SA_MC, projected_SA (vector), Radii (vector), db (+1, -1), dg_secondary (+1,-1)]
     */
    
    /*
     The volume and surface areas provided to this function are printed out in angstrom units
     */
    
    double Volume_MC, SA_MC, N;
    Volume_MC = mc_volume_SA[0];
    SA_MC = mc_volume_SA[1];
    N = (Volume_MC*1e-30) * Rho * avogadro ; //This volume needs to be in m3 for this to be valid since Rho is in moles/m3
    //Not printing out rNhere because this conversion will be automatically made by the add_to_N function  in the Get_FreeEnergyProfile.cpp file
    
    fprintf(VolumeSAouput,"%10.10f\t%10.20f\t%10.20f\t", N, Volume_MC, SA_MC);
    for(size_t i=0; i<compcluster_projected_SA.size(); i++)
    {
        fprintf(VolumeSAouput, "%10.20f\t", compcluster_projected_SA[i]);
    }
    for(size_t i=0; i<Radii.size(); i++)
    {
        if(i==Radii.size()-1)
        {
            fprintf(VolumeSAouput, "%10.20f\t", Radii[i]);
        }
        else
        {
            fprintf(VolumeSAouput, "%10.20f\t", Radii[i]);
        }
    }
    
    fprintf(VolumeSAouput, "%d\t%d\n", db, dg_secondary);

}

//The one in use for the manager spherical caps class
void add_Volume_SA_parallel(std::vector<double> Nparticles_global_array, std::vector<double> Radii_global_array, std::vector<double> Volume_global_array, std::vector<double>SA_global_array, std::vector<double>proj_SA_global_array, std::vector<int> clstr_centre_location_modifier,  int n_patches, FILE* VolumeSAouput)
{
     int n_unique_patches = ((n_patches-1)/2) + 1;
     int vol_size = (int)Volume_global_array.size();
     bool volume_size_comparison = (Nparticles_global_array.size() == Volume_global_array.size() && Volume_global_array.size() == SA_global_array.size()) ;
    
    bool clstr_identifier_size_comparison = (Radii_global_array.size() == clstr_centre_location_modifier.size() && (int)Radii_global_array.size() == (vol_size * n_unique_patches));
    
    bool proj_SA_size_comparison = ((int)proj_SA_global_array.size() == n_patches*vol_size) ;
    
    if(volume_size_comparison && clstr_identifier_size_comparison && proj_SA_size_comparison)
    {
        for(int i=0; i<vol_size; i++)
        {
            fprintf(VolumeSAouput,"%10.10f\t%10.10f\t%10.10f\t", Nparticles_global_array[i], Volume_global_array[i], SA_global_array[i]);
            int k_start = i * n_patches;
            int k_end = k_start + n_patches;
            for(int k=k_start; k<k_end; k++)
            {
                fprintf(VolumeSAouput, "%10.10f\t", proj_SA_global_array[k]);
            }
            int j_start = (int)i * n_unique_patches;
            int j_end = j_start + n_unique_patches;
            for(int j=j_start; j<j_end; j++)
            {
                fprintf(VolumeSAouput, "%d\t", clstr_centre_location_modifier[j]);
            }
            for(int j=j_start; j<j_end; j++)
            {
                fprintf(VolumeSAouput, "%10.10f\t", Radii_global_array[j]);
            }
            
            
            fprintf(VolumeSAouput, "\n");
        }
    }
    else
    {
        printf("Sizes of global arrays are not correct. volume_size_comparison=%d\t clstr_identifier_size_comparison=%d\t proj_SA_size_comparison=%d\n", volume_size_comparison, clstr_identifier_size_comparison, proj_SA_size_comparison);
        printf("n_particles=%d volume=%d SA=%d radii=%d proj_SA=%d clstr_centre=%d\n", (int)Nparticles_global_array.size(), (int)Volume_global_array.size(), (int)SA_global_array.size(), (int)Radii_global_array.size(), (int)proj_SA_global_array.size(), (int)clstr_centre_location_modifier.size());
        
        abort();
    }
    
    
}

void add_Volume_SA_parallel(std::vector<double> Nparticles_global_array, std::vector<double> Radii_global_array, std::vector<double> Volume_global_array, std::vector<double>SA_global_array, std::vector<double> proj_SA_global_array, std::vector<int> db_global_array, std::vector<int> dg_secondary_global_array , FILE* VolumeSAouput)
{
    int size = (int)Volume_global_array.size();
    
    printf("Volume size=%d\t proj_sa_size=%d\t radii_size=%d\n", size, (int)proj_SA_global_array.size(), (int)Radii_global_array.size());
    if(!((int)SA_global_array.size() == size && (int)db_global_array.size() == size && (int)dg_secondary_global_array.size() == size && (int)Nparticles_global_array.size() == size))
    {
        printf("Sizes of one of N=%d, SA=%d, db=%d or dg=%d is not the same as V global array",(int)Nparticles_global_array.size(), (int)SA_global_array.size(), (int)db_global_array.size(), (int)dg_secondary_global_array.size());
        abort();
    }
    else
    {
        for(size_t i=0; i<Volume_global_array.size(); i++)
        {
            fprintf(VolumeSAouput,"%10.10f\t%10.20f\t%10.20f\t", Nparticles_global_array[i], Volume_global_array[i], SA_global_array[i]);
            
            int k_start = (int)i * 5;
            int k_end = k_start+5;
            for(int k=k_start; k<k_end; k++)
            {
                fprintf(VolumeSAouput, "%10.20f\t", proj_SA_global_array[k]);
            }
            
            int j_start = (int)i * 3;
            int j_end = j_start+3;
            for(int j=j_start; j<j_end; j++)
            {
                fprintf(VolumeSAouput, "%10.20f\t", Radii_global_array[j]);
                
            }
            fprintf(VolumeSAouput, "%d\t%d\n", db_global_array[i], dg_secondary_global_array[i]);

        }
    }
    
}

void add_Volume_SA_spherocylinder_parallel(std::vector<int> dB_global_array, std::vector<double> Number_particles_global, std::vector<double> Volume_global_array, std::vector<double> SA_global_array, std::vector<double> projected_global_SA_array, std::vector<double> Rg_global_array, std::vector<double> Rb_global_array, std::vector<double> cyl_length_global_array, std::vector<double> chord_length_global_array, FILE* outputfile)
{
    int size = (int)Volume_global_array.size();
    
    if(!((int)SA_global_array.size() == size && (int)dB_global_array.size() == size && (int)Number_particles_global.size() == size))
    {
        printf("Sizes of one of N=%d, SA=%d, db=%d is not the same as V global array",(int)Number_particles_global.size(), (int)SA_global_array.size(), (int)dB_global_array.size());
        abort();
    }
    else
    {
        for(size_t i=0; i<Volume_global_array.size(); i++)
        {
            fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%d\t%10.10f\t%10.10f\t%10.10f\t" , Rg_global_array[i], cyl_length_global_array[i], chord_length_global_array[i], Rb_global_array[i], dB_global_array[i], Number_particles_global[i], Volume_global_array[i], SA_global_array[i]);
            
            int k_start = (int)i * 3;
            int k_end = k_start+3;
            for(int k=k_start; k<k_end; k++)
            {
                fprintf(outputfile, "%10.10f\t", projected_global_SA_array[k]);
            }
            fprintf(outputfile, "\n");
        }
        
    }
  
}


