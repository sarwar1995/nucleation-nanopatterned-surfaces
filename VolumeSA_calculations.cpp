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
