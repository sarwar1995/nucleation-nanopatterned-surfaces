//
//  Get_FreeEnergyProfile.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/16/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include "FreeEnergy.hpp"
#include "miscformulas.hpp"
#include "Get_FreeEnergyProfile.hpp"

GetFreeEnergyProfile::GetFreeEnergyProfile(FILE* input, FILE* output, int min_N, int max_N, double d_N, double rho, double mu, double sigma, double t, double theta_g, double theta_b)
: Nmin(min_N), Nmax (max_N), dN(d_N), input_file(input), output_file(output), Rho(rho), Mu(mu), Sigma(sigma), T(t), theta_cg(theta_g), theta_cb(theta_b)
{
    len_N  = (int) ((Nmax-Nmin)/dN) + 1;
    NArray_Gmin.resize(len_N, 0.0);    //Vector of minimum G value for each integer N value
    NArray_confs.resize (len_N , 0.0);   //Number of distinct points (i.e. clusters) for each integer value of N
    NArray_quant.resize (len_N); //[N, Gmin, Rg, Rb, Rg_secondary, db, dg_secondary, V, SA]. Quantities such as G, Rg, Rb etc. correponding to the minimum G cluster for each integer value of N
}

GetFreeEnergyProfile::~GetFreeEnergyProfile()
{
    
}


void GetFreeEnergyProfile::ReadProfileSphericalCaps()
{
    double N, Volume, SA;
    double projected_SA_single_cap ;
    int db, dg_secondary;
    
    std::vector<double> projected_SA(5, 0.0);   //This is obviously only for the 5 patch case with 3 pseudo-indipendent radii.
    std::vector<double> Radii (3, 0.0);
    /* projected_SA = [centre_good, bad_left, bad_right, good_left, good_right]
     Radii        = [centre_good, bad, secondary_good]*/
    double G;
    
    //NArray_quant will have 9 values only for the all spherical caps case
    for (int i=0 ; i<len_N; i++)
    {
        NArray_quant[i] = std::vector<double> (9,0.0);
        if(i==0)
        {
            NArray_quant[i][0] = Nmin;
        }
        else
        {
            NArray_quant[i][0]=NArray_quant[i-1][0]+dN ;
        }
    }
    
    format = "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf%d\t%d\n" ;
    if(input_file == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    int fscanf_result = fscanf(input_file, format, &N, &Volume, &SA, &projected_SA[0], &projected_SA[1], &projected_SA[2], &projected_SA[3], &projected_SA[4], &Radii[0], &Radii[1], &Radii[2], &db, &dg_secondary);
    
    while(fscanf_result != EOF)
    {
        Volume = Volume * 1e-30;    //This converts the volume from Ang^3 to m^3
        SA = SA* 1e-20;
        for(size_t y=0; y<projected_SA.size(); y++)
        {
            projected_SA[y] = projected_SA[y]*1e-20 ;
        }
        
        if(Radii[1]==0.0 && Radii[2] == 0.0) //Just a single good cap case
        {
            projected_SA_single_cap = projected_SA[0];
            if(projected_SA[1] != 0.0 || projected_SA[2] != 0.0)
            {
                printf("projected SA of bad and good secondary are not zero for Rb=Rg_second=0.0, projected_SA[1]=%10.10f\t projected_SA[2]=%10.10f\n", projected_SA[1], projected_SA[2]);
                abort();
            }
            G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, projected_SA_single_cap, theta_cg, T);
            add_to_N(N,  G, Radii[0], Radii[1], Radii[2], db, dg_secondary, Volume, SA, dN, Nmin, len_N, NArray_Gmin, NArray_confs,  NArray_quant);
        }
        
        else
        {
            G = free_energy(Rho, Mu, Sigma, Volume, SA, projected_SA, theta_cg, theta_cb, T);
            add_to_N(N,  G, Radii[0], Radii[1], Radii[2], db, dg_secondary, Volume, SA, dN, Nmin, len_N, NArray_Gmin, NArray_confs,  NArray_quant);
        }
        fscanf_result = fscanf(input_file, format, &N, &Volume, &SA, &projected_SA[0], &projected_SA[1], &projected_SA[2], &projected_SA[3], &projected_SA[4], &Radii[0], &Radii[1], &Radii[2], &db, &dg_secondary);
    }
    print_NGDataFile (NArray_quant, output_file);
    fclose(output_file);
    fclose(input_file);
}

void GetFreeEnergyProfile::ReadProfileSpherocylinder()
{
    //Input file should be in the format:
    // Rg, cyl_length, N, Volume, SA, ProjSA
    format = "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n";
    double Rg, cyl_length, N, Volume, SA, ProjSA;
    double G;
    //NArray_quant will have 9 values only for the all spherical caps case
    for (int i=0 ; i<len_N; i++)
    {
        NArray_quant[i] = std::vector<double> (7,0.0);
        if(i==0)
        {
            NArray_quant[i][0] = Nmin;
        }
        else
        {
            NArray_quant[i][0]=NArray_quant[i-1][0]+dN ;
        }
    }
    if(input_file == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    int fscanf_result = fscanf(input_file, format, &Rg, &cyl_length, &N, &Volume, &SA, &ProjSA);
    
    while(fscanf_result != EOF)
    {
        Volume = Volume * 1e-30;    //This converts the volume from Ang^3 to m^3
        SA = SA * 1e-20;
        ProjSA = ProjSA * 1e-20;
        G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, ProjSA, theta_cg, T);
        add_to_N_spherocylinder (N, G, Rg, cyl_length, Volume,  SA, ProjSA,  dN,  Nmin,  len_N, NArray_Gmin, NArray_confs, NArray_quant);
        fscanf_result = fscanf(input_file, format, &Rg, &cyl_length, &N, &Volume, &SA, &ProjSA);
    }
    print_NGDataFile_spherocylinder (NArray_quant, output_file);
    fclose(output_file);
    fclose(input_file);
}


int main(int argc, const char * argv[])
{
    /* Energy related variables*/
    double Rho;         //Moles per m3
    double Mu;
    double Sigma;       //Liquid-crystal surface tension
    double T;           //Kelvin
    double theta_cg, theta_cb;
    
    int Nmin, Nmax;
    double dN;
    FILE* input_file;
    FILE* output_file;
    
    Nmin        = atoi(argv[1]);
    dN          = atof(argv[2]);
    Nmax        = atoi(argv[3]);
    Rho         = atof(argv[4]);   //Moles per m3
    Mu          = atof(argv[5]);
    Sigma       = atof(argv[6]);   //Liquid-crystal surface tension
    T           = atof(argv[7]);  //Kelvin
    theta_cg    = atof(argv[8]);
    theta_cb    = atof(argv[9]);
    input_file  = fopen(argv[10], "r");
    output_file = fopen(argv[11], "w");
    
    GetFreeEnergyProfile get_free_energy_profile (input_file, output_file, Nmin, Nmax, dN, Rho, Mu, Sigma, T, theta_cg, theta_cb);
    get_free_energy_profile.ReadProfileSphericalCaps();
    
    
    return 0;
}


//printf("%10.3f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%10.20f\t%d\t%d\n", N, Volume, SA, projected_SA[0], projected_SA[1], projected_SA[2], projected_SA[3], projected_SA[4], Radii[0], Radii[1], Radii[2], db, dg_secondary);
