//
//  FreeEnergy.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/7/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "FreeEnergy.hpp"

double free_energy (double Rho, double Mu, double Sigma, double volume, double SA, std::vector<double>& projected_SA, double theta_good, double theta_bad, double T)
{
    /* Here projected_SA has the projected surface area lying exclusively on the [0]: good, [1]:bad on left and [2]:bad on right patches, [3]: Good left and [4]: Good right patches */
    double G;
    double Vcomp = -1*volume* Rho * Mu;
    double Scomp = Sigma*SA;
    
    printf("volume = %10.10f\t Rho = %10.10f\t Mu=%10.10f\t Vcomp=%10.10f\t Scomp=%10.10f\n", volume*1e30, Rho, Mu, Vcomp/(kb*T), Scomp/(kb*T));
    
    double theta_patch ;
    double Sproj_comp = 0;
    
//    double SA_proj_good = projected_SA[0];
//    double SA_proj_bad_1 = projected_SA[1];
//    double SA_proj_bad_2 = projected_SA[2];
    
    int p_bad_start = 1;
    int p_good_start = 3;
    for(int p=0; p<projected_SA.size(); p++)
    {
        if (p==0)
        {
            theta_patch = theta_good;
        }
        else if (p==p_bad_start || p-p_bad_start == 1)
        {
            theta_patch = theta_bad;
            if(p-p_bad_start == 1){p_bad_start = p_bad_start + 4;}
        }
        else if(p==p_good_start || p-p_good_start == 1)
        {
            theta_patch = theta_good;
            if(p-p_good_start == 1){p_good_start = p_good_start + 4;}
        }
        else
        {
            printf("All the patches must be covered \n"); abort();
        }
        Sproj_comp = Sproj_comp + cos(theta_patch)*projected_SA[p];
        printf("p = %d theta_patch=%10.10f\t projected_SA[p]=%10.10f\t Sproj_comp=%10.10f\n", p, theta_patch, projected_SA[p], Sproj_comp);
    }
    
    Sproj_comp = -Sigma*Sproj_comp;
      //-Sigma*(cos(theta_good)*SA_proj_good + cos(theta_bad)*SA_proj_bad_1 + cos(theta_bad)*SA_proj_bad_2);
    G = Vcomp + Scomp + Sproj_comp;
    G = G/(kb*T) ;
    return G;
    
}

double free_energy_singlecap(double Rho, double Mu, double Sigma, double volume, double SA, double projected_SA, double theta, double T)
{
    double G;
    double Vcomp = -1*volume* Rho * Mu;
    double Scomp = Sigma*SA;
    double Sproj_comp = -Sigma*(cos(theta)*projected_SA);
    G = Vcomp + Scomp + Sproj_comp;
    G = G/(kb*T) ;
    return G;
}

double free_energy_spherocylinder (double Rho, double Mu, double Sigma, double volume, double SA, std::vector<double>& projected_SA, double theta_good, double theta_bad, double T)
{
    double G;
    double Vcomp = -1*volume* Rho * Mu;
    double Scomp = Sigma*SA;
    
    double theta_patch ;
    double Sproj_comp = 0;
    for(int p=0; p<projected_SA.size(); p++)
    {
        if (p==0 || p==2)
        {
            theta_patch = theta_bad;
        }
        else if (p==1)
        {
             theta_patch = theta_good;
        }
        else
        {
            printf("There are more than 3 projected SAs for spherocylinder \n"); abort();
        }
        Sproj_comp = Sproj_comp + cos(theta_patch)*projected_SA[p];
    }
    
    Sproj_comp = -Sigma*Sproj_comp;
    //-Sigma*(cos(theta_good)*SA_proj_good + cos(theta_bad)*SA_proj_bad_1 + cos(theta_bad)*SA_proj_bad_2);
    G = Vcomp + Scomp + Sproj_comp;
    G = G/(kb*T) ;
    return G;
}
