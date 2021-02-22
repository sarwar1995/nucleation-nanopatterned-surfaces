//
//  TestSphericalCap.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/12/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "MC.hpp"
#include <chrono>
#include "FreeEnergy.hpp"

//Define these but CONSIDER taking input
//double Rho = 5.4944e4 ; moles/m3
//double Mu = 641.4332101443232;
//double Sigma = 25.19e-03; // +- 2.95 Remains unchanged with Temperature
//double avogadro = 6.022140857e23;
//double kb = 1.38064852e-23;
//double T = 170.0;

using namespace std;

int main(int argc, const char * argv[])
{
    double Rho, Mu, Sigma, T;
    double Rgmax, d_Rg, theta_cg;
    double Rg;
    double z_surface = 0.0;
    double Volume, SA, projected_SA, G, N;
    double dN, Nmin;
    FILE* outputfile;
    Rgmax       = atof(argv[1]);
    d_Rg        = atof(argv[2]);
    theta_cg    = atof(argv[3]);
    Rho         = atof(argv[4]);
    Mu          = atof(argv[5]);
    Sigma       = atof(argv[6]);
    T           = atof(argv[7]);
    Nmin        = atof(argv[8]);
    dN          = atof(argv[9]);
    outputfile  = fopen(argv[10], "w");
    
    printf("Rho=%10.10f Mu=%10.10f Sigma=%10.10f T=%10.10f\n", Rho, Mu, Sigma, T);
    int len_Rg = int((Rgmax-0.0)/d_Rg) + 1;
    for(int i = 0; i<len_Rg; i++) // len_Rg
    {
        printf("i=%d\t",i);
        Rg = 0.0 + i*d_Rg ;
        std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
        centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        printf("projected_rg = %10.10f\n",projected_rg);
        
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        Volume = GoodCap.getVolume();
        Volume = Volume * 1e-30; //This converts the volume from Ang^3 to m^3
        SA = GoodCap.getSA();
        projected_SA = GoodCap.projected_SA();
        SA = SA* 1e-20; projected_SA = projected_SA * 1e-20;
        G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, projected_SA, theta_cg, T);
        N = Volume*Rho*avogadro ;   //This volume needs to be in m3 for this to be valid since Rho is in moles/m3
        double rNhere = round_nearN (N , dN);   //This rounds to nearest N.
        int indN = (int) ((rNhere-Nmin)/dN);
        int N_integer = (int) indN*dN ;
        fprintf(outputfile,"%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n",N_integer, Rg, projected_rg, G, Volume*1e30);
    }
    return 0;
}
