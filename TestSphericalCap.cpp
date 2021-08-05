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
    int n_points;
    int aSeed[3];
    double delta;
    
    Rgmax       = atof(argv[1]);
    d_Rg        = atof(argv[2]);
    theta_cg    = atof(argv[3]);
    Rho         = atof(argv[4]);
    Mu          = atof(argv[5]);
    Sigma       = atof(argv[6]);
    T           = atof(argv[7]);
    Nmin        = atof(argv[8]);
    dN          = atof(argv[9]);
    n_points    = atoi(argv[10]);
    n_points    = n_points*1e06;
    aSeed[0]    = atoi(argv[11]);       //Seed in x-direction for MC
    aSeed[1]    = atoi(argv[12]);       //Seed in y-direction for MC
    aSeed[2]    = atoi(argv[13]);       //Seed in z-direction for MC
    delta       = atof(argv[14]);
    outputfile  = fopen(argv[15], "w");
    
    std::vector<std::vector<double> > box(3);
    for(size_t i=0; i<3; i++)
    {
        box[i].resize(2);
    }
    box[0][0] = -1.5*Rgmax;   box[0][1] = 1.5*Rgmax;
    box[1][0] = -1.5* Rgmax;   box[1][1] = 1.5* Rgmax;
    box[2][0] = z_surface;            box[2][1] = (2.0* Rgmax);
    
    // Normal constructor for random points
    //MC mc_engine (n_points, box, aSeed);
    //Specific constructor for lattic points
    MC mc_engine (n_points, box, aSeed, 1);
    Shape* Shape_ptr;
    std::vector<double> mc_volume_SA(2,0.0);
    double box_volume = mc_engine.box_volume();
    double density = ((double)n_points/box_volume);
    printf("n_points = %d\t box_volume = %10.5f\t density = %10.5g\n", n_points, box_volume, density);
    printf("Rho=%10.10f Mu=%10.10f Sigma=%10.10f T=%10.10f\n", Rho, Mu, Sigma, T);
    int len_Rg = int((Rgmax-0.0)/d_Rg) + 1;
    for(int i = 0; i<len_Rg; i++) // len_Rg
    {
        Rg = 0.0 + i*d_Rg ;
        std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
        centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        //printf("i=%d\t projected_rg = %10.10f\n",i, projected_rg);
        
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        Shape_ptr = &GoodCap;
//        mc_volume_SA = mc_engine.calc_volume_SA(Shape_ptr, delta);
        
        Volume =  GoodCap.getVolume(); //mc_volume_SA[0];
        Volume = Volume * 1e-30; //This converts the volume from Ang^3 to m^3
        SA =  GoodCap.getSA(); //mc_volume_SA[1] ;
        projected_SA = GoodCap.projected_SA();
        SA = SA* 1e-20; projected_SA = projected_SA * 1e-20;
        G = free_energy_singlecap(Rho, Mu, Sigma, Volume, SA, projected_SA, theta_cg, T);
        N = Volume*Rho*avogadro ;   //This volume needs to be in m3 for this to be valid since Rho is in moles/m3
        double rNhere = round_nearN (N , dN);   //This rounds to nearest N.
        int indN = (int) ((rNhere-Nmin)/dN);
        int N_integer = (int) indN*dN ;
        fprintf(outputfile,"%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n",N_integer, Rg, projected_rg, G, Volume*1e30, SA*1e20, projected_SA*1e20);
    }
    return 0;
}
