//
//  TestOnlyCap.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/15/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC.hpp"
#include <chrono>
#include "FreeEnergy.hpp"

int main(int argc, const char * argv[])
{
    
    double R, theta, delta, points_density;
    int n_points_max;
    int aSeed[3];
    R               = atof(argv[1]);
    theta           = atof(argv[2]);
    n_points_max    = atoi(argv[3]);
//    points_density  = atof(argv[4]);
    delta           = atof(argv[4]);
    aSeed[0]    = atoi(argv[5]);       //Seed in x-direction for MC
    aSeed[1]    = atoi(argv[6]);       //Seed in y-direction for MC
    aSeed[2]    = atoi(argv[7]);       //Seed in z-direction for MC
    FILE* output = fopen(argv[8], "w");
    
    std::vector<std::vector<double> > box;
    box.resize(3);
    for(size_t i=0; i<3; i++)
    {
        box[i].resize(2);
    }

    double z_surface = 0.0;
    std::vector<double> centre(3);
    centre[0] = 0.0; centre[1] = 0.0; centre[2] = z_surface;
    double projected_rg = R*sin(theta);
    printf("projected_rg = %10.10f\n",projected_rg);
    
    double box_dim = (theta <= (pi/2)) ? 1.5*projected_rg : 2*R ;
    double box_z = 2*R;
    box[0][0] = -box_dim;   box[1][0] = -box_dim;  box[2][0] = 0.0 ;
    box[0][1] = box_dim;   box[1][1] = box_dim;  box[2][1] = box_z;
    for(size_t i =0; i<box.size(); i++)
    {
        for(size_t j =0; j<box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,box[i][j]);
        }
    }
    
    double box_volume = 2*box_dim * 2*box_dim * box_z;
    printf("box_volume = %10.10f\n", box_volume);
    Spherical_cap Cap (centre, projected_rg, theta, z_surface);
    double Volume = Cap.getVolume();
    double SA = Cap.getSA();
    Shape* Cluster_shape_ptr = &Cap ;
    std::vector<double> mc_volume_SA(2,0.0);
    double SA_projected = Cap.projected_SA();
    printf("Analytical: \n");
    printf("V =%10.15f SA=%10.15f SA_projected=%10.15f\n", Volume, SA, SA_projected);
    
    double delta_V, delta_SA;
    for(int n=1; n<=n_points_max; n++)
    {
        int n_points = n * 1e06 ;
        printf("n_points = %d\n",n_points);
        double density = (double) (n_points/box_volume) ;
        printf("points denstiy = %10.10f\n", density);
        MC mc_engine (n_points, box, aSeed);
        mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
        double Volume_MC = mc_volume_SA[0];
        double SA_MC = mc_volume_SA[1];
        delta_V = Volume_MC - Volume;
        delta_SA = SA_MC - SA;
        printf("MC: \n");
        fprintf(output, "%10.15f\t%10.15f\t%10.15f\t%10.15f\t%10.15f\n", density, Volume_MC, SA_MC, delta_V, delta_SA);
    }
  
    fclose(output);
}
