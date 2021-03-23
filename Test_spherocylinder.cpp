//
//  Test_spherocylinder.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/17/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "SpheroCylinder.hpp"
#include "FreeEnergy.hpp"

using namespace std;

int main(int argc, const char * argv[]) {

    //Test a single spherocylinder volume and surface area with analytical and MC values.
    //Density
    double Rho ;
    
    //Variables for spherocylinder
    double d_cyl_length;
    double cyl_max; //Maximum length of the cylinder
    
    //Variables for caps
    double d_Rg;
    double Rg_max;
    double Rb_max = 0.0;
    double theta_cg;
    double pG_width_x, pG_width_y;
    double Mu;
    double Sigma;       //Liquid-crystal surface tension
    double T;           //Kelvin
    
    
    pG_width_x  = atof(argv[1]);     //x width of the good patches (same for all good patches) (Angstroms)
    pG_width_y  = atof(argv[2]);     //y width of the good patches (same for all good patches) (Angstroms)
    theta_cg    = atof(argv[3]);       //Good patch contact angle
    d_Rg        = atof(argv[4]);       //increments in good patch cap radius (Angstroms i.e. d_Rg=0.05 A)
    d_cyl_length= atof(argv[5]);
    Rg_max      = atof(argv[6]);
    cyl_max     = atof(argv[7]);
    Rho         = atof(argv[8]);
    Mu          = atof(argv[9]);
    Sigma       = atof(argv[10]);
    T           = atof(argv[11]);
    
    FILE* outputfile = fopen(argv[12], "w");
    FILE* Gdatafile = fopen(argv[13], "w");
    
    
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1; //This 1 is added because for loop is i<len_Rg
    int len_cyl = (int) ((cyl_max-0.0)/d_cyl_length) + 1;
    double z_surface = 0.0;
    /* Setting the centreal good patch */
    std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
    centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
    std::vector<double> dim_Good(2,0.0);
    
    dim_Good[0] = pG_width_x; dim_Good[1] = pG_width_y; //Common dimensions for all Good patches
    Patch good_patch (theta_cg, centre_good, dim_Good);
    
    /* Setting the surface */
    std::vector<Patch> list_of_patches {good_patch};
    std::vector<std::vector<double> > orientations {centre_good};
    
    Stripes stripes (list_of_patches, orientations, z_surface);
    
    /* Setting the height of the box to R_max */
    double box_height = (Rg_max >= Rb_max) ? 2.0*Rg_max : 2.0*Rb_max ;
    //(theta_cg < pi/2.0) ? Rg_max : 2.0*Rg_max ;
    
    /* Calculating box. Here box[i] = 2-dimensional and represents the boundaries of the box in each direction. */
    stripes.calc_box(box_height);
    
    Shape* Cluster_shape_ptr;
    
    
    double Rg;
    double AnalytVolume, AnalytSA, AnalytProjSA;
    double Volume, SA, ProjSA;
    double N;
    double cyl_length;
    double G;
    
    printf("Rho =  %10.10f\n", Rho);
    printf("Mu =  %10.10f\n", Mu);
    printf("Sigma =  %10.10f\n", Sigma);
    printf("T =  %10.10f\n", T);
    
    for(int i = 0 ; i< len_Rg; i++) // 1
    {
        
        Rg = 0.0 + i*d_Rg ;
        printf("Rg =  %10.10f\n", Rg);
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        cyl_length = 0.0;
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        Cluster_shape_ptr= &GoodCap;
        std::vector<double> shape_xy_spread = Cluster_shape_ptr->xy_spread();
        std::vector<double> cPatchBounds = good_patch.patch_boundaries();
        //printf("shape_xy_spread = %10.10f\t%10.10f\t%10.10f\t%10.10f\n", shape_xy_spread[0], shape_xy_spread[1], shape_xy_spread[2], shape_xy_spread[3]);
        if(cPatchBounds[0] < shape_xy_spread[0] && cPatchBounds[1] > shape_xy_spread[1]) //Only spherical cap within the boundary
        {
            if((cPatchBounds[0] > shape_xy_spread[0] - d_Rg*sin(theta_cg)) && (cPatchBounds[0] <= shape_xy_spread[0]) && (cPatchBounds[1] < shape_xy_spread[1] + d_Rg*sin(theta_cg)) && (cPatchBounds[1] >= shape_xy_spread[1])  )
            {
                printf("Inside spherocylinder patch \n");
                for(int j = 0 ; j< len_cyl; j++) //j=0 is just the spherical cap that touches the boundary
                {
                    cyl_length = 0.0 + j*d_cyl_length;
                    std::vector<double> normal_to_the_circle({0.0 , 1.0, 0.0});
                    
                    SpheroCylinder sphero_cylinder (GoodCap, cyl_length, normal_to_the_circle, z_surface);
                    AnalytVolume = sphero_cylinder.AnalytVolumeWithPlaneIntersect();
                    AnalytSA = sphero_cylinder.AnalytSurfAreaWithPlaneIntersect();
                    AnalytProjSA = sphero_cylinder.AnalytProjSurfAreaWithPlaneIntersect();
                    Volume = AnalytVolume;
                    SA = AnalytSA;
                    ProjSA = AnalytProjSA;
                    
                    N = (Volume * 1e-30) * Rho * avogadro ;
                    G = free_energy_singlecap(Rho, Mu, Sigma, Volume*1e-30, SA*1e-20, ProjSA*1e-20, theta_cg, T);
                    fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, N, AnalytVolume, AnalytSA, AnalytProjSA);
                    fprintf(Gdatafile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\n", N, G, Rg, cyl_length);
                }
            }
            else
            {
                AnalytVolume = GoodCap.getVolume();
                AnalytSA = GoodCap.getSA();
                AnalytProjSA = GoodCap.projected_SA();
                
                Volume = AnalytVolume;
                SA = AnalytSA;
                ProjSA = AnalytProjSA;
                N = (Volume * 1e-30) * Rho * avogadro ;
                
                G = free_energy_singlecap(Rho, Mu, Sigma, Volume*1e-30, SA*1e-20, ProjSA*1e-20, theta_cg, T);
                
                fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, N, AnalytVolume, AnalytSA, AnalytProjSA);
                fprintf(Gdatafile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\n", N, G, Rg, cyl_length);
            }
        }
        
        else if (cPatchBounds[0] > shape_xy_spread[0] && cPatchBounds[1] < shape_xy_spread[1])
        {
            //Not allowing rg to grow beyond the good patch boundary
            break;
        }
        
    
    }
    
}

