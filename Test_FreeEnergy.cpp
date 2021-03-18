//
//  Test_FreeEnergy.cpp
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
    double Rg, Rb, Rg_secondary, theta_cg, theta_cb, pB_width_x, pB_width_y, pG_width_x, pG_width_y, box_height, delta;
    double Rho, Mu, Sigma, T;
    int aSeed[3];
    int db_list, dg_secondary_list;
    int n_points;
    Rg          = atof(argv[1]);
    Rb          = atof(argv[2]);
    Rg_secondary = atof(argv[3]);
    theta_cg     = atof(argv[4]);
    theta_cb     = atof(argv[5]);
    db_list      = atoi(argv[6]);
    dg_secondary_list = atoi(argv[7]);
    pG_width_x   = atof(argv[8]);
    pG_width_y   = atof(argv[9]);
    pB_width_x   = atof(argv[10]);
    pB_width_y   = atof(argv[11]);
    box_height  = atof(argv[12]);
    n_points    = atoi(argv[13]);
    n_points    = n_points*1e06;
    aSeed[0]    = atoi(argv[14]);       //Seed in x-direction for MC
    aSeed[1]    = atoi(argv[15]);       //Seed in y-direction for MC
    aSeed[2]    = atoi(argv[16]);       //Seed in z-direction for MC
    delta       = atof(argv[17]);
    Rho         = atof(argv[18]);   //Moles per m3
    Mu          = atof(argv[19]);
    Sigma       = atof(argv[20]);   //Liquid-crystal surface tension
    T           = atof(argv[21]);  //Kelvin
    double dN   = atof(argv[22]);
//    std::string tag (argv[20]);
//    std::string surf_fileName = tag + "_surfpoints.txt";
//    std::string surf_fileName_1 = tag + "-suface_points_goodcap.txt";
//    std::string surf_fileName_2 = tag + "-suface_points_badleftcap.txt";
//    std::string surf_fileName_3 = tag + "-suface_points_badrightcap.txt";
//
//    FILE* output_surface_file   = fopen(surf_fileName.c_str(), "w");
//    FILE* output_surface_file_1 = fopen(surf_fileName_1.c_str(), "w");
//    FILE* output_surface_file_2 = fopen(surf_fileName_2.c_str(), "w");
//    FILE* output_surface_file_3 = fopen(surf_fileName_3.c_str(), "w");
    double z_surface = 0.0;
    
    printf("n_points = %d\n", n_points);
    /* Setting the centreal good patch */
    std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
    centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
    std::vector<double> dim_Good(2,0.0);
    std::vector<double> dim_Bad(2,0.0);
    
    dim_Good[0] = pG_width_x; dim_Good[1] = pG_width_y; //Common dimensions for all Good patches
    dim_Bad[0] = pB_width_x; dim_Good[1] = pB_width_y; //Common dimensions for all Bad patches
    Patch good_patch (theta_cg, centre_good, dim_Good);
    
    /* Setting the secondry good patches */
    std::vector<double> centre_good_left(3); //(-2a, 0, 0)
    std::vector<double> centre_good_right(3);//(2a, 0, 0)
    centre_good_left[0] = -(pG_width_x + pB_width_x); centre_good_left[1] = 0.0; centre_good_left[2] = z_surface;
    centre_good_right[0] = pG_width_x + pB_width_x; centre_good_right[1] = 0.0; centre_good_right[2] = z_surface;
    Patch good_patch_left (theta_cg, centre_good_left, dim_Good);
    Patch good_patch_right (theta_cg, centre_good_right, dim_Good);
    
    /* Setting the bad patches */
    std::vector<double> centre_bad_left(3); //(-a, 0, 0)
    std::vector<double> centre_bad_right(3);//(a, 0, 0)
    centre_bad_left[0] = -(pG_width_x + pB_width_x)/2.0; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
    centre_bad_right[0] = (pG_width_x + pB_width_x)/2.0; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
    Patch bad_patch_left (theta_cb, centre_bad_left, dim_Bad);
    Patch bad_patch_right (theta_cb, centre_bad_right, dim_Bad);
    
    /* Setting the surface */
    std::vector<Patch> list_of_patches {good_patch, bad_patch_left, bad_patch_right, good_patch_left, good_patch_right};
    std::vector<std::vector<double> > orientations {centre_good, centre_bad_left, centre_bad_right, centre_good_left, centre_good_right};
    
    Stripes stripes (list_of_patches, orientations, z_surface);
    
    /* Calculating box. Here box[i] = 2-dimensional and represents the boundaries of the box in each direction. */
    stripes.calc_box(box_height);
    for(size_t i =0; i<stripes.box.size(); i++)
    {
        for(size_t j =0; j<stripes.box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
        }
    }
    
    /* Setting the MC engine */
    MC mc_engine (n_points, stripes.box, aSeed);
    
    
    printf("box_volume = %10.5f\n", mc_engine.box_volume());
    
    double projected_rg = Rg*sin(theta_cg); 
    Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
    
    double projected_rb_min = sqrt(projected_rg*projected_rg - ((pG_width_x*pG_width_x)/4.0));
    double Rb_min = projected_rb_min/(double)sin(theta_cb) ; //Be careful to not have a theta that is 0 or 180
    printf("Rb_min = %10.10f\n", Rb_min);
    
    double projected_rb = Rb * (double)sin(theta_cb);
    double inside_sq = (projected_rb*projected_rb) - (projected_rg*projected_rg) + ((pG_width_x*pG_width_x)/4.0);
    if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
    double d_B = db_list * sqrt(inside_sq) ;
    printf("d_B = %10.10f\n", d_B);
    
    double dG1 = pB_width_x - d_B ; /* b - dB1*/
    double projected_rg_secondary_min = sqrt(projected_rb*projected_rb - dG1*dG1);
    printf("projected_rg_secondary_min = %10.10f\n", projected_rg_secondary_min);
    
    double Rg_secondary_min = (projected_rg_secondary_min)/(double)sin(theta_cg);
    printf("Rg_secondary_min = %10.10f\n", Rg_secondary_min);
    
    
    /* Symmetric caps on either side of central patch */
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0);
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0;
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface;
    c_bad_left[0] = -(pG_width_x/2.0) - d_B;
    c_bad_right[0] = (pG_width_x/2.0) + d_B;
    
    std::vector<double> c_good_left(3,0.0), c_good_right(3, 0.0);
    double projected_rg_secondary = Rg_secondary * (double)sin(theta_cg);

    double inside_sq_Rg_second = (projected_rg_secondary*projected_rg_secondary) - (projected_rb*projected_rb) + (dG1*dG1);
    if(abs(inside_sq_Rg_second)<1e-10 && inside_sq_Rg_second < 0.0) {inside_sq = 0.0;}
    double d_G1 = dg_secondary_list * sqrt(inside_sq_Rg_second);
    c_good_left[0] = -(pB_width_x + (pG_width_x/2.0)) - d_G1;
    c_good_right[0] = (pB_width_x + (pG_width_x/2.0)) + d_G1;
    
    Spherical_cap BadCap_left(c_bad_left, projected_rb, theta_cb, z_surface);
    Spherical_cap BadCap_right(c_bad_right, projected_rb, theta_cb, z_surface);
    Spherical_cap GoodCap_left(c_good_left, projected_rg_secondary, theta_cg, z_surface);
    Spherical_cap GoodCap_right(c_good_right, projected_rg_secondary, theta_cg, z_surface);
    std::vector<Spherical_cap> capsList = {GoodCap, BadCap_left, BadCap_right, GoodCap_left, GoodCap_right};
    
    Composite_cluster Comp_cluster (capsList, stripes);
    std::vector<double> compcluster_projected_SA(capsList.size(),0.0);
    compcluster_projected_SA = Comp_cluster.projected_SAs();
    printf("compcluster_projected_SA size = %d\n", (int)compcluster_projected_SA.size());
    
    Shape* Cluster_shape_ptr = &Comp_cluster;

    std::vector<double> mc_volume_SA(2,0.0);
    mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
    //mc_engine.print_surf_points(output_surface_file);

    double Volume_MC = mc_volume_SA[0];
    Volume_MC = Volume_MC * 1e-30;
    double SA_MC = mc_volume_SA[1];
    SA_MC = SA_MC * 1e-20;
    for(size_t y=0; y<compcluster_projected_SA.size();y++){compcluster_projected_SA[y] = compcluster_projected_SA[y] * 1e-20;}
    double G = free_energy(Rho, Mu, Sigma, Volume_MC, SA_MC, compcluster_projected_SA, theta_cg, theta_cb, T);
    double N = Volume_MC*Rho*avogadro ;   //This volume needs to be in m3 for this to be valid since Rho is in moles/m3
    double rNhere = round_nearN (N, dN);
    printf("G = %10.15f N_actual=%10.2f N_round=%10.2f V=%10.15f SA=%10.15f\n",G ,N, rNhere, Volume_MC*1e30, SA_MC*1e20);
    printf("c_bad_left = %10.15f c_bad_right=%10.15f\n", c_bad_left[0], c_bad_right[0]);
    printf("compcluster_projected_SA[0], [1], [2] = [%10.15f %10.15f %10.15f]\n",compcluster_projected_SA[0]*1e20, compcluster_projected_SA[1]*1e20, compcluster_projected_SA[2]*1e20);
    double z_centre_bad = 0.0 - Rb*cos(theta_cb);
    double z_centre_good = 0.0 - Rg*cos(theta_cg);
    printf("z_centre_bad = %10.15f z_centre_good=%10.15f\n",z_centre_bad,z_centre_good);
    printf("projected_rg = %10.15f projected_rb = %10.15f\n",projected_rg, projected_rb);
   
    
    
    /* This part is for printing surface points for individual caps in the composite cluster*/
//    std::vector<double> mc_volume_SA_2(2,0.0);
//    MC mc_engine_2 (n_points, stripes.box, aSeed);
//    Shape* Goodcap_shape_ptr = &GoodCap;
//    Shape* BadLeftcap_shape_ptr = &BadCap_left;
//    Shape* BadRightcap_shape_ptr = &BadCap_right;
//    mc_volume_SA_2 = mc_engine_2.calc_volume_SA(Goodcap_shape_ptr, delta);
//    printf("Good cap V=%10.15f SA=%10.15f\n", mc_volume_SA_2[0], mc_volume_SA_2[1]);
//    mc_engine_2.print_surf_points(output_surface_file_1);
//
//    mc_volume_SA_2 = mc_engine_2.calc_volume_SA(BadLeftcap_shape_ptr, delta);
//    printf("Bad left cap V=%10.15f SA=%10.15f\n", mc_volume_SA_2[0], mc_volume_SA_2[1]);
//    mc_engine_2.print_surf_points(output_surface_file_2);
//
//    mc_volume_SA_2 = mc_engine_2.calc_volume_SA(BadRightcap_shape_ptr, delta);
//    printf("Bad rigt cap V=%10.15f SA=%10.15f\n", mc_volume_SA_2[0], mc_volume_SA_2[1]);
//    mc_engine_2.print_surf_points(output_surface_file_3);
//
//    fclose(output_surface_file);
//    fclose(output_surface_file_1);
//    fclose(output_surface_file_2);
//    fclose(output_surface_file_3);
    
}
