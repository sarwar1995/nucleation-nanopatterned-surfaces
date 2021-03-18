//
//  BoxSizePrediction.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/22/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"


int main(int argc, char * argv[])
{
    double pG_width_x, pG_width_y, pB_width_x, pB_width_y, theta_cg, theta_cb;
    //int num_patches;
    double Rg_max, Rb_max;
    double z_surface = 0.0;
    
    pG_width_x  = atof(argv[1]);     //x width of the good patches (same for all good patches) (Angstroms)
    pG_width_y  = atof(argv[2]);     //y width of the good patches (same for all good patches) (Angstroms)
    pB_width_x  = atof(argv[3]);     //x width of the good patches (same for all good patches) (Angstroms)
    pB_width_y  = atof(argv[4]);     //y width of the good patches (same for all good patches) (Angstroms)
    theta_cg    = atof(argv[5]);       //Good patch contact angle
    theta_cb    = atof(argv[6]);       //Bad patch contact angle
    Rg_max      = atof(argv[7]);       //maximum limit of good patch cap.    (Angstroms)
    Rb_max      = atof(argv[8]);       //maximum limit of bad patch cap      (Angstroms)

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
    
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0);
    std::vector<double> c_good_left(3,0.0), c_good_right(3, 0.0);
    
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0; c_good_left[1] = 0.0; c_good_right[1] = 0.0;
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface; c_good_left[2] = z_surface; c_good_right[2] = z_surface;
    
    double projected_rg, projected_rb, projected_rg_secondary;
    
    projected_rg = Rg_max * (double)sin(theta_cg);
    
    double projected_rb_min = sqrt(projected_rg*projected_rg - ((pG_width_x*pG_width_x)/4.0));
    double Rb_min = projected_rb_min/(double)sin(theta_cb) ; //Be careful to not have a theta that is 0 or 180
    if(Rb_max > Rb_min)
    {
        projected_rb = Rb_max * (double)sin(theta_cb);
    }
    else
    {
        printf("Rb_min = %10.10f\n", Rb_min);
        projected_rb = Rb_min * (double)sin(theta_cb);
        printf("projected_rb = %10.10f\n", projected_rb);
    }
    

    double inside_sq = (projected_rb*projected_rb) - (projected_rg*projected_rg) + ((pG_width_x*pG_width_x)/4.0);
    printf("inside_sq = %10.10f\n", inside_sq);
    
    if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
    
    double d_B = 1.0 * sqrt(inside_sq) ;
    printf("d_B = %10.10f\n", d_B);
    /* Symmetric caps on either side of central patch */
    c_bad_left[0] = -(pG_width_x/2.0) - d_B;
    c_bad_right[0] = (pG_width_x/2.0) + d_B;
    
     printf("c_left = %10.10f c_right = %10.10f\n", c_bad_left[0], c_bad_right[0]);
    printf("projected_rg = %10.10f\n", projected_rg);
    printf("projected_rb = %10.10f\n", projected_rb);
   
    
    //    if(c_bad_left[0] >= 0.0 || c_bad_right[0] <= 0.0)
    //    {
    //        printf("cb1 = %10.10f\tcb2=%10.10f\n",c_bad_left[0],c_bad_right[0]);
    //        abort();
    //    }
    //
    //    if(!(projected_rg < projected_rb - c_bad_left[0] && projected_rg > projected_rb + c_bad_left[0] && projected_rg < projected_rb + c_bad_right[0] && projected_rg > projected_rb - c_bad_right[0]))
    //    {
    //        printf("Project_Rg=%10.10f Project_Rb=%10.10f dB=%10.10f dB_list=%d c_bad_left[0]=%10.10f c_bad_right[0]=%10.10f inside_sq=%10.20f\n", projected_rg, projected_rb, d_B, 1, c_bad_left[0], c_bad_right[0], inside_sq);
    //        abort();
    //    }
    
  
    
    /*In here we are looking at only the minimum of Rg secondary because Rg would only increase upto the box boundary */
    double dG1 = pB_width_x - d_B ; /* b - dB1*/
    double projected_rg_secondary_min = sqrt(projected_rb*projected_rb - dG1*dG1);
    double Rg_secondary_min = projected_rg_secondary_min/(double)sin(theta_cg) ; //Be careful to not have a theta that is 0 or 180
    
    projected_rg_secondary = Rg_secondary_min * (double)sin(theta_cg);
    
    printf("projected_rg_secondary = %10.10f\n", projected_rg_secondary);
    double inside_sq_Rg_second = (projected_rg_secondary*projected_rg_secondary) - (projected_rb*projected_rb) + (dG1*dG1);
    
    if(abs(inside_sq_Rg_second)<1e-10 && inside_sq_Rg_second < 0.0) {inside_sq_Rg_second = 0.0;}
    double d_G1 = 1 * sqrt(inside_sq_Rg_second);
    
    c_good_left[0] = -(pB_width_x + (pG_width_x/2.0)) - d_G1;
    c_good_right[0] = (pB_width_x + (pG_width_x/2.0)) + d_G1;
    
    printf("c_good_left = %10.10f c_good_right = %10.10f\n", c_good_left[0], c_good_right[0]);
    
    //printf("projected_rg_second=%10.10f\t c_good_left=%10.10f\t c_good_right=%10.10f\t", projected_rg_secondary, c_good_left[0], c_good_right[0]);
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> stripes_box_breach; //array of (0,1) to check box surface breach
    bool good_left_and_bad_left = ((projected_rb - (c_bad_left[0] - c_good_left[0]) < projected_rg_secondary) && projected_rg_secondary < projected_rb + (c_bad_left[0] - c_good_left[0]));
    
    bool good_right_and_bad_right = (projected_rb - (c_good_right[0] - c_bad_right[0]) < projected_rg_secondary && projected_rg_secondary < projected_rb + (c_good_right[0] - c_bad_right[0]) ) ;
    
    printf("good_left_and_bad_left = %d good_right_and_bad_right=%d Before clstrs\n", good_left_and_bad_left, good_right_and_bad_right);
    
    
//    if(!(good_left_and_bad_left && good_right_and_bad_right))
//    {
//        printf("good_left_bad_left =%d\t good_right_and_bad_right=%d\n",good_left_and_bad_left, good_right_and_bad_right);
//        if(good_left_and_bad_left != good_right_and_bad_right)
//        {
//            printf("bad_left = %10.5f\t good_left = %10.5f\t bad_right = %10.5f\t good_right = %10.5f\tn", c_bad_left[0], c_good_left[0], c_bad_right[0], c_good_right[0]);
//            printf("Symmetry is broken\n"); abort();
//        }
//        abort();
//    }
    
    Shape* Cluster_shape_ptr;
    
    Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
    Spherical_cap BadCap_left(c_bad_left, projected_rb, theta_cb, z_surface);
    Spherical_cap BadCap_right(c_bad_right, projected_rb, theta_cb, z_surface);
    Spherical_cap GoodCap_left(c_good_left, projected_rg_secondary, theta_cg, z_surface);
    Spherical_cap GoodCap_right(c_good_right, projected_rg_secondary, theta_cg, z_surface);
    
     printf("after caps\n");
    
    std::vector<Spherical_cap> capsList_good_second = {GoodCap, BadCap_left, BadCap_right, GoodCap_left, GoodCap_right};
    
    Composite_cluster Comp_cluster_good_second (capsList_good_second, stripes);
         printf("after cluster\n");
    Cluster_shape_ptr = &Comp_cluster_good_second;
    
//    stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
//    stripes_box_breach = stripes.monitor_box_breach(Cluster_shape_ptr);
    
    printf("Before spreads calculation\n");
    std::vector<double> cluster_xy_spread = Cluster_shape_ptr -> xy_spread();
    std::vector<double> cluster_3d_spread = Cluster_shape_ptr -> threeDim_spread();
    
    for(size_t i = 0; i<cluster_xy_spread.size(); i++)
    {
        printf("%d, xy-bound=%10.5f\n", (int)i, cluster_xy_spread[i]);
    }
    
    for(size_t i = 0; i<cluster_3d_spread.size(); i++)
    {
        printf("%d, 3d-bound=%10.5f\n", (int)i, cluster_3d_spread[i]);
    }
    
}
