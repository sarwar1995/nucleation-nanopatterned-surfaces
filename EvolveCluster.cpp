//
//  EvolveCluster.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "EvolveCluster.hpp"
EvolveCluster::EvolveCluster(){
    
}

EvolveCluster::~EvolveCluster()
{
    
}

EvolveCluster::EvolveCluster(CheckBoundary* boundary_check,  MC* mc, std::vector<double> volume_SA_mc, std::vector<double> maximum_limits, std::vector<double> increments, double theta_good, double theta_bad, std::vector<double> patch_widths, double z, double Delta, double rho):
check_boundary(boundary_check), mc_engine(mc), theta_cg(theta_good) ,theta_cb(theta_bad), z_surface(z), mc_volume_SA(volume_SA_mc), delta(Delta), Rho(rho)
{
    
 // maximum_limits = [Rgmax, Rbmax , Cylmax]
 // increments = [dRg, dRb , d_cyl_length]
 // patch_widths = [pG_width_x, pG_width_y, pB_width_x, pB_width_y]
    
    Rb_max = maximum_limits[1];
    d_Rb = increments[1];
    pG_width_x = patch_widths[0];
    pG_width_y = patch_widths[1];
    
    printf("Rb_max = %10.10f\t d_Rb= %10.10f\t  pG_width_x= %10.10f\t  pG_width_y= %10.10f\n", Rb_max, d_Rb, pG_width_x, pG_width_y);
    printf("Evolve cluster initialized\n");

}



void EvolveCluster::EvolveBadCapWithSpherocylinder (SpheroCylinder* spherocylinder, Shape* Cluster_shape_ptr, double cyl_length, double d_chord_length, double Rg, FILE* outputfile, FILE* output_points_file)
{
    printf("Inside evolve bad cap with spherocylinder\n");
    double chord_length_max = cyl_length;
    double chord_length;
    double Rb, projected_rb_min, projected_rb, Rb_min;
    int len_Rb;
    double dB ; //This is the distance between bad cap centre and good patch width
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0);
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0;
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface;
    
    for(int j=0; j <= chord_length_max; j++)
    {
        chord_length = 0.0 + j * d_chord_length;
        projected_rb_min = 0.5 * chord_length ;
        
        Rb_min = projected_rb_min/(double)sin(theta_cb) ; //Be careful to not have a theta that is 0 or 180
        len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;
        std::vector<int> dB_list = {-1, 1};
        // 1: centre on bad patches
        // -1: centres on good patches
        for(int k=0; k<2; k++)
        {
            for(int l=0; l<len_Rb; l++)
            {
                Rb = Rb_min + l * d_Rb ;
                printf("cyl length =%10.10f Rb = %10.10f\t dB_sign=%d\n", cyl_length, Rb, dB_list[k]);
                projected_rb = Rb * sin(theta_cb);
                
                int dB_sign = dB_list[k];
                double inside_sq = projected_rb*projected_rb - (0.5 * chord_length *0.5 * chord_length);
                if(abs(inside_sq)<1e-10 && inside_sq < 0.0) {inside_sq = 0.0;}
                
                dB = dB_sign * sqrt(inside_sq);
                c_bad_left[0] = -(pG_width_x/2.0) - dB;
                c_bad_right[0] = (pG_width_x/2.0) + dB;

                bool check_bad_cap = check_boundary->CheckSpherocylinderBadPatch(cyl_length, dB_sign, pG_width_x, c_bad_left, c_bad_right, projected_rb);

                if(!check_bad_cap)
                {
                    break;
                }
                else
                {
                    
                    
                    Spherical_cap BadCap_left (c_bad_left, projected_rb, theta_cb, z_surface);
                    Spherical_cap BadCap_right (c_bad_right, projected_rb, theta_cb, z_surface);
                    Spherical_cap* bad_cap_left_ptr = &BadCap_left;
                    Spherical_cap* bad_cap_right_ptr = &BadCap_right;
                    
                    std::vector<Spherical_cap*> caps ({bad_cap_left_ptr, bad_cap_right_ptr});
                    
                    
                    //std::vector<Spherical_cap*> caps = get_bad_caps(c_bad_left, c_bad_right, projected_rb);
                    Spherocylinder_cap_composite spherocylinder_composite (caps, spherocylinder);
                    
                    Cluster_shape_ptr = &spherocylinder_composite;
                    check_boundary->ManageBoxBreach(Cluster_shape_ptr);
                    
                    //Ideally here will also check for crossing of the patch boundary and what to do when that happens but for now
                    //implementing this for an infinite bad patch size.
                    mc_volume_SA = mc_engine->calc_volume_SA(Cluster_shape_ptr, delta);
                
                    double Volume = mc_volume_SA[0]; //GoodCap.getVolume();
                    double SA = mc_volume_SA[1];
                    
                    //This AnalytProjSurfAreaWithPlaneIntersect() is incorrect
                    double projected_SA = spherocylinder_composite.projected_SA();
                    
                    double N = (Volume * 1e-30) * Rho * avogadro ;
                    
                    if(dB_sign == -1 && l == len_Rb-2)
                    {
                        printf("inside printig surface points\n");
                        mc_engine->print_surf_points(output_points_file);
                        
                    }
                    fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%d\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, Rb, dB_sign, N, Volume, SA, projected_SA);
                }
            }
        }
        
    }
}


std::vector<Spherical_cap*> EvolveCluster::get_bad_caps(std::vector<double> c_bad_left, std::vector<double> c_bad_right, double projected_rb)
{
    std::vector<Spherical_cap*> result(2);
    Spherical_cap BadCap_left (c_bad_left, projected_rb, theta_cb, z_surface);
    Spherical_cap BadCap_right (c_bad_right, projected_rb, theta_cb, z_surface);
    
    result[0] = &BadCap_left;
    result[1] = &BadCap_right;
    return result;
    
}

