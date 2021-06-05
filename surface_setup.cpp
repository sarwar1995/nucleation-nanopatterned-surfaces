//
//  surface_setup.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "surface_setup.hpp"

SurfaceSetup::SurfaceSetup()
{
    
}

SurfaceSetup::~SurfaceSetup()
{
    
}

SurfaceSetup::SurfaceSetup(int num_patches, double z, double good_patch_x_width, double good_patch_y_width, double bad_patch_x_width, double bad_patch_y_width, double theta_good, double theta_bad): z_surface(z)
{
    /* Creating a patch factory will be beneficial */
    /* num_patches have to be odd, for a central good patch and surrounding bad/good patches type stripes pattern */
    if (num_patches % 2 == 0) {printf("Number of patches is not odd. You have %d patches\n", num_patches); throw num_patches;}
    else
    {
        double patch_centre_increment = (good_patch_x_width + bad_patch_x_width)/2.0 ;
        /* Setting the centreal good patch */
        std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
        centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
        std::vector<double> dim_Good(2,0.0);
        std::vector<double> dim_Bad(2,0.0);
        dim_Good[0] = good_patch_x_width; dim_Good[1] = good_patch_y_width; //Common dimensions for all Good patches
        dim_Bad[0] = bad_patch_x_width; dim_Good[1] = bad_patch_y_width; //Common dimensions for all Bad patches
        Patch good_patch (theta_good, centre_good, dim_Good);
        
        if(list_of_patches.empty() && orientations.empty())
        {
            list_of_patches.push_back(good_patch);
            orientations.push_back(centre_good);
            
        }
        else
        {
            throw "Patch list and orientations not empty during surface setup.";
        }
//        std::vector<Patch> list_of_patches {good_patch, bad_patch_left, bad_patch_right};
//        std::vector<std::vector<double> > orientations {centre_good, centre_bad_left, centre_bad_right};

        
        /* Setting the bad patches */
        std::vector<double> patch_centre_left(3); //(-a, 0, 0)
        std::vector<double> patch_centre_right(3);//(a, 0, 0)

        if (bad_patch_x_width == -1 && bad_patch_y_width == -1)
        {
//            int is_semi_infinite = 1;
//            int position_flag_left = -1;
//            int position_flag_right = 1;
//            Patch bad_patch_left (is_semi_infinite, position_flag_left);
//            Patch bad_patch_right (is_semi_infinite, position_flag_right);
//            centre_bad_left[0] = -(good_patch_x_width + bad_patch_x_width)/2.0; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
//            centre_bad_right[0] = (good_patch_x_width + bad_patch_x_width)/2.0; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
//
        }
        else
        {
            /* Implement a more general case of multiple patches perhaps using Patch factory if needed */
            int n_surrounding_patches = num_patches -  1 ;
            int n_patches_per_side_of_good = (n_surrounding_patches/2) ;
            for(int j=0; j<n_patches_per_side_of_good; j++) //Here j=0 corresponds to the patch next to the central good patch.
            {
                int patch_centre_modifier = j+1;
                patch_centre_left[0] = -patch_centre_modifier * patch_centre_increment; patch_centre_left[1] = 0.0; patch_centre_left[2] = z_surface;
                patch_centre_right[0] = patch_centre_modifier * patch_centre_increment; patch_centre_right[1] = 0.0; patch_centre_right[2] = z_surface;
                if (j%2== 0)
                {
                    Patch jth_bad_patch_left (theta_bad, patch_centre_left, dim_Bad);
                    Patch jth_bad_patch_right (theta_bad, patch_centre_right, dim_Bad);
                    list_of_patches.push_back(jth_bad_patch_left);
                    list_of_patches.push_back(jth_bad_patch_right);
                    orientations.push_back(patch_centre_left);
                    orientations.push_back(patch_centre_right);
                }
                else
                {
                    Patch jth_good_patch_left (theta_good, patch_centre_left, dim_Good);
                    Patch jth_good_patch_right (theta_good, patch_centre_right, dim_Good);
                    list_of_patches.push_back(jth_good_patch_left);
                    list_of_patches.push_back(jth_good_patch_right);
                    orientations.push_back(patch_centre_left);
                    orientations.push_back(patch_centre_right);
                }
            }
        }
        
    }    
}


Stripes SurfaceSetup::create_new_stripes()
{
    Stripes new_stripes (list_of_patches, orientations, z_surface);
    return(new_stripes);
}


//Code from main_parallel

/* Setting the central good patch */
//std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
//centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
//std::vector<double> dim_Good(2,0.0);
//std::vector<double> dim_Bad(2,0.0);
//
//dim_Good[0] = good_patch_x_width; dim_Good[1] = good_patch_y_width; //Common dimensions for all Good patches
//dim_Bad[0] = bad_patch_x_width; dim_Good[1] = bad_patch_y_width; //Common dimensions for all Bad patches
//Patch good_patch (theta_good, centre_good, dim_Good);
//
///* Setting the secondry good patches */
//std::vector<double> centre_good_left(3); //(-2a, 0, 0)
//std::vector<double> centre_good_right(3);//(2a, 0, 0)
//centre_good_left[0] = -(good_patch_x_width + bad_patch_x_width); centre_good_left[1] = 0.0; centre_good_left[2] = z_surface;
//centre_good_right[0] = good_patch_x_width + bad_patch_x_width; centre_good_right[1] = 0.0; centre_good_right[2] = z_surface;
//Patch good_patch_left (theta_good, centre_good_left, dim_Good);
//Patch good_patch_right (theta_good, centre_good_right, dim_Good);
//
///* Setting the bad patches */
//std::vector<double> centre_bad_left(3); //(-a, 0, 0)
//std::vector<double> centre_bad_right(3);//(a, 0, 0)
//centre_bad_left[0] = -(good_patch_x_width + bad_patch_x_width)/2.0; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
//centre_bad_right[0] = (good_patch_x_width + bad_patch_x_width)/2.0; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
//Patch bad_patch_left (theta_bad, centre_bad_left, dim_Bad);
//Patch bad_patch_right (theta_bad, centre_bad_right, dim_Bad);
//
///* Setting the surface */
//std::vector<Patch> list_of_patches {good_patch, bad_patch_left, bad_patch_right, good_patch_left, good_patch_right};
//std::vector<std::vector<double> > orientations {centre_good, centre_bad_left, centre_bad_right, centre_good_left, centre_good_right};
//
//Stripes stripes (list_of_patches, orientations, z_surface);
