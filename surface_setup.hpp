//
//  surface_setup.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef surface_setup_hpp
#define surface_setup_hpp

#include <stdio.h>
#include <stdlib.h>
#include "Stripes.hpp"

#define INFINITE_PATCH 300
class SurfaceSetup {
    
public:
    SurfaceSetup();
    ~SurfaceSetup();
    
    SurfaceSetup(int num_patches, double z_surface, double good_patch_x_width, double good_patch_y_width, double bad_patch_x_width, double bad_patch_y_width, double theta_good, double theta_bad, int last_patch_isinfinite);
    

    Stripes create_new_stripes (); //This accessor can return the more general surface object in future.
    
protected:
    std::vector<Patch> list_of_patches;
    std::vector<std::vector<double> > orientations;
    double z_surface;
    //Stripes stripes; //This can be a more general surface object in future.
    
    int num_patches;
    double theta_bad, theta_good;
    std::vector<double> patch_centre_left; //(-a, 0, 0)
    std::vector<double> patch_centre_right;//(a, 0, 0)
    
    double good_patch_x_width, bad_patch_x_width;
    double good_patch_y_width, bad_patch_y_width;
    std::vector<double> dim_Good;
    std::vector<double> dim_Bad;
    std::vector<double> dim_Infinite;
    double patch_centre_increment;
    int n_surrounding_patches, n_patches_per_side_of_good;
    
    void compute_infinite_last_patch_centres(int, int);
    void compute_normal_patch_centre(int);
    void get_symmetric_patch();
    void add_patch(int, int);
    
    
};


#endif /* surface_setup_hpp */
