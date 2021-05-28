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

class SurfaceSetup {
    
public:
    SurfaceSetup();
    ~SurfaceSetup();
    
    SurfaceSetup(int num_patches, double z_surface, double good_patch_x_width, double good_patch_y_width, double bad_patch_x_width, double bad_patch_y_width, double theta_good, double theta_bad);
    

    Stripes create_new_stripes (); //This accessor can return the more general surface object in future.
    
protected:
    std::vector<Patch> list_of_patches;
    std::vector<std::vector<double> > orientations;
    //Stripes stripes; //This can be a more general surface object in future.
    
    
};


#endif /* surface_setup_hpp */
