//
//  Stripes.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Stripes.hpp"


Stripes::Stripes()
{
    
}

Stripes::Stripes(std::vector<Patch>& list_of_patches, std::vector<std::vector<double> >& patches_orientaions, double z):Surface(list_of_patches, patches_orientaions)
{
    central_patch=list_of_patches[0]; //This is assuming that list of patches[0] is the central good patch. Do an error check for this.
    std::vector<double> origin {0.0, 0.0, 0.0};
    if(!central_patch.check_patch_centre(origin))
    {
        throw std::invalid_argument( "zeroth patch in list is not at origin" );
    }
    z_wall = z;
    num_patches = (int)list_of_patches.size();
    n_surrounding_patches = num_patches -  1 ;
    n_patches_per_side_of_good = (n_surrounding_patches/2) ;
    n_unique_patches = n_patches_per_side_of_good + 1;
    
    bounds_crossed.resize(num_patches, 0);
}

Stripes::~Stripes()
{
    Patches.clear();
    orientations.clear();
}

//int Stripes::valid_point(std::vector<double> & point)
//{
//    if(point[3] >= z_wall)
//    {
//        return 1;
//    }
//    else
//    {
//        return 0;
//    }
//}

void Stripes::initial_box(double initial_box_dim)
{
    std::vector<double> central_patch_centre = central_patch.get_patch_centre(); //This will be origin in most cases.
    double centre_x = central_patch_centre[0];
    double centre_y = central_patch_centre[1];
    if(box.empty())
    {
        box.resize(3);
        for(size_t i=0; i<3; i++)
        {
            box[i].resize(2);
        }
        box[0][0] = centre_x - (initial_box_dim/2.0);
        box[0][1] = centre_x + (initial_box_dim/2.0);
        box[1][0] = centre_y - (initial_box_dim/2.0);
        box[1][1] = centre_y + (initial_box_dim/2.0);
        box[2][0] = z_wall;
        box[2][1] = initial_box_dim;
    }
    else
    {
        throw "Box is not empty before initial box has been calculated in Stripes\n" ;
    }
}

void Stripes::calc_box(double box_z){
    box.resize(3);
    for(size_t i=0; i<3; i++)
    {
        box[i].resize(2);
    }
//    double patch_x = Patches[0].get_x();    //Assuming all patches have the same dims, not equal x and y.
//    double patch_y = Patches[0].get_y();
    std::vector<double> cPatchBounds = central_patch.patch_boundaries();
    box[0][0] = cPatchBounds[0];
    box[0][1] = cPatchBounds[1];
    box[1][0] = cPatchBounds[2];
    box[1][1] = cPatchBounds[3];
    for(size_t i=0; i<Patches.size(); i++)
    {
        if(Patches[i].isEqual(central_patch)){continue;}
        else
        {
            std::vector<double> iPatchBounds;
            iPatchBounds = Patches[i].patch_boundaries();
            if(iPatchBounds[0] <= box[0][0]){box[0][0] = iPatchBounds[0];}
            if(iPatchBounds[1] >= box[0][1]){box[0][1] = iPatchBounds[1];}
            if(iPatchBounds[2] <= box[1][0]){box[1][0] = iPatchBounds[2];}
            if(iPatchBounds[3] >= box[1][1]){box[1][1] = iPatchBounds[3];}
        }
    }
    
    //Increasing the box size slightly to allow for obtuse spherical caps to fit within the box. This is only done for the case when the boundary patches are bad i.e. (patches.size()-1) is not divisible by 4
    
    //double buffer_x  = 0.1*(box[0][1] - box[0][0]) ;          //This buffer is a temporary fix to avoid box breach even before the patch boundary is crossed, as can happen in the case of obtuse angled caps on bad patches.
//    box[0][0] -= buffer;
//    box[0][1] += buffer;
//    box[1][0] -= buffer;
//    box[1][1] += buffer;
//
    
    box[2][0] = z_wall;
    box[2][1] = box_z;
}

std::vector<int> Stripes::monitor_cluster_spread(Shape* cluster)
{
    /*
     1:if that patch boundary is crossed (symmetrically. If it is crossed in one direction then it must also be crossed in the other)
    0:if not
     */
    /* 0: central patch, 1: left bad, 2:right bad, 3: left left good, 4: right right good and so on...*/
    std::vector<double> shape_xy_spread = cluster->xy_spread();
//    printf("Cluster xy spread = [%10.5f %10.5f][%10.5f %10.5f]\n", shape_xy_spread[0], shape_xy_spread[1], shape_xy_spread[2], shape_xy_spread[3]);
    
    //check central patch bound
    std::vector<double> cPatchBounds = central_patch.patch_boundaries();
    
    /* Only checking for x bounds as y-bounds for Stripes are expected to always be satisfied i.e. infinite in y-direction at least for stripes*/
    if(cPatchBounds[0] <= shape_xy_spread[0] && cPatchBounds[1] >= shape_xy_spread[1])
    {
        bounds_crossed[0] = 0;
    }
    else
    {
        if(cPatchBounds[0] > shape_xy_spread[0] && cPatchBounds[1] < shape_xy_spread[1])
        {
            bounds_crossed[0] = 1;
        }
        else
        {
            printf("Both boundaries of the central patch are not crossed. Symmetry error.\n");
            abort();
        }
    }
   
    if(Patches.size() >= 3) //i.e. if it is 3 or more
    {
        /* Always odd number of patches for the stripes case. Leftmost = size-2. Rightmost = size - 1*/
        for(size_t i=1; i<Patches.size(); i++)
        {
             std::vector<double> iPatchBounds = Patches[i].patch_boundaries();
            //For a given bounds_crossed. For patches beside the central patch we only care about crossing of the boundary in one direction because the centres will also shift and the other boundary is already crossed ex: in case of the left bad patch.
            // if i is odd (left patches) - compare patchbounds[0] with spread[0]
            // if i is even (right patches) - compare patchbounds[1] with spread[1]
            if(i%2 == 0)
            {
                if(iPatchBounds[1] >= shape_xy_spread[1])
                {
                    bounds_crossed[(int)i] = 0;
                }
                else{bounds_crossed[(int)i] = 1;}
            }
            else
            {
                if(iPatchBounds[0] <= shape_xy_spread[0])
                {
                    bounds_crossed[(int)i] = 0;
                }
                else{bounds_crossed[(int)i] = 1;}
            }

        }
        
    }
    return bounds_crossed;
}


bool Stripes::is_bounds_crossed(int unique_patch_id)
{
    int symmetric_left_patch = 2*unique_patch_id - 1;
    int symmetric_right_patch = symmetric_left_patch + 1;
    if(bounds_crossed[symmetric_left_patch] != bounds_crossed[symmetric_right_patch])
    {
        throw "bounds being crossed assymmetrically\n";
    }
    else if(bounds_crossed[symmetric_left_patch] == 1 && bounds_crossed[symmetric_right_patch] == 1 )
    { return true; }
    else {return false;}
}

bool Stripes::surface_bounds_breach (Shape* cluster)
{
    bounds_crossed = monitor_cluster_spread(cluster);
    if(Patches.size()>2)
    {
        if (bounds_crossed[(int) Patches.size() - 1] !=  bounds_crossed[(int) Patches.size() - 2])
        {printf ("bouds crossed calculation within surface bounds breached is not symmetric\n"); abort();}
        else if(bounds_crossed[(int) Patches.size() - 1] == 1 &&  bounds_crossed[(int) Patches.size() - 2] == 1)
        {return true;}
        else
        {return false;}
    }
    return true;
}

std::vector<int> Stripes::monitor_box_breach(Shape* cluster)
{
    /* Checking whether the cluster breaches any of the 5 planes of the box i.e. 2 y-z, 2 x-z and one top x-y plane. The bottom surface is always supposed to be satisfied due to the way the cluster is set. */
    
    std::vector<int> box_breach (5,0.0);
    std::vector<double> cluster_threeDim_spread = cluster->threeDim_spread();
    if (cluster_threeDim_spread[0] < box[0][0]){ box_breach[0] = 1 ;}
    if (cluster_threeDim_spread[1] > box[0][1]){ box_breach[1] = 1 ;}
    if (cluster_threeDim_spread[2] < box[1][0]){ box_breach[2] = 1 ;}
    if (cluster_threeDim_spread[3] > box[1][1]){ box_breach[3] = 1 ;}
    if (cluster_threeDim_spread[4] > box[2][1]){ box_breach[4] = 1 ;}
    return box_breach;
}


void Stripes::print_all_patch_bounds()
{
    std::vector<double> patch_bounds;
    for(size_t i=0; i<Patches.size() ; i++)
    {
        patch_bounds = Patches[i].patch_boundaries();
        printf("%zu th patch bounds = [%10.3f %10.3f][%10.3f %10.3f]\n", i, patch_bounds[0], patch_bounds[1], patch_bounds[2], patch_bounds[3]);
    }
}
