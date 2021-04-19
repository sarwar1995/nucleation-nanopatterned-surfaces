//
//  DynamicBox.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "DynamicBox.hpp"

DynamicBox::DynamicBox()
{
    
}

DynamicBox::~DynamicBox()
{
    box.clear();
    dimensions.clear();
}

DynamicBox::DynamicBox(std::vector<std::vector<double>> existing_box, double length)
:box(existing_box)
,fixed_expansion_length(length)
{
    dimensions.resize(3);
    dimensions[0] = box[0][1] - box[0][0];
    dimensions[1] = box[1][1] - box[1][0];
    dimensions[2] = box[2][1] - box[2][0];
}


void DynamicBox::updateDimensions()
{
    dimensions[0] = box[0][1] - box[0][0];
    dimensions[1] = box[1][1] - box[1][0];
    dimensions[2] = box[2][1] - box[2][0];
}

int DynamicBox::add_fix_box(int direction_to_expand_in, double density_of_points)
{
    //This function adds a fixed size box symmetrically (except in z) to the existing box in the direction that expansion needs to happen in.
    // direction_to_expand_in = 0 (expand along x-axis)
    // direction_to_expand_in = 1 (expand along y-axis)
    // direction_to_expand_in = 2 (expand along z-axis)
    
    double volume_added_per_additional_box;
    int particles_added_per_additional_box;
    
    double existing_lower_bound = box[direction_to_expand_in][0];
    double existing_upper_bound = box[direction_to_expand_in][1];
    
    double new_lower_bound, new_upper_bound;
    if(direction_to_expand_in == 2)
    {
        new_lower_bound = existing_lower_bound; //Because lower bound of z is the wall which stays the same
        new_upper_bound = existing_upper_bound + fixed_expansion_length ;
        volume_added_per_additional_box = dimensions[0] * dimensions[1] * fixed_expansion_length ;
    }
    else
    {
        new_lower_bound = existing_lower_bound - fixed_expansion_length ;
        new_upper_bound = existing_upper_bound + fixed_expansion_length ;
        
        volume_added_per_additional_box = 1.0;
        for(int i=0; i<3; i++)
        {
            if(i != direction_to_expand_in)
            {
                volume_added_per_additional_box *= dimensions[i] ;
                
            }
            else
            {
                volume_added_per_additional_box *= fixed_expansion_length;
            }
        }
        
    }
    box[direction_to_expand_in][0] = new_lower_bound;
    box[direction_to_expand_in][1] = new_upper_bound;
    
    particles_added_per_additional_box = (int) (volume_added_per_additional_box * density_of_points) ;
    updateDimensions();
    return particles_added_per_additional_box;
    
}


std::vector<int> DynamicBox::CheckBoxBreach(Shape* cluster)
{
    /* Checking whether the cluster breaches any of the 5 planes of the box i.e. 2 y-z, 2 x-z and one top x-y plane. The bottom surface is always supposed to be satisfied due to the way the cluster is set. */
    std::vector<int> box_breach (5,0.0);
    std::vector<double> cluster_threeDim_spread = cluster->threeDim_spread();
    if (cluster_threeDim_spread[0] < box[0][0]){ box_breach[0] = 1 ;} //yz-left
    if (cluster_threeDim_spread[1] > box[0][1]){ box_breach[1] = 1 ;} //yz-right
    if (cluster_threeDim_spread[2] < box[1][0]){ box_breach[2] = 1 ;} //xz-bottom
    if (cluster_threeDim_spread[3] > box[1][1]){ box_breach[3] = 1 ;} //xz-top
    if (cluster_threeDim_spread[4] > box[2][1]){ box_breach[4] = 1 ;} //xy-top
    return box_breach;
}


void DynamicBox::print_box()
{
    for(size_t i =0; i<box.size(); i++)
    {
        for(size_t j =0; j<box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j, box[i][j]);
        }
    }
}
