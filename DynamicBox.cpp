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


int DynamicBox::add_fix_box(int direction_to_expand_in, double density_of_points)
{
    //This function adds a fixed size box symmetrically (except in z) to the existing box in the direction that expansion needs to happen in.
    // direction_to_expand_in = 0 (expand along x-axis)
    // direction_to_expand_in = 1 (expand along y-axis)
    // direction_to_expand_in = 2 (expand along z-axis)
    
    double volume_added;
    int particles_added;
    
    double existing_lower_bound = box[direction_to_expand_in][0];
    double existing_upper_bound = box[direction_to_expand_in][1];
    
    double new_lower_bound, new_upper_bound;
    if(direction_to_expand_in == 2)
    {
        new_lower_bound = existing_lower_bound; //Because lower bound of z is the wall which stays the same
        new_upper_bound = existing_upper_bound + fixed_expansion_length ;
        volume_added = dimensions[0] * dimensions[1] * fixed_expansion_length ;
    }
    else
    {
        new_lower_bound = existing_lower_bound - fixed_expansion_length ;
        new_upper_bound = existing_upper_bound + fixed_expansion_length ;
        
        volume_added = 1.0;
        for(int i=0; i<3; i++)
        {
            if(i != direction_to_expand_in)
            {
                volume_added *= dimensions[i] ;
                
            }
            else
            {
                volume_added *= fixed_expansion_length;
            }
        }
        
    }
    box[direction_to_expand_in][0] = new_lower_bound;
    box[direction_to_expand_in][1] = new_upper_bound;
    
    particles_added = (int) (volume_added * density_of_points) ;
    return particles_added;
    
}
