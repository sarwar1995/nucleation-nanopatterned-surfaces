//
//  Patch.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Patch.hpp"


Patch::Patch(){
}

Patch::~Patch(){
}

Patch::Patch(double theta_c, std::vector<double>& centre, std::vector<double>& dims){
    dimensions = dims;
    contact_angle = theta_c;
    patch_centre.resize(centre.size());
    patch_centre[0] = centre[0]; //X
    patch_centre[1] = centre[1]; //Y
}

std::vector<double> Patch::patch_boundaries() //0:Xmin, 1:Xmax, 2:Ymin, 3:Ymax
{
    std::vector<double> result(4,0.0);
    result[0] = patch_centre[0] - (0.5*dimensions[0]); //left bound
    result[1] = patch_centre[0] + (0.5*dimensions[0]); // right bound
    result[2] = patch_centre[1] - (0.5*dimensions[1]);
    result[3] = patch_centre[1] + (0.5*dimensions[1]);
    return result;
}

bool Patch::isEqual(Patch & B)
{
    return patch_centre == B.patch_centre ;
}

bool Patch::check_patch_centre(std::vector<double> point)
{
    return patch_centre == point;
}

