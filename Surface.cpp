//
//  Surface.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Surface.hpp"
Surface::Surface()
{
    
}
Surface::Surface(std::vector<Patch>& list_of_patches, std::vector<std::vector<double> >& patches_orientaions)
{
    Patches = list_of_patches;
    orientations = patches_orientaions;
    box.resize(0);
}

Surface::~Surface(){
    Patches.clear();
    orientations.clear();
    box.clear();
}

