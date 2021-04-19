//
//  CheckBoundary.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef CheckBoundary_hpp
#define CheckBoundary_hpp
#include <stdio.h>
#include <stdlib.h>
#include "Surface.hpp"
#include "Shape.hpp"
#include "DynamicBox.hpp"
#include <math.h>


class CheckBoundary{

public:
    CheckBoundary();
    ~CheckBoundary();
    CheckBoundary(Surface*, MC* , DynamicBox*, double);
    void ManageBoxBreach(Shape*);
    
    std::vector<int> CheckPatchBoundaries();
    bool CheckSpherocylinderBadPatch(double cyl_length, int dB_sign, double patch_width, std::vector<double> centre_left, std::vector<double> centre_right, double projected_radius); //If this returns true the current cap size is good otherwise break the loop for bad patch cap growth

protected:
    Surface* surface;
    MC* mc_engine;
    Shape* cluster;
    DynamicBox* dynamic_box;
    double point_density;
    std::vector<int> FindBoxBreach(std::vector<int> box_breach);
    void ExpandBox (std::vector<int> find_box_breach);
};




#endif /* CheckBoundary_hpp */
