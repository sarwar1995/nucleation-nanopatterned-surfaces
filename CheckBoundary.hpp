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
#include <vector>
#include "Surface.hpp"
#include "Shape.hpp"
#include "DynamicBox.hpp"
//#include "MC.hpp"
#include "MC_parallel.hpp"
#include <math.h>


class CheckBoundary{
public:
    CheckBoundary();
    ~CheckBoundary();
    CheckBoundary(Surface*, MC* , DynamicBox*, double);
    void ManageBoxBreach(Shape*);
    int ManageBoxBreach_and_calc_Vol_SA(Shape*);
    
    std::vector<int> CheckPatchBoundaries();
    bool CheckSpherocylinderBadPatch(double cyl_length, int dB_sign, double patch_width, std::vector<double> centre_left, std::vector<double> centre_right, double projected_radius); //If this returns true the current cap size is good otherwise break the loop for bad patch cap growth
    
    bool CheckSphericalCapsOnInfiniteBadPatch(); //Only condition that needs to be checked for an infinite bad patch is that bad patch centres and less than and greater than zero respectively.
    
    double get_point_density() {return point_density;}
    
    void reset_n_breach_cycles();
protected:
    Surface* surface;
    MC* mc_engine;
    Shape* cluster;
    DynamicBox* dynamic_box;
    double point_density;
    std::vector<int> FindBoxBreach (std::vector<int>);
    void ExpandBox (std::vector<int> find_box_breach);
    void ExpandBox_and_calc_Vol_SA (std::vector<int> find_box_breach);
    
    //Calculating whether the box was breached at least once in the recursion;
    int n_breach_cycles;
};




#endif /* CheckBoundary_hpp */
