//
//  SphericalCapDivided.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef SphericalCapDivided_hpp
#define SphericalCapDivided_hpp

#include <stdio.h>
#include <vector>
#include "Shape.hpp"
#include "Spherical_cap.hpp"
using namespace std;

class Spherical_cap_divided:public Shape
{
    
public:
    Spherical_cap_divided();
    ~Spherical_cap_divided();
    
    Spherical_cap_divided(std::vector<std::vector<double>> centres, double r, std::vector<double> the_normal, double z_surface);

    int isInside (std::vector<double>&); //Tells whether a point is inside that shape
    int nearSurface(std::vector<double>&, double);
    
protected:
    std::vector<std::vector<double>> divided_cap_centres;

    Sphere sphere_left;
    Sphere sphere_right;
    
    int normal_direction;
    std::vector<double> normal;
    std::vector<double> bounds_along_normal;
    double radius;
    double z_wall;
    
    int InsideSphereLeft (std::vector<double>& point);
    int InsideSphereRight (std::vector<double>& point);
};
#endif /* SphericalCapDivided_hpp */
