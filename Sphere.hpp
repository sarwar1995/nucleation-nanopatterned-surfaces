//
//  Sphere.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Sphere_hpp
#define Sphere_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include "Shape.hpp"

#define pi 3.141592653589793

class Sphere:public Shape {
public:
    Sphere();
    ~Sphere();
    Sphere(double, std::vector<double>&);
    
    inline void calc_volume();
    inline void calc_SA();
    int isInside (std::vector<double>&); 
    int nearSurface(std::vector<double>&, double);
    double volume;
    double SA;
    double centre[3];
    double radius ;
    double centre_point_dist(std::vector<double>&);
protected:
    //double centre_point_dist(std::vector<double>&);
   
    
};
#endif /* Sphere_hpp */

