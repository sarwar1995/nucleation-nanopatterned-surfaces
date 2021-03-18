//
//  TestHourglass.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef Cylinder_hpp
#define Cylinder_hpp

#include "Shape.hpp"
#include "Circle.hpp"
#include <string>
#include "Surface.hpp"
#include "miscformulas.hpp"
#include <stdio.h>

class Cylinder:public Shape
{
public:
    Cylinder();
    ~Cylinder();
    Cylinder(double, std::vector<std::vector<double> >&); //radius and centres of the two base circles
    
    int isWithin (std::vector<double>&); //Checks whether a point is in between the two circular planes of the cylinder
    int isInside (std::vector<double>&);
    int nearSurface(std::vector<double>&, double);
    
    //These three functions are for cylinder with axis parallel to the xy plane.
    double AnalytVolumeWithPlaneIntersect (double z); //plane parallel to xy axis i.e. z=constant
    double AnalytSurfAreaWithPlaneIntersect (double z);
    double AnalytProjSurfAreaWithPlaneIntersect (double z); 
    
protected:
    double rc;
    std::vector<std::vector<double> > centres; //Size of this must be 2. The two centres of the circles
    std::vector<double> normalized_axis;
    double length;
    Circle base_circle;
    Circle top_circle;
};


#endif /* Cylinder_hpp */
