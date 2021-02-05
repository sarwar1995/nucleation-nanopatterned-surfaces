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
    
    int isWithin (std::vector<double>&);
    int isInside (std::vector<double>&);
    int nearSurface(std::vector<double>&, double);
protected:
    double rc;
    std::vector<std::vector<double> > centres; //Size of this must be 2
    
};


#endif /* Cylinder_hpp */
