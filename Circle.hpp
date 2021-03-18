//
//  Circle.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/6/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Circle_hpp
#define Circle_hpp

#include <stdio.h>
#include "Sphere.hpp"
#include "Plane.hpp"
#include "miscformulas.hpp"

class Circle{

public:
    Circle();
    ~Circle();
    Circle(std::vector<double>& , double);  //Point radius form with normal in the z-direction i.e. the circle is in the xy plane
    Circle(std::vector<double>& , double, std::vector<double>&);  //Point radius form with Normal to the plane of the circle
    Circle(Sphere&, Plane&);
    Circle(Sphere&, Sphere&);
    
    std::vector<double> get_centre(){return(centre);}
    double get_radius() { return(radius) ;};
    std::vector<double> get_normal(){return(normal);}
    
    double area();
    double intersection_area(Circle);
    double sector_area (double);
    double segment_area (double);
    double segment_perimeter (double);
    double arc_length (double);
    
protected:
    double radius;
    std::vector<double> centre;
    std::vector<double> normal;
    
    
};
#endif /* Circle_hpp */

