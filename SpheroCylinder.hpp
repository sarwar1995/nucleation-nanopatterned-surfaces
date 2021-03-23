//
//  SpheroCylinder.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef SpheroCylinder_hpp
#define SpheroCylinder_hpp

#include <stdio.h>
#include "Spherical_cap.hpp"
#include "Circle.hpp"
#include "Cylinder.hpp"
#include "miscformulas.hpp"
#include "SphericalCapDivided.hpp"

class SpheroCylinder: public Shape
{
public:
    //Constructors and Destructors
    SpheroCylinder ();
    ~SpheroCylinder ();
    SpheroCylinder (Spherical_cap, double, std::vector<double>, double); //A spherical cap on the xy plane, the length of the cylinder in between and the normal/axis of the cylinder or the normal to the cutting plane
    
    //Concrete methods from Shape class
    int isInside (std::vector<double>&); //Tells whether a point is inside that shape
    int nearSurface(std::vector<double>&, double);
    std::vector<double> xy_spread();
    std::vector<double> threeDim_spread();
    
    std::vector<double> get_normal() {return normal;}
    
    //These three functions are for cylinder with axis parallel to the xy plane.
    double AnalytVolumeWithPlaneIntersect (); //plane parallel to xy axis i.e. z=constant
    double AnalytSurfAreaWithPlaneIntersect ();
    double AnalytProjSurfAreaWithPlaneIntersect ();
    
    
protected:
    Spherical_cap spherical_cap;
    Spherical_cap_divided spherical_cap_divided;
    Cylinder cylinder;
    double cylinder_length;
    std::vector<double> normal;
    double z_wall;
};


#endif /* SpheroCylinder_hpp */
