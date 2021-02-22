//
//  Spherical_cap.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Spherical_cap_hpp
#define Spherical_cap_hpp

#include <stdio.h>
#include "Sphere.hpp"
#include "Shape.hpp"
#include "Circle.hpp"
#include <string>
#include "Surface.hpp"
#include "miscformulas.hpp"

using namespace std;
class Spherical_cap:public Shape
{
public:
    Spherical_cap();
    ~Spherical_cap();
    //Spherical_cap(Sphere&);
    Spherical_cap(Sphere&, std::vector<double>&, double, double);
    Spherical_cap(std::vector<double>&, double, double, double);
    Spherical_cap(std::vector<double>&, double, double, double, int);
    
    //Functions that determine a points location
    int isInside (std::vector<double>&); //Any point is inside the spherical cap
    int nearSurface(std::vector<double>&, double);
 
    inline double sphere_radius(){return (sphere.radius);};
    std::vector<double> sphere_centre();
    //void print_radius_proj(){printf("radius_proj = %10.10f\n",radius_proj);};
    
    double getVolume();
    double getSA();
    double getHeight();
    double projected_SA();
    std::vector<double> xy_spread();
    std::vector<double> threeDim_spread();
        
    //The sphere
    Sphere sphere;
    Circle get_circle();
    
protected:
    double z_wall;
    //Surface surface;
    double radius_proj;
    std::vector<double> center_proj;
    double theta_c;
    int alternate;
};
#endif /* Spherical_cap_hpp */
