//
//  Spherocylinder_cap_composite.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef Spherocylinder_cap_composite_hpp
#define Spherocylinder_cap_composite_hpp

#include <stdio.h>
#include <vector>
#include "SpheroCylinder.hpp"
#include "Spherical_cap.hpp"
//#include "CompositeShape.hpp"


class Spherocylinder_cap_composite: public Shape {
    
public:
    Spherocylinder_cap_composite();
    ~Spherocylinder_cap_composite();
    Spherocylinder_cap_composite(std::vector<Spherical_cap*>, SpheroCylinder*, Surface*);
    
    int isInside (std::vector<double>&);//Tells whether a point is inside that shape
    int nearSurface(std::vector<double>&, double);
    
    std::vector<double> xy_spread();
    std::vector<double> threeDim_spread();
    std::vector<double> projected_SA();
    void oppositeBoundCross ();
    
protected:
    Surface* surface_ptr;
    std::vector<Spherical_cap*> list_of_spherical_caps;
    SpheroCylinder* sphero_cylinder;
    std::vector<Shape*> list_of_shapes;
};


#endif /* Spherocylinder_cap_composite_hpp */
