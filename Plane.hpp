//
//  Plane.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Plane_hpp
#define Plane_hpp

#include <stdio.h>
#include <vector>
#include "Shape.hpp"

class Plane{

public:
    Plane();
    ~Plane();
    Plane(std::vector<double>&, std::vector<double>&) ;
    
    std::vector<double> getNormal(){return(normal);} ;
    double getConstant() ;
    
protected:
    std::vector<double> point ;
    std::vector<double> normal ;
    
};
#endif /* Plane_hpp */
