//
//  Plane.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Plane.hpp"

//P: Point N: Normal
Plane::Plane(std::vector<double>& P, std::vector<double>& N)
{
    point.push_back(P[0]);
    point.push_back(P[1]);
    point.push_back(P[2]);
    
    normal.push_back(N[0]);
    normal.push_back(N[1]);
    normal.push_back(N[2]);
}

double Plane::getConstant()
{
    double c = -1*(normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2]);
    return(c);
}
