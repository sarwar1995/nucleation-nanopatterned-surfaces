//
//  HourglassCluster.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef HourglassCluster_hpp
#define HourglassCluster_hpp

#include <stdio.h>
#include "Shape.hpp"
#include <iostream>

class Hourglass:public Shape
{
public:
    Hourglass();
    ~Hourglass();
    Hourglass(std::vector<Shape*>&);
    
    int isInside (std::vector<double>&);
    int nearSurface(std::vector<double>&, double);
protected:
    std::vector<Shape*> Components ;
};

#endif /* HourglassCluster_hpp */
