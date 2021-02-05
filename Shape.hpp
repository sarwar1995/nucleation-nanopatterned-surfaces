//
//  Shape.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Shape_hpp
#define Shape_hpp

#include <stdio.h>
#include <vector>


class Shape{
public:
    Shape();
    ~Shape();
    virtual int isInside (std::vector<double>&){return 0;} //Tells whether a point is inside that shape
    virtual int nearSurface(std::vector<double>&, double){return 0;}
    virtual std::vector<double> xy_spread()
    {std::vector<double> zero_vec(4,0.0); return zero_vec;}
    
    virtual std::vector<double> threeDim_spread() {std::vector<double> zero_vec(5,0.0); return zero_vec;}
    
    
};
#endif /* Shape_hpp */
