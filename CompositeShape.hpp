//
//  CompositeShape.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef CompositeShape_hpp
#define CompositeShape_hpp

#include <stdio.h>
#include "Shape.hpp"

class CompositeShape: public Shape {
    
public:
    CompositeShape();
    ~CompositeShape();
    CompositeShape(std::vector<Shape*>);
    
    int isInside (std::vector<double>&);//Tells whether a point is inside that shape
    int nearSurface(std::vector<double>&, double);
    
    std::vector<double> xy_spread();
    std::vector<double> threeDim_spread();

protected:
    std::vector<Shape*> list_of_shapes;
    
};
#endif /* CompositeShape_hpp */
