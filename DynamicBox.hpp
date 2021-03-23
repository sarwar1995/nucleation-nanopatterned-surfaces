//
//  DynamicBox.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef DynamicBox_hpp
#define DynamicBox_hpp

#include <stdio.h>
#include <vector>
#include "Shape.hpp"

class DynamicBox{
public:
    DynamicBox();
    ~DynamicBox();
    DynamicBox(std::vector<std::vector<double>> existing_box, double length);
    int add_fix_box(int direction_to_expand_in, double density_of_points);
    
    std::vector<std::vector<double>> get_box(){return box;}
    
    std::vector<int> CheckBoxBreach(Shape* cluster);
    
    void print_box();
    
    
protected:
    std::vector<std::vector<double>> box;
    double fixed_expansion_length;
    std::vector<double> dimensions;
//    x_dim, y_dim, z_dim;
};

#endif /* DynamicBox_hpp */
