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
#include <mpi.h>

class DynamicBox{
public:
    DynamicBox();
    ~DynamicBox();
    DynamicBox(std::vector<std::vector<double>> existing_box, double length);
    int add_fix_box(int direction_to_expand_in, double density_of_points);
    int add_virtual_fix_box(int direction_to_expand_in, double density_of_points);
    
    void reset_virtual_box();
    std::vector<std::vector<double>> get_box(){return box;}
    std::vector<std::vector<double>> get_virtual_box(){return virtual_box;}
    
    std::vector<int> CheckBoxBreach(Shape* cluster);
    std::vector<int> CheckVirtualBoxBreach(Shape* cluster);

    
    void updateDimensions();
    void updateVirtualDimensions();
    
    void print_box();
    
    
protected:
    int myRank, nProcs;
    std::vector<std::vector<double>> box;
    std::vector<std::vector<double>> virtual_box; //A box that gets reset to original box. Used for the memory-efficient version of add and loop points function of MC
    double fixed_expansion_length;
    std::vector<double> dimensions;
    std::vector<double> virtual_dimensions;
//    x_dim, y_dim, z_dim;
};

#endif /* DynamicBox_hpp */
