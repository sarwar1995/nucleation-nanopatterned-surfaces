//
//  DynamicBox_parallel.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/28/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef DynamicBox_parallel_hpp
#define DynamicBox_parallel_hpp

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
    
    std::vector<std::vector<double>> get_box(){return box;}
    
    std::vector<int> CheckBoxBreach(Shape* cluster);
    
    void print_box();
    
    
protected:
    MPI_Comm branch_comm;
    int branch_rank, branch_size;
    
    std::vector<std::vector<double>> box;
    double fixed_expansion_length;
    std::vector<double> dimensions;
    //    x_dim, y_dim, z_dim;
};

#endif /* DynamicBox_parallel_hpp */
