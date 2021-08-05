//
//  Stripes.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Stripes_hpp
#define Stripes_hpp
#include "Surface.hpp"
#include <stdexcept>

#include <stdio.h>
class Stripes:public Surface{
public:
    Stripes();
    Stripes(std::vector<Patch>&, std::vector<std::vector<double> >&, double);
    ~Stripes();
    
    std::vector<int> monitor_cluster_spread(Shape*);
    
    //This function tells whether the stripe bound of stripe with index int i is crossed
    bool is_bounds_crossed(int);
    
    std::vector<int> monitor_box_breach(Shape*); //std::vector<int>
    
    bool surface_bounds_breach (Shape*);
    
    void calc_box(double);
    void initial_box(double);
    inline double box_volume() {return (box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0]) ;}
    
    //Getters
    double get_zwall (){return z_wall;}
    int get_n_unique_patches (){return n_unique_patches;}
    
    void print_all_patch_bounds();
    
private:
    Patch central_patch;    //The patch that has its centre at the origin
    double z_wall;
    int num_patches;
    int n_surrounding_patches;
    int n_patches_per_side_of_good;
    int n_unique_patches;
    
    
    std::vector<int> bounds_crossed;
};


#endif /* Stripes_hpp */
