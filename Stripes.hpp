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
    
    std::vector<int> monitor_box_breach(Shape*);
    
    //void evolve_cluster(double, double);
    //int valid_point(std::vector<double>&);
    void calc_box(double);
    std::vector<std::vector<double> > box;
    inline double box_volume() {return (box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0]) ;}
private:
    Patch central_patch;    //The patch that has its centre at the origin
    double z_wall;
};


#endif /* Stripes_hpp */
