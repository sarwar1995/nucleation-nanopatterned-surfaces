//
//  MC.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 8/3/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef MC_hpp
#define MC_hpp

#include <stdio.h>
#include <vector>
#include <set>
#include <random>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include "CellList.hpp"
#include "Shape.hpp"
#include "miscformulas.hpp"
#include <stdexcept>

class MC {
public:
    MC();
    ~MC();
    MC(int, std::vector<std::vector<double> >&, int[3]);
    void generate_points();
    inline double box_volume() {return BoxVolume;}
    
    //Sets of Interior and surface points: Using sets to improve uniqueness checking. For the continuos update of volume and surface area by just checking the surface points. Only works for growth of the same type of cluster.
    std::set<std::vector<double> > interior_points;
    std::set<std::vector<double> > surface_points;
    
    
    //Volume and SA calculation
    std::vector<double> calc_volume_SA(Shape*, double); //Returns volume and surface area
    std::vector<double> update_volume_SA(Shape*, CellList*, double);
    
    //Identify interior and surface points
    void IdentifyPoints (Shape*, double);

    //Miscellaneous I/O
    void print_points(FILE*);
    void print_surf_points(FILE*);
    
    //Getters
    std::vector<std::vector<double> > get_points() {return x_points;}
    std::vector<std::vector<double> > get_box() {return box;}
    
    
private:
    int seed[3] ;
    int n_points;
    int n_inside, n_near_surf;
    std::vector<std::vector<double> > box;
    std::vector<std::vector<double> > x_points;
    double BoxVolume;
    
    //Volume and SA given interior and surface points
    std::vector<double> getMeasures (int, int, double);

};

#endif /* MC_hpp */
