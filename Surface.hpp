//
//  Surface.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Surface_hpp
#define Surface_hpp

#include <stdio.h>
#include "Patch.hpp"
#include "Shape.hpp"

class Surface {
public:
    Surface();
    Surface(std::vector<Patch>&, std::vector<std::vector<double> >&);
    ~Surface();
    //virtual void evolve_cluster();                  //Grows cluster starting at the central good patch
    //virtual int valid_point (std::vector<double>&) = 0;  //A valid point above the surface
    
    
    virtual std::vector<int> monitor_cluster_spread(Shape*){std::vector<int> zero_vec(1,-1); return zero_vec;}
    virtual void calc_box(double){};
    std::vector<std::vector<double> > box;
    
    std::vector<Patch> Patches;   //List of all patches. [0]:central patch, [1]:left bad patch, [2]:right bad patch
    //[3]: left left good patch [4]: right right good patch and so on...
    
protected:
    
    std::vector<std::vector<double> > orientations; //List of patch centres
    
};




#endif /* Surface_hpp */
