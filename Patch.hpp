//
//  Patch.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 9/25/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Patch_hpp
#define Patch_hpp

#include <stdio.h>
#include <vector>
//This is a rectangular patch. Can be made into a base class for more varied patch shapes
class Patch
{
public:
    Patch();
    ~Patch();
    Patch(double, std::vector<double>&, std::vector<double>&);
    double get_x(){return dimensions[0];}
    double get_y(){return dimensions[1];}
    std::vector<double> get_patch_centre() {return patch_centre;}
    std::vector<double> patch_boundaries(); //0:Xmin, 1:Xmax, 2:Ymin, 3:Ymax
    bool isEqual (Patch&);
    bool check_patch_centre(std::vector<double>);
    
protected:
    double contact_angle;
    std::vector<double> patch_centre;
    std::vector<double> dimensions; //0: X-dim 1:Y-dim
    
    
};

#endif /* Patch_hpp */
