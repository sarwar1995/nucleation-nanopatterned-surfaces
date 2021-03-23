//
//  EvolveCluster.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#ifndef EvolveCluster_hpp
#define EvolveCluster_hpp

#include <stdio.h>
#include <vector>
#include "MC.hpp"
#include "CheckBoundary.hpp"
#include "Shape.hpp"
#include "Spherocylinder_cap_composite.hpp"
#include "FreeEnergy.hpp"

class EvolveCluster{
public:
    EvolveCluster();
    ~EvolveCluster();
    EvolveCluster(CheckBoundary* boundary_check,  MC* mc, std::vector<double> volume_SA_mc, std::vector<double> maximum_limits, std::vector<double> increments, double theta_good, double theta_bad, std::vector<double> patch_widths, double z, double Delta, double rho);
    
    
    void EvolveGoodCap ();
    
    void EvolveBadCapWithSpherocylinder (SpheroCylinder* spherocylinder, Shape* Cluster_shape_ptr, double cyl_length, double d_chord_length, double Rg, FILE* outputfile, FILE* output_points_file);
    
    std::vector<Spherical_cap*> get_bad_caps(std::vector<double> c_bad_left, std::vector<double> c_bad_right, double projected_rb);
    
protected:
    CheckBoundary* check_boundary;
    MC* mc_engine;
    std::vector<double> mc_volume_SA;
    double theta_cb, theta_cg;
    double Rb_max, d_Rb, pG_width_x, pG_width_y, z_surface;
    double delta, Rho;
};


#endif /* EvolveCluster_hpp */
