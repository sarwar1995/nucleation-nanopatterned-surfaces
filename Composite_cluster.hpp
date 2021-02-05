//
//  Composite_cluster.hpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#ifndef Composite_cluster_hpp
#define Composite_cluster_hpp

#include <stdio.h>
#include "Sphere.hpp"
#include "Spherical_cap.hpp"
#include "Circle.hpp"
#include "Surface.hpp"

class Composite_cluster:public Shape {
public:
    Composite_cluster();
    ~Composite_cluster();
    
    Composite_cluster(std::vector<Spherical_cap>&, Surface&); //vector of planes is boundaries. Plane is wall
    
    int isInside(std::vector<double>&) ; //Checks whether a point lies inside the cluster or not
    int nearSurface(std::vector<double>&, double); //Checks if a point is whithin a region of width 2*delta of the surface
    inline int num_caps() {return (int)list_of_spherical_caps.size();}
    bool similar_shapes (Composite_cluster&, Composite_cluster&);
    
    std::vector<double> xy_spread(); //Spread of the cluster in x-y plane
    std::vector<double> threeDim_spread(); //Maximum spread of the cluster in the box to check for box breach
    
    std::vector<double> projected_SAs();
    
//    void calc_volume(MC&);
//    void calc_SA(MC&);
    
    //void calc_intersections();
    
//    std::vector<double> get_two_intersect_measures(Sphere&, Sphere&, int , int, Plane&, double, double);
//    //part of B1 B2 intersection that is outside G.
//    //Not the actual measures of the three intersection region above the wall
//    std::vector<double> get_three_intersect_measures();
//    std::vector<double> get_net_measures();
//    std::vector<double> get_proj_SA();      //[0]=good (centre) {[1]=bad 1 [2] = bad 2}surrounding two
//
    
protected:
    std::vector<Spherical_cap> list_of_spherical_caps; //[Good, Bad left, Bad right]
    Surface surface;
    
    std::vector<double> Cb1, Cg, Cb2;               // Sphere centres
    double Rb1, Rg, Rb2 ;                           // Sphere radii
    double d_B1_G, d_B2_G, d_B1_B2;                 // Distance between centres
    double R_int_B1_G, R_int_B2_G, R_int_B1_B2;     // Radius of intersection between spheres
    std::vector<double> vec_B1G, vec_B2G, vec_B1B2; // Vector joining centres
    std::vector<double> normal_B1G, normal_B2G, normal_B1B2;  //Normals joining centres
    
    std::vector<std::vector<int> >intersections;    // Matrix of ints showing pairwise intersections, -1: self, 1: intersects, 0: no intersection
    int three_way_intersect;                        // Tell wheter three way interesct is happening
    int three_way_type;                             //0: completely inside good, 1:partially inside
    
    
};

#endif /* Composite_cluster_hpp */
