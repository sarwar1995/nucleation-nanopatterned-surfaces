//
//  SpheroCylinder.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "SpheroCylinder.hpp"

SpheroCylinder::SpheroCylinder()
{}

SpheroCylinder::~SpheroCylinder ()
{}

SpheroCylinder::SpheroCylinder (Spherical_cap cap, double length, std::vector<double>& normal, double z)
{
    //This constructor assumes that the cutting plane passes through the centre of the sphere.
    z_wall = z;
    spherical_cap = cap;
    double radius_cylinder_base = spherical_cap.sphere_radius();
    std::vector<std::vector<double> > centres_cylinder_base(2);
    std::vector<double> sphere_centre = spherical_cap.sphere_centre();
    std::vector<double> centre_shift_parallel_normal = scalar_mult_to_vector ( (length/2.0), normal) ;
    std::vector<double> centre_shift_antiparallel_normal = scalar_mult_to_vector ( -1*(length/2.0), normal) ;

    centres_cylinder_base[0] = add_vectors(sphere_centre, centre_shift_antiparallel_normal);
    centres_cylinder_base[1] =  add_vectors(sphere_centre, centre_shift_parallel_normal);
    
    Cylinder thecylinder (radius_cylinder_base, centres_cylinder_base);
    cylinder = thecylinder;
}


int SpheroCylinder::isInside (std::vector<double>& point)
{
    double z_coord_point = point[2];
    if(spherical_cap.isInside(point)==1 || (cylinder.isInside(point)==1 && z_coord_point >= z_wall) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
    
}


int SpheroCylinder::nearSurface(std::vector<double>& point, double delta)
{
    if (spherical_cap.nearSurface(point, delta) == 1 || (cylinder.nearSurface(point, delta)==1 && point[2] >= z_wall))
    {
        return 1;
        
    }
    else
    {
        return 0;
        
    }
    
}

