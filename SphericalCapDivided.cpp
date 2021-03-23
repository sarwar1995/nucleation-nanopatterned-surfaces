//
//  SphericalCapDivided.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/21/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "SphericalCapDivided.hpp"

Spherical_cap_divided::Spherical_cap_divided()
{
    
}
Spherical_cap_divided::~Spherical_cap_divided()
{
    
}

Spherical_cap_divided::Spherical_cap_divided(std::vector<std::vector<double>> centres, double r, std::vector<double> the_normal, double z_surface): divided_cap_centres(centres), radius(r), z_wall(z_surface), normal(the_normal)
{
    bounds_along_normal.resize(2); //Basically the lower and upper bound for the plane that has divided the spherical cap.
    for(int i=0 ; i<(int)the_normal.size(); i++)
    {
        if(the_normal[i] != 0.0){
            normal_direction = i;
            bounds_along_normal[0] = centres[0][normal_direction];
            bounds_along_normal[1] = centres[1][normal_direction];
        }
    }
    
    
    Sphere sphere_1 (radius, centres[0]);
    Sphere sphere_2 (radius, centres[1]);
    sphere_left = sphere_1;
    sphere_right = sphere_2;
}

int Spherical_cap_divided::InsideSphereLeft (std::vector<double>& point)
{
    return (sphere_left.isInside(point) == 1 && point[2] >= z_wall && point[normal_direction] <= bounds_along_normal[0]);
}

int Spherical_cap_divided::InsideSphereRight (std::vector<double>& point)
{
    return (sphere_right.isInside(point) == 1 && point[2] >= z_wall && point[normal_direction] >= bounds_along_normal[1]);
}


int Spherical_cap_divided::isInside (std::vector<double>& point)
{
    
    if (InsideSphereLeft(point) || InsideSphereRight(point))
    {
        return 1;
    }
    else
    {
        return 0;
    }
    
}


int Spherical_cap_divided::nearSurface(std::vector<double>& point, double delta)
{
    if (
        (sphere_left.nearSurface(point, delta) == 1 && point[2] >= z_wall && point[normal_direction] <=  bounds_along_normal[0])
        || (sphere_right.nearSurface(point, delta) == 1 && point[2] >= z_wall && point[normal_direction] >=  bounds_along_normal[1])
        )
    {
        return 1;
        
    }
    else
    {
        return 0;
        
    }
    
    
}
