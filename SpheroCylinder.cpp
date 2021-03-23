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

SpheroCylinder::SpheroCylinder (Spherical_cap cap, double length, std::vector<double> normal_to_the_circle, double z):
z_wall(z), spherical_cap(cap), cylinder_length(length), normal (normal_to_the_circle)
{
    //This constructor assumes that the cutting plane passes through the centre of the sphere.
    double radius_cylinder_base = spherical_cap.sphere_radius();
    std::vector<std::vector<double> > centres_cylinder_base(2);
    std::vector<double> sphere_centre = spherical_cap.sphere_centre();
    
    std::vector<double> centre_shift_parallel_normal = scalar_mult_to_vector ( (length/2.0), normal) ;
    std::vector<double> centre_shift_antiparallel_normal = scalar_mult_to_vector ( -1*(length/2.0), normal) ;
   
    
    centres_cylinder_base[0] = add_double_vectors(sphere_centre, centre_shift_antiparallel_normal);
    centres_cylinder_base[1] =  add_double_vectors(sphere_centre, centre_shift_parallel_normal);
    Cylinder thecylinder (radius_cylinder_base, centres_cylinder_base);
    cylinder = thecylinder;
    
    
    Spherical_cap_divided divided_cap (centres_cylinder_base, radius_cylinder_base, normal, z);
    spherical_cap_divided = divided_cap;
}


int SpheroCylinder::isInside (std::vector<double>& point)
{
    double z_coord_point = point[2];
    if(spherical_cap_divided.isInside(point)==1 || (cylinder.isInside(point)==1 && z_coord_point >= z_wall) )
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
    if (spherical_cap_divided.nearSurface(point, delta) == 1 || (cylinder.nearSurface(point, delta)==1 && point[2] >= z_wall))
    {
        return 1;
        
    }
    else
    {
        return 0;
        
    }
    
}


std::vector<double> SpheroCylinder::threeDim_spread()
{
    std::vector<double> result (5,0.0); //[Xmin, Xmax, Ymin, Ymax, Zmax]
    std::vector<double> spherical_cap_spread = spherical_cap.threeDim_spread();
    int cylinder_expansion_direction = 1; //The default value for cylinder expanding in the y-direction.
    for(int i=0; i<(int) normal.size(); i++)
    {
        if(normal[i] != 0.0)
        {
            cylinder_expansion_direction = i;
            
        }
    }
    
    result = spherical_cap_spread;
    if(cylinder_expansion_direction == 0)
    {
        result[0] = result[0] - (cylinder_length/2.0);
        result[1] = result[1] + (cylinder_length/2.0);
        printf("The normal direction is x!! Should  not be the case\n");
        abort();
        
    }
    if(cylinder_expansion_direction == 1)
    {
        result[2] = result[2] - (cylinder_length/2.0);
        result[3] = result[3] + (cylinder_length/2.0);
        
    }
            
     if(cylinder_expansion_direction == 2)
     {
            printf("The normal direction is z!! Should  not be the case in 3d spread\n");
            printf("normal size = %d\t normal = %10.10f\t%10.10f\t%10.10f\n", (int)normal.size(), normal[0], normal[1], normal[2]);
            abort();
     }
//    printf("After spherocylinder three Dim spread\n");
//    printf("spherocylinder three Dim spread result = [ %10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f ]\n", result[0], result[1], result[2], result[3], result[4]);
    return result;
}


std::vector<double> SpheroCylinder::xy_spread()
{
    std::vector<double> result (4, 0.0); //[Xmin, Xmax, Ymin, Ymax]
    std::vector<double> spherical_cap_xy_spread = spherical_cap.xy_spread();
    int cylinder_expansion_direction = 1; //The default value for cylinder expanding in the y-direction.
    for(int i=0; i<normal.size(); i++)
    {
        if(normal[i] != 0.0) {cylinder_expansion_direction = i;}
    }
    result = spherical_cap_xy_spread;
    if(cylinder_expansion_direction == 0)
    {
        result[0] = result[0] - (cylinder_length/2.0);
        result[1] = result[1] + (cylinder_length/2.0);
        printf("The normal direction is x!! Should  not be the case\n");
        abort();
        
    }
    if(cylinder_expansion_direction == 1)
    {
        result[2] = result[2] - (cylinder_length/2.0);
        result[3] = result[3] + (cylinder_length/2.0);
    }
    if(cylinder_expansion_direction == 2)
    {
        printf("The normal direction is z!! Should  not be the case in xy spread\n");
        printf("normal size = %d\t normal = %10.10f\t%10.10f\t%10.10f\n", (int)normal.size(), normal[0], normal[1], normal[2]);
        abort();
    }
    
    return result;
}


double SpheroCylinder::AnalytVolumeWithPlaneIntersect ()
{
    double volume;
    double cyl_vol = cylinder.AnalytVolumeWithPlaneIntersect(z_wall);
    double sph_vol = spherical_cap.getVolume();
    printf("cyl_vol = %10.10f\t sph_vol = %10.10f\n", cyl_vol, sph_vol);
    volume = cyl_vol + sph_vol;
    //cylinder.AnalytVolumeWithPlaneIntersect(z_wall) +
    return volume;
    
}

double SpheroCylinder::AnalytSurfAreaWithPlaneIntersect ()
{
    double SA;
    SA = cylinder.AnalytSurfAreaWithPlaneIntersect(z_wall) + spherical_cap.getSA();
    return SA;
}

double SpheroCylinder::AnalytProjSurfAreaWithPlaneIntersect ()
{
    double proj_SA;
    proj_SA = cylinder.AnalytProjSurfAreaWithPlaneIntersect(z_wall) + spherical_cap.projected_SA();
    return proj_SA;
}
