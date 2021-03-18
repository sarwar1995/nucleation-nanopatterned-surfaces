//
//  Spherical_cap.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Spherical_cap.hpp"


Spherical_cap::Spherical_cap()
{
    //Default
}

Spherical_cap::Spherical_cap (Sphere& S, std::vector<double>& proj_centre, double theta, double z)
{
    z_wall = z;
    sphere = S;
    center_proj = proj_centre ;
    theta_c = theta;
    alternate = 0;
}

Spherical_cap::Spherical_cap(std::vector<double>& proj_centre, double proj_radius, double theta, double z)
{
    z_wall = z;
    center_proj = proj_centre ;
    radius_proj = proj_radius;
    theta_c = theta;
    double r = proj_radius/(double)sin(theta);
    std::vector<double> c(3,0.0);
    c[0] = proj_centre[0]; c[1] = proj_centre[1];
    c[2] = proj_centre[2] - r*cos(theta) ; //Here 0.0 is z_surface = 0.0 unless otherwise stated
    Sphere s (r, c);
    sphere = s;
    alternate = 0;
}

Spherical_cap::Spherical_cap(std::vector<double>& proj_centre, double proj_radius, double theta, double z,  int alt)
{
    z_wall = z;
    center_proj = proj_centre ;
    radius_proj = proj_radius;
    theta_c = theta;
    double r = proj_radius/(double)sin(theta);
    std::vector<double> c(3,0.0);
    c[0] = proj_centre[0]; c[1] = proj_centre[1];
    
    alternate = alt;
    printf("alternate = %d\n", alternate);
    /* This alternate is here for the case of hourglass nuclei where the second cap is reversed, hence the z-coordinate of the centre is reversed if alternate is true */
    c[2] = (alternate == 0) ? (proj_centre[2] - r*cos(theta)) : (proj_centre[2] + r*cos(theta)) ;
    Sphere s (r, c);
    sphere = s;
}


Spherical_cap::~Spherical_cap()
{
    
}

Circle Spherical_cap::get_circle()
{
    Circle proj_circle (center_proj , radius_proj);
    return proj_circle;
}

//double Spherical_cap::get_theta_c(string cap_type)
//{
//    double dist = plane_centre_dist();
//    double radius;
//    sphere.getRadius(radius);
//    double theta_c = asin(dist/radius);
//    string small ("small");
//    if(small.compare(cap_type))
//    {return(theta_c);}
//    else{return(pi - theta_c);}
//}

double Spherical_cap::getHeight()
{
    double r = sphere.radius;
    return (r - r*cos(theta_c));
}

double Spherical_cap::getVolume()
{
    double sphere_vol = sphere.volume;
    double result = sphere_vol*potency_factor(theta_c);
    return result;
}

double Spherical_cap::getSA()
{
    double result = 2*pi*sphere.radius*getHeight() ;
    return result;
}

double Spherical_cap::projected_SA()
{
    Circle circle = get_circle();
    return (circle.area());
}

std::vector<double> Spherical_cap::xy_spread()
{
    std::vector<double> result (4, 0.0); //[Xmin, Xmax, Ymin, Ymax]
    result[0] = center_proj[0] - radius_proj;
    result[1] = center_proj[0] + radius_proj;
    result[2] = center_proj[1] - radius_proj;
    result[3] = center_proj[1] + radius_proj;
    return result;
}

std::vector<double> Spherical_cap::threeDim_spread()
{
    /* Only use for non hourglass like composite clusters */
    std::vector<double> result (5, 0.0); //[Xmin, Xmax, Ymin, Ymax, Zmax]
    
    if(theta_c >= (pi/2.0))
    {
        if(!(sphere.centre[2] >= z_wall))
        {throw std::invalid_argument("Invalid sphere centre for obtuse angle in Shperical_cap.cpp"); abort();}
        else
        {
            result[0] = sphere.centre[0] - sphere.radius;
            result[1] = sphere.centre[0] + sphere.radius;
            result[2] = sphere.centre[1] - sphere.radius;
            result[3] = sphere.centre[1] + sphere.radius;
            result[4] = sphere.centre[2] + sphere.radius;
            
        }
    }
    else
    {
        if(!(sphere.centre[2] < z_wall))
        {throw std::invalid_argument("Invalid sphere centre for acute angle in Shperical_cap.cpp"); abort();}
        else
        {
            result[0] = center_proj[0] - radius_proj;
            result[1] = center_proj[0] + radius_proj;
            result[2] = center_proj[1] - radius_proj;
            result[3] = center_proj[1] + radius_proj;
            result[4] = sphere.centre[2] + sphere.radius;
        }
        
    }
    return result; 
}

int Spherical_cap::isInside (std::vector<double>& point)
{
    if(alternate == 0)
    {
        if (sphere.isInside(point) == 1 && point[2] >= z_wall)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (sphere.isInside(point) == 1 && point[2] <= z_wall)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

int Spherical_cap::nearSurface(std::vector<double>& point, double delta){
    if(alternate == 0)
    {
        if (sphere.nearSurface(point, delta) == 1 && (point[2] >= z_wall))
        {return 1;}
        else{return 0;}
    }
    else
    {
        if (sphere.nearSurface(point, delta) == 1 && (point[2] <= z_wall))
        {return 1;}
        else{return 0;}
    }
}

std::vector<double> Spherical_cap::sphere_centre()
{
    std::vector<double>  result(3,0.0);
    result[0] = sphere.centre[0];
    result[1] = sphere.centre[1];
    result[2] = sphere.centre[2];
    return result;
}

//double Spherical_cap::plane_centre_dist()
//{
//    std::vector<double> normal = Wall.getNormal();
//    std::vector<double> point(3,0.0);
//    sphere.getCentre(point);
//    double constant = Wall.getConstant();
//    double dist = point_plane_dist(normal, point, constant);
//    return(abs(dist));
//}

//int Spherical_cap::plane_intersection()
//{
//    int yes;
//    double radius;
//    sphere.getRadius(radius);
//    if(plane_centre_dist() > radius)
//    {
//        yes = 1;
//    }
//    else{yes = 0;}
//    return(yes);
//}

