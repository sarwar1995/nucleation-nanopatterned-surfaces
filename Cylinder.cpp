//
//  TestHourglass.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "Cylinder.hpp"

Cylinder::Cylinder()
{
    
}

Cylinder::~Cylinder()
{
    
}

Cylinder::Cylinder(double radius, std::vector<std::vector<double> >& circle_centres)
{
    rc = radius;
    if((int)circle_centres.size() != 2)
    {
        throw std::invalid_argument("Size of centres should be 2");
    }
    else
    {
        centres = circle_centres;
        
    }
    normalized_axis.resize(3);
    for (size_t i=0; i<centres[0].size(); i++)
    {
        normalized_axis[i] = centres[1][i] - centres[0][i];
    }
    
    length = vector_norm(normalized_axis);
    for (size_t i=0; i<centres[0].size(); i++)
    {
        normalized_axis[i] = normalized_axis[i]/length;
    }
    
    Circle a_base_circle (centres[0], rc, normalized_axis); //The normalized_axis of the cylinder is the circle's normal
    Circle a_top_circle (centres[1], rc, normalized_axis);
    base_circle = a_base_circle;
    top_circle = a_top_circle;
    
}

int Cylinder::isWithin(std::vector<double> & point)
{
    int yes = 0;
    std::vector<double> pointMinusP1 (3, 0.0);
    std::vector<double> pointMinusP2 (3, 0.0);
    std::vector<double> Axis (3, 0.0);
    for (size_t i=0; i<point.size(); i++)
    {
        pointMinusP1[i] = point[i] - centres[0][i];
        pointMinusP2[i] = point[i] - centres[1][i];
        //Axis[i] = centres[1][i] - centres[0][i];
    }
    double dotProd_with_p1 = dotprod(pointMinusP1, normalized_axis);
    double dotProd_with_p2 = dotprod(pointMinusP2, normalized_axis);
    if(dotProd_with_p1 >= 0 && dotProd_with_p2 <= 0)
    {yes = 1;}
    else{yes = 0;}
    return yes;
}


int Cylinder::nearSurface(std::vector<double> & point,  double delta)
{
    std::vector<double> centre_near_wall = centres[0]; //p1
    std::vector<double> centre_near_free = centres[1]; //p2
    double point_axis_distance = point_line_dist(point, centre_near_wall, centre_near_free);
    int pointWithin = isWithin(point);
    if((pointWithin == 1) && (point_axis_distance <= rc + delta) && (rc - delta <= point_axis_distance))
    {
        return 1;
    }
    else{return 0;}    
}

int Cylinder::isInside(std::vector<double> & point)
{
    int result = 0;
    std::vector<double> centre_near_wall = centres[0]; //p1
    std::vector<double> centre_near_free = centres[1]; //p2
    
    double point_axis_distance = point_line_dist(point, centre_near_wall, centre_near_free);
    int pointWithin = isWithin(point);
    if(pointWithin == 1 && point_axis_distance < rc )
    {
        result = 1;
    }
    else{result = 0;}
    return result ;
    
}


double Cylinder::AnalytVolumeWithPlaneIntersect(double z)
{
    //Here the height from the circle centre to the plane (intersecting chord) must be positive for acute angled caps and negative for obtuse ones. i.e. if circle centre is below z than h should be positive.
    double volume;
    std::vector<double> circle_centre = base_circle.get_centre();
    double h = z - circle_centre[2];
    //centres[0][2];
    double segment_area = base_circle.segment_area(h);
    volume = segment_area * length;
    return volume;
}

double Cylinder::AnalytSurfAreaWithPlaneIntersect(double z)
{
    double SA;
    double h = z - centres[0][2];
    double segment_perimeter = base_circle.segment_perimeter(h);
    SA = segment_perimeter * length;
    return SA;
}

double Cylinder::AnalytProjSurfAreaWithPlaneIntersect(double z)
{
    double proj_sa;
    double h = z - centres[0][2];
    double theta = acos(h/rc); //Here rc is the radius of the base circles which is the radius of the spherical cap that the circles belong to
    proj_sa = 2.0 * rc * sin(theta) * length ;
    return proj_sa;
}

