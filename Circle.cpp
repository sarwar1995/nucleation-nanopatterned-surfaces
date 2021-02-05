//
//  Circle.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/6/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Circle.hpp"
Circle::Circle()
{
    
}

Circle::~Circle()
{
    
}

Circle::Circle(std::vector<double>& Point , double Radius)
{
    radius = Radius;
    centre.resize(Point.size());
    for(size_t i =0 ; i<Point.size(); i++)
    {
        centre[i] = Point[i];
    }
    normal.resize(3);
    normal[0] = 0; normal[1] = 0; normal[2] = 1;
}


double Circle::area()
{
    return pi*radius*radius;
}

double Circle::intersection_area (Circle circle_2)
{
    std::vector<double> circle_2_centre;
    double R2, intersect_area;
    circle_2_centre = circle_2.get_centre();
    R2 = circle_2.get_radius();
    double R2_sq = R2*R2;
    double R1 = radius, R1_sq = R1*R1;
    std::vector<double> d_vec;
    subtract_vectors(centre, circle_2_centre, d_vec);
    double d = vector_norm(d_vec);
    if(d >= R1+R2)
    {return 0.0 ;}
    else
    {
        double d_sq = d*d;
        double h1 = (d_sq - R2_sq + R1_sq)/(2*d);
        double h2 = (d_sq + R2_sq - R1_sq)/(2*d);
        double area_segment_1 = segment_area(h1);
        double area_segment_2 = circle_2.segment_area(h2);
        intersect_area = area_segment_1 + area_segment_2;
        return intersect_area;
        
    }
}

double Circle::sector_area (double theta)
{
    return radius*radius*(theta/2.0) ;
}

double Circle::segment_area (double h)
{
    double Rsquare = radius*radius;
    double theta = 2 * acos(h/radius);
    double area_triangle = h*sqrt(Rsquare - h*h);
    double area_seg = sector_area (theta) - area_triangle;
    return area_seg;
}

//Circle::Circle(Sphere& sphere, Plane& plane)
//{
//    double R;
//    std::vector<double> C(3,0.0);
//    sphere.getRadius(R);
//    sphere.getCentre(C);
//
//    //std::vector<double> normal;
//    double Constant;
//    Constant = plane.getConstant();
//    normal.resize(3);
//    normal = plane.getNormal();
//
//    double dist_plane_centre = point_plane_dist(normal, C,  Constant) ;
//    if(dist_plane_centre < R || -R < dist_plane_centre)
//    {
//        radius = sqrt(R*R - dist_plane_centre*dist_plane_centre);
//        for(size_t i =0 ; i<3; i++)
//        {
//            centre[i] = C[i] - normal[i]*dist_plane_centre ;
//        }
//    }
//    else
//    {
//        radius = -1.0;
//        centre[0] = 0.0; centre[1] = 0.0; centre[2] = 0.0;
//        printf("No circle of intersection possible");
//    }
//
//}
//
//Circle::Circle(Sphere& sphere_1, Sphere& sphere_2)
//{
//    double R1, R2;
//    std::vector<double> C1, C2;
//    //std::vector<double> normal_12(3,0.0);
//    sphere_1.getRadius(R1); sphere_2.getRadius(R2);
//    sphere_1.getCentre(C1); sphere_2.getCentre(C2);
//    std::vector<double> dist_12;
//    subtract_vectors(C1, C2, dist_12);
//    double d12 = vector_norm(dist_12) ;
//    double R_int_12 = radius_of_intersection(R1, R2, d12);
//
//    double x_first_two = x_intersection (R1, R2, d12) ;
//    normal.resize(3);
//    for(size_t i = 0; i<3; i++)
//    {
//        normal[i] = dist_12[i]/d12 ;
//        centre[i] = C1[i] + x_first_two*normal[i];
//    }
//    radius = R_int_12;
//}


