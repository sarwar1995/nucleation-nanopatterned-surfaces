//
//  Composite_cluster.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include "Composite_cluster.hpp"
#include "Integrals.hpp"


Composite_cluster::Composite_cluster(){
    
}

Composite_cluster::~Composite_cluster(){
    
}

Composite_cluster::Composite_cluster(std::vector<Spherical_cap>& capsList, Surface& surf) //std::vector<Plane>& boundaries,d
{
 
    list_of_spherical_caps = capsList;
    surface = surf;
    //0:Good , 1:Bad left 1, 2:Bad right
    Spherical_cap Cap_B1 = capsList[1];
    Spherical_cap Cap_G = capsList[0];
    Spherical_cap Cap_B2 = capsList[2];
    
//    sphere_B1.getCentre(Cb1);
//    sphere_G.getCentre(Cg);
//    sphere_B2.getCentre(Cb2);
//
//    sphere_B1.getRadius(Rb1);
//    sphere_G.getRadius(Rg);
//    sphere_B2.getRadius(Rb2);
//
//    subtract_vectors(Cb1, Cg, vec_B1G);        //B1 to G
//    subtract_vectors(Cg, Cb2, vec_B2G);        //G to B2
//    subtract_vectors(Cb1, Cb2, vec_B1B2);      //B1 to B2
//
//
//    d_B1_G = sqrt(vec_B1G[0]*vec_B1G[0] + vec_B1G[1]*vec_B1G[1] + vec_B1G[2]*vec_B1G[2]);
//    d_B2_G = sqrt(vec_B2G[0]*vec_B2G[0] + vec_B2G[1]*vec_B2G[1] + vec_B2G[2]*vec_B2G[2]);
//    d_B1_B2 = sqrt(vec_B1B2[0]*vec_B1B2[0] + vec_B1B2[1]*vec_B1B2[1] + vec_B1B2[2]*vec_B1B2[2]);
//
//
//    R_int_B1_G = radius_of_intersection(Rb1, Rg, d_B1_G); // of two spheres
//    R_int_B2_G = radius_of_intersection(Rg, Rb2, d_B2_G);
//    R_int_B1_B2 = radius_of_intersection(Rb1, Rb2, d_B1_B2);
//
//    for(size_t i=0; i<3; i++)
//    {
//        normal_B1G[i] = vec_B1G[i]/d_B1_G ;
//        normal_B2G[i] = vec_B2G[i]/d_B2_G ;
//        normal_B1B2[i] = vec_B1B2[i]/d_B1_B2 ;
//    }
//    std::vector<double> plane_point, plane_normal;
//    double plane_constant;
//    plane_normal = wall.getNormal();
//    plane_constant = wall.getConstant();
//    plane_point[0] = -(plane_constant/plane_normal[0]) ; plane_point[1] = 0.0; plane_point[2] = 0.0;
//    Plane Wall(plane_point, plane_normal);
}

// Check whether a point lies inside the boudaries of the cluster
int Composite_cluster::isInside(std::vector<double>& point)
{
    //Write a test to check this i.e. if the point is inside any of the caps then it is inside the composite cluster
    for (size_t i=0; i<list_of_spherical_caps.size(); i++)
    {
        if (list_of_spherical_caps[i].isInside(point) == 1)
        {
            return 1;
        }
        else{continue;}
    }
    return 0;
}

int Composite_cluster::nearSurface(std::vector<double> & point, double delta)
{
    for (size_t i=0; i<list_of_spherical_caps.size(); i++)
    {
        Spherical_cap cap_i = list_of_spherical_caps[i];
        if (cap_i.nearSurface(point, delta) == 1)
        {
            int outOfOthers = 1;
            for (size_t j=0; j<list_of_spherical_caps.size(); j++)
            {
                if(i != j)
                {
                    Spherical_cap cap_j = list_of_spherical_caps[j];
                    if (cap_j.isInside(point) == 1){outOfOthers = 0;}
                    else{continue;}
                }
                else{continue;}
                
            }
            if(outOfOthers==1)
            {return 1;}
            else{return 0;} //continue;
        }
        else{continue;}
    }
    return 0;  //If it gets out of the for loop without returning 1 within it then returns 0
}

//Here defined as clusters having the same number of constituent spherical caps
bool Composite_cluster::similar_shapes (Composite_cluster& cluster_1, Composite_cluster& cluster_2){
    return (cluster_1.num_caps() == cluster_2.num_caps());
}


std::vector<double> Composite_cluster::xy_spread()
{
    std::vector<double> thisSpread (4,0.0);
    std::vector<double> spread_good = list_of_spherical_caps[0].xy_spread();
    std::vector<double> spread_bad_left = list_of_spherical_caps[1].xy_spread();
    std::vector<double> spread_bad_right = list_of_spherical_caps[2].xy_spread();
  
    thisSpread[0] = std::min({spread_good[0],spread_bad_left[0], spread_bad_right[0]});
    thisSpread[1] = std::max({spread_good[1],spread_bad_left[1], spread_bad_right[1]});
    thisSpread[2] = std::min({spread_good[2],spread_bad_left[2], spread_bad_right[2]});
    thisSpread[3] = std::max({spread_good[3],spread_bad_left[3], spread_bad_right[3]});
    return thisSpread;
}

std::vector<double> Composite_cluster::threeDim_spread()
{
    //To check the breach of 5 box planes, excluding the wall surface
    //[Xmin, Xmax, Ymin, Ymax, Zmax]
    //Start with central caps spread as the spread of cluster and then modify based on other caps
    std::vector<double> thisSpread = list_of_spherical_caps[0].threeDim_spread();
    
    for(size_t i=1; i<list_of_spherical_caps.size(); i++)
    {
        std::vector<double> iCap3DSpread = list_of_spherical_caps[i].threeDim_spread();
        if(iCap3DSpread[0] <= thisSpread[0]) { thisSpread[0] = iCap3DSpread[0]; }
        if(iCap3DSpread[1] >= thisSpread[1]) { thisSpread[1] = iCap3DSpread[1]; }
        if(iCap3DSpread[2] <= thisSpread[2]) { thisSpread[2] = iCap3DSpread[2]; }
        if(iCap3DSpread[3] >= thisSpread[2]) { thisSpread[3] = iCap3DSpread[3]; }
        if(iCap3DSpread[4] >= thisSpread[4]) { thisSpread[4] = iCap3DSpread[4]; }
    }
    return thisSpread;
}


std::vector<double> Composite_cluster::projected_SAs()
{
    std::vector<double> projected_areas;
    std::vector<Circle> circles(list_of_spherical_caps.size());
    for(size_t i=0; i<circles.size(); i++)
    {
        circles[i] = list_of_spherical_caps[i].get_circle();
    }
    
    /* For central good circle's projected area*/
    double patch_x = surface.Patches[0].get_x();
    double h_g_onb1 = (patch_x/2.0); //Height is the distance of chord from centre
    double h_g_onb2 = h_g_onb1; //Symmetric
    double segment_area_G_onB1 = circles[0].segment_area(h_g_onb1);
    double segment_area_G_onB2 = circles[0].segment_area(h_g_onb2);
    double proj_area_good = circles[0].area() - segment_area_G_onB1 - segment_area_G_onB2;
    
    /* For left bad circle's projected area*/
    std::vector<double> c_b1 = circles[1].get_centre();
    double h_b1_ong = (-patch_x/2.0) - c_b1[0];
    double segment_area_B1_onG = circles[1].segment_area(h_b1_ong);
    double proj_area_bad_1 = (circles[1].area() - segment_area_B1_onG);
    
    /* For right bad circle's projected area*/
    std::vector<double> c_b2 = circles[2].get_centre();
    double h_b2_ong = c_b2[0] - (patch_x/2.0);
    double segment_area_B2_onG = circles[2].segment_area(h_b2_ong);
    double proj_area_bad_2 =  (circles[2].area() - segment_area_B2_onG);
    
//    double lens_good_badleft = circles[0].intersection_area(circles[1]);
//    double lens_good_badright = circles[0].intersection_area(circles[2]);
    projected_areas.push_back(proj_area_good);
    projected_areas.push_back(proj_area_bad_1);
    projected_areas.push_back(proj_area_bad_2);
    //    return (area_good + area_bad_left + area_bad_right - lens_good_badleft - lens_good_badright);

    return projected_areas;
}


//void Composite_cluster::calc_intersections()
//{
//    intersections.resize(3);
//    //intersections[0][i] : Good with others
//    //intersections[1][i] : Bad 1 with others
//    //intersections[2][i] : Bad 2 with others
//
//    for(size_t i=0; i<3; i++)
//    {
//        intersections[i].resize(3);
//        for(size_t j=0; j<3; j++)
//        {
//            if(i==j){intersections[i][j] = -1 ;}
//
//        }
//    }
//
//    if(isnan(R_int_B1_G)){intersections[0][1] = 0; intersections[1][0] = 0;}
//    else{ intersections[0][1] = 1; intersections[1][0] = 1; }
//
//    if(isnan(R_int_B2_G)){intersections[0][2] = 0; intersections[2][0] = 0;}
//    else{intersections[0][2] = 1; intersections[2][0] = 1;}
//
//    if(isnan(R_int_B1_B2)){intersections[1][2] = 0; intersections[2][1] = 0;}
//    else{intersections[1][2] = 1; intersections[2][1] = 1;}
//
//    //Three_way_intersect: only if all three at least intersect in pairs
//    if(intersections[0][1] && intersections[0][2] && intersections[1][2])
//    {
//        //Intersection of Rb1 and Rb2
//        double x_first_two = x_intersection (Rb1, Rb2, d_B1_B2) ;
//        std::vector<double> Pb1b2(3,0.0);       //Point of intersection of B1 and B2. Also the centre of the intersected circle
//
//        for(size_t i = 0; i<3; i++)
//        {
//            Pb1b2[i] = Cb1[i] + x_first_two*normal_B1B2[i];
//        }
//
//        Circle circle_B1B2 (Pb1b2, R_int_B1_B2);       //Circle of intersection of B1 and B2
//        Plane plane_B1B2 (Pb1b2, normal_B1B2);         //Plane of intersection of B1 and B2
//        Circle circle_G_B1B2 (sphere_G, plane_B1B2) ;  //Circle of intersection of G and Plane of B1 and B2
//
//        if(circle_G_B1B2.get_radius() == -1)           //G does not intersect with Plane B1B2
//        {three_way_intersect = 0;}
//        else
//        {
//            double R = circle_B1B2.get_radius(), r = circle_G_B1B2.get_radius();
//            std::vector<double> centre_2 =  circle_G_B1B2.get_centre();
//            double chord_Cb1b2_Cg = circle_circle_int_radii(R, r, Pb1b2, centre_2);
//
//            std::vector<double> dist_Cb1b2_Cg(3,0.0);
//            subtract_vectors(Pb1b2, centre_2, dist_Cb1b2_Cg);
//            double d_Cb1b2_Cg = sqrt(dist_Cb1b2_Cg[0]*dist_Cb1b2_Cg[0] + dist_Cb1b2_Cg[1]*dist_Cb1b2_Cg[1] + dist_Cb1b2_Cg[2]*dist_Cb1b2_Cg[2]);    //distance between the two centres on the plane of intersection of B1 and B2
//
//
//            if(isnan(chord_Cb1b2_Cg))                           //Completely inside or outside
//            {
//                if(d_Cb1b2_Cg < (std::max(R,r)- std::min(R,r))) //Completely inside
//                {
//                    three_way_intersect = 1;
//                    three_way_type = 0;
//                }
//                else
//                {three_way_intersect = 0;}                  //Completely outside
//
//            }
//            else
//            {
//                if(chord_Cb1b2_Cg == 0.0)                       //touching or same circle
//                {
//                    if(d_Cb1b2_Cg == std::max(R,r)- std::min(R,r)) //Inside with touching
//                    {
//                        three_way_intersect = 1;
//                        three_way_type = 0;
//                    }
//                    else{three_way_intersect = 0;}
//                }
//                else{three_way_intersect = 1; three_way_type = 1;}
//            }
//        }
//    }
//    else
//    {
//        three_way_intersect = 0;
//    }
//}


//std::vector<double> Composite_cluster::get_two_intersect_measures(Sphere& sphere_1, Sphere& sphere_2, int ind_1, int ind_2, Plane& wall, double theta_c1, double theta_c2)
//{
//    double theta_11, theta_12, theta_21, theta_22;  //Angles with X axis
//    double phi_11, phi_12, phi_21, phi_22 ;         //Angles with Z axis
//    double d1, d2;
//
//    std::vector<double> Measures (2,0.0);       //[Volume , Surface Area]
//    if(intersections[ind_1][ind_2] == 0)
//    {
//        Measures[0] = 0.0 ;  Measures[1] = 0.0 ;
//    }
//    else
//    {
//        Circle circle_12 (sphere_1, sphere_2);                     //Circle of intersection of sphere 1 and 2
//        Circle circle_proj_1 (sphere_1, wall);                     //Projected circle of intersection of sphere 1 & wall
//        Circle circle_proj_2 (sphere_2, wall);                     //Projected circle of intersection of sphere 2 & wall
//        std::vector<double> centre_circle = circle_12.get_centre();
//        std::vector<double> circle_normal = circle_12.get_normal();
//        double circle_radius = circle_12.get_radius();
//
//        Plane plane_of_intersection (centre_circle, circle_normal);
//
//        std::vector<std::vector<double> > points_XZ_plane (2);  //Top and bottom point of circle in XZ plane
//        points_XZ_plane[0].resize(3);
//        points_XZ_plane[1].resize(3);
//
//        std::vector<double> perp_normal_12 (3,0.0);          //Normal perpendicular to the circle normal. In XZ plane
//        perp_normal_12[0] = -circle_normal[2];              //x component of circle normal is z component of perp
//        perp_normal_12[1] = 0.0;
//        perp_normal_12[2] = circle_normal[0];
//
//
//        points_XZ_plane[0][0] = centre_circle[0] + (circle_radius)*perp_normal_12[0];
//        points_XZ_plane[0][1] = 0.0;
//        points_XZ_plane[0][2] = centre_circle[2] + (circle_radius)*perp_normal_12[2];
//
//        points_XZ_plane[1][0] = centre_circle[0] - (circle_radius)*perp_normal_12[0];
//        points_XZ_plane[1][1] = 0.0;
//        points_XZ_plane[1][2] = centre_circle[2] - (circle_radius)*perp_normal_12[2];
//
//        if(points_XZ_plane[0][2] > wall_z && points_XZ_plane[1][2] > wall_z)
//        {
//            Measures[0] = volume_lens();
//            Measures[1] = surface_area_lens();
//
//        }
//        else
//        {
//            if(points_XZ_plane[0][2] < wall_z && points_XZ_plane[1][2] < wall_z)
//            {
//                Measures[0] = 0.0 ;  Measures[1] = 0.0 ;
//            }
//            else
//            {
//                std::vector<double> top_z = points_XZ_plane[0][2] > wall_z ? points_XZ_plane[0] : points_XZ_plane[1] ;
//                std::vector<double> centre_1, centre_2, d_centre;
//                double radius_1, radius_2;
//                std::vector<double> top_z_minus_c_1, top_z_minus_c_2;
//                sphere_1.getCentre(centre_1); sphere_1.getRadius(radius_1);
//                sphere_2.getCentre(centre_2); sphere_1.getRadius(radius_2);
//
//                subtract_vectors(centre_1, top_z, top_z_minus_c_1);
//                subtract_vectors(centre_2, top_z, top_z_minus_c_2);
//                subtract_vectors(centre_1, centre_2, d_centre);
//                double dist_centre = vector_norm(d_centre);
//                double x_int_spheres = x_intersection (radius_1,radius_2, dist_centre);
//                d1 = abs(x_int_spheres);                            //lower limit for R integration
//                d2 = abs(dist_centre - x_int_spheres);
//
//                std::vector<double> z_axis (3,0.0);
//                z_axis[2] = 1;
//                double dot_z_1 = dotprod(z_axis, top_z_minus_c_1);
//                double dot_z_2 = dotprod(z_axis, top_z_minus_c_2);
//
//                if(dist_centre < 0){printf("dist centres is less than zero"); abort();}
//
//                phi_11 = acos(dot_z_1/(vector_norm(top_z_minus_c_1))) ;
//                phi_12 = theta_c1;
//                phi_21 = acos(dot_z_2/(vector_norm(top_z_minus_c_2))) ;
//                phi_22 = theta_c2;
//
//                std::vector<double> wall_normal = wall.getNormal();
//                double constant_1 = plane_of_intersection.getConstant();
//                double constant_wall = wall.getConstant();
//                std::vector<double> three_plane_int_point = three_plane_intersection(circle_normal, wall_normal, constant_1, constant_wall);
//
//                if(0<x_int_spheres && x_int_spheres<dist_centre)
//                {
//
//                }
//                else
//                {
//                    if(x_int_spheres < 0)
//                    {
//                        phi_11 = acos(dot_z_1/(vector_norm(top_z_minus_c_1))) ;
//                        phi_12 = -1*theta_c1;
//                    }
//                    else if (x_int_spheres > dist_centre)
//                    {
//
//                    }
//
//                }
//
//                //Projected circles on the XY plane
//                std::vector<double> dist_proj_centre, circle_centre_1, circle_centre_2;
//                circle_centre_1 = circle_proj_1.get_centre(); circle_centre_2 = circle_proj_2.get_centre();
//                subtract_vectors(circle_centre_1, circle_centre_2, dist_proj_centre);
//                double d_proj_centre = vector_norm(dist_proj_centre);
//                double proj_circle_intersection = x_intersection (circle_proj_1.get_radius(), circle_proj_2.get_radius(), d_proj_centre);
//
//                if((0<proj_circle_intersection && proj_circle_intersection<d_proj_centre) || (d_proj_centre<proj_circle_intersection && proj_circle_intersection<0))
//                {}
//                else
//                {}
//
//                double d_int_from_centre_1 = proj_circle_intersection;                      //which direction is the intersect from centre 1
//                double d_int_from_centre_2 = -(d_proj_centre - proj_circle_intersection);   //which direction is the intersect from centre 2
//                theta_11 = acos(d_int_from_centre_1/circle_proj_1.get_radius());
//                theta_12 = 2*pi - theta_11;
//                theta_21 = acos(d_int_from_centre_2/circle_proj_2.get_radius());
//                theta_22 = 2*pi - theta_21;
//
////                double partial_volume_1 = Volume_integral (radius_1, d1, theta_11, theta_12,  phi_11,  phi_12);
////                double partial_volume_2 = Volume_integral (radius_1, d1, theta_11, theta_12,  phi_11,  phi_12);
//                Measures[0];
//
//            }
//        }
//
//
//    }
//
//}
//

//std::vector<double> Composite_cluster::get_three_intersect_measures()
//{
//    std::vector<double> three_way_measures(2,0.0);
//    if(three_way_intersect)
//    {
//        if(three_way_type == 0)         //completely inside meaning no volume outside of good
//        {
//            return(three_way_measures);
//        }
//        else
//        {
//
//        }
//    }
//    else
//    {
//        return(three_way_measures);
//    }
//}


//std::vector<double> Composite_cluster::get_net_measures()
//{
//    std::vector<double> net_measures(2,0.0);               //[0]=volume [1] = surface area
//    double V_b1, V_g , V_b2;                        //single cap volumes and surface areas
//    double V_b1_b2, V_g_b1 , V_g_b2, V_g_b1_b2;     //double and triple intersection volumes and surface areas
//    double S_b1, S_g , S_b2;
//    double S_b1_b2, S_g_b1 , S_g_b2, S_g_b1_b2;
//
//    V_b1 = Cap_B1.getVolume();  S_b1 = Cap_B1.getSA();
//    V_b2 = Cap_B2.getVolume();  S_b2 = Cap_B2.getSA();
//    V_g = Cap_G.getVolume();    S_g = Cap_G.getSA();
//
//    double theta_c1, theta_c2;
//
//    std::vector<double> measures_b1_b2 = get_two_intersect_measures(sphere_B1, sphere_B2, 1, 2, Wall, theta_c1, theta_c2);
//    std::vector<double> measures_g_b1 = get_two_intersect_measures(sphere_B1, sphere_G, 1, 0, Wall, theta_c1, theta_c2);
//    std::vector<double> measures_g_b2 = get_two_intersect_measures(sphere_G, sphere_B2, 0, 2, Wall, theta_c1, theta_c2);
//
//    std::vector<double> measures_three = get_three_intersect_measures();
//
//    V_b1_b2 = measures_b1_b2[0];    V_g_b1 = measures_g_b1[0];    V_g_b2 = measures_g_b2[0];
//    S_b1_b2 = measures_b1_b2[1];    S_g_b1 = measures_g_b1[1];    S_g_b2 = measures_g_b2[1];
//    V_g_b1_b2 = measures_three[0];
//    S_g_b1_b2 = measures_three[1];
//
//    net_measures[0] = V_b1 + V_b2 + V_g - V_b1_b2 - V_g_b1 - V_g_b2 + V_g_b1_b2;
//    net_measures[1] = S_b1 + S_b2 + S_g - S_b1_b2 - S_g_b1 - S_g_b2 + S_g_b1_b2;
//}

//std::vector<double> Composite_cluster::get_proj_SA()
//{
//
//    Circle circle_B1 (sphere_B1, Wall);
//    Circle circle_G  (sphere_G,  Wall);
//    Circle circle_B2 (sphere_B2, Wall);
//    std::vector<double> proj_SAs(3, 0.0);
//    double area_on_B1, area_on_B2, area_on_G;
//    double rb1, rb2, rg;
//    std::vector<double> cb1(3,0.0), cb2(3,0.0), cg(3,0.0);
//    double r_inter_b1g, r_inter_b2g, r_inter_b1b2;
//    rb1 = circle_B1.get_radius();   cb1 = circle_B1.get_centre();
//    rb2 = circle_B2.get_radius();   cb2 = circle_B2.get_centre();
//    rg  = circle_G.get_radius();    cg = circle_G.get_centre();
//    if(cg[0]!=0.0){printf("cg is not at origin"); abort();}
//
//    r_inter_b1g = circle_circle_int_radii(rb1, rg, cb1, cg);
//    r_inter_b2g = circle_circle_int_radii(rg, rb2, cg, cb2);
//    r_inter_b1b2 = circle_circle_int_radii(rb1, rb2, cb1, cb2);
//    double cdist_b1g = abs(0.0 - cb1[0]);
//    double cdist_b2g = abs(cb2[0] - 0.0);
//    double cdist_b1b2 = abs(cb2[0] - cb1[0]);
//
//    double d_b1g, d_b2g, d_b1b2;
//    double area_lens_b1g, area_lens_b2g, area_lens_b1b2;
//
//    if(Rb1 == 0.0 && Rb2 == 0.0)
//    {
//        area_lens_b1g = 0.0; area_lens_b2g = 0.0; area_lens_b1b2 = 0.0;
//        area_on_B1 = 0.0; area_on_B2 = 0.0; area_on_G = pi*rg*rg ;
//    }
//    else
//    {
//        if(!isnan(r_inter_b1g) && r_inter_b1g != 0.0 && !isnan(r_inter_b2g) && r_inter_b2g != 0.0)
//        {
//            d_b1g = x_intersection(rb1, rg, cdist_b1g);
//            d_b2g = x_intersection(rg, rb2, cdist_b2g);
//
//            area_lens_b1g = area_circular_segment (rb1, d_b1g) + area_circular_segment (rg, (cdist_b1g - d_b1g));
//            area_lens_b2g = area_circular_segment (rg, d_b2g) + area_circular_segment (rb2, (cdist_b2g - d_b2g));
//
//            if(!isnan(r_inter_b1b2) && r_inter_b1b2 != 0.0)
//            {
//                d_b1b2 = x_intersection(rb1, rb2, cdist_b1b2);
//                area_lens_b1b2 = area_circular_segment(rb1, d_b1b2) + area_circular_segment (rb2, (cdist_b1b2 - d_b1b2));;
//            }
//            else
//            {
//                area_lens_b1b2 = 0.0;
//            }
//
//            area_on_B1 = pi*rb1*rb1 - area_lens_b1g;
//            area_on_B2 = pi*rb2*rb2 - area_lens_b2g;
//            area_on_G = pi*rg*rg - area_lens_b1g - area_lens_b2g + area_lens_b1b2;
//        }
//        else{printf("B1 and G or B2 and G do not intersect"); abort();}
//
//    }
//
//    proj_SAs[0] = area_on_G;
//    proj_SAs[1] = area_on_B1;
//    proj_SAs[2] = area_on_B2;
//
//    //area_circular_segment (double R, double d)
//}
