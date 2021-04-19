//
//  Spherocylinder_cap_composite.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "Spherocylinder_cap_composite.hpp"


Spherocylinder_cap_composite::Spherocylinder_cap_composite()
{
    
}


Spherocylinder_cap_composite::~Spherocylinder_cap_composite()
{
    
}

Spherocylinder_cap_composite::Spherocylinder_cap_composite(std::vector<Spherical_cap*> caps, SpheroCylinder* spherocylinder, Surface* SurfacePtr):
list_of_spherical_caps(caps),
sphero_cylinder(spherocylinder),
surface_ptr(SurfacePtr)
{
    Shape* sphero_cylinder_ptr = (Shape*) spherocylinder;
    list_of_shapes.push_back(sphero_cylinder_ptr);
    
    for (size_t i=0; i<caps.size(); i++)
    {
        Shape* cap_i_ptr = (Shape*) caps[i] ;
        list_of_shapes.push_back(cap_i_ptr);
    }
    
}

void Spherocylinder_cap_composite::oppositeBoundCross()
{
    double single_patch_left_bound = surface_ptr->Patches[0].patch_boundaries()[0];
    double single_patch_right_bound = surface_ptr->Patches[0].patch_boundaries()[1];
}

std::vector<double> Spherocylinder_cap_composite::projected_SA()
{
    std::vector<double> result (3, 0.0);

    double left_circle_h_on_right, right_circle_h_on_left;
    std::vector<double> projected_areas;
    std::vector<Circle> circles(list_of_spherical_caps.size());
    
    double spherocylinder_proj_area = sphero_cylinder->AnalytProjSurfAreaWithPlaneIntersect();
    
    //This is currently assuming there is a single good patch surrounded by infinite bad patches.
    double single_patch_left_bound = surface_ptr->Patches[0].patch_boundaries()[0];
    double single_patch_right_bound = surface_ptr->Patches[0].patch_boundaries()[1];
    for(size_t i=0; i<circles.size(); i++)
    {
        circles[i] = list_of_spherical_caps[i]->get_circle();
        
    }
    
    //This is to check if left circle is bleeding onto the right patch and is beyond the right circle's boundary and vice versa
    double oppositeBoundDist_left_circle = (single_patch_right_bound - circles[0].get_centre()[0]);
    double h_badcap_oppositeBound_left_circle = circles[0].get_radius() - oppositeBoundDist_left_circle ;
    double oppositeBoundDist_right_circle = circles[1].get_centre()[0] - single_patch_left_bound ;
    double h_badcap_oppositeBound_right_circle = circles[1].get_radius() - oppositeBoundDist_right_circle ;
    
    
    double proj_area_left, proj_area_right;
    
    left_circle_h_on_right = single_patch_left_bound - circles[0].get_centre()[0];
    right_circle_h_on_left = circles[1].get_centre()[0] - single_patch_right_bound;
    double h_badcap_left_circle = circles[0].get_radius() - (-1*left_circle_h_on_right);
    double h_badcap_right_circle = circles[1].get_radius() - (-1*right_circle_h_on_left);
    
    bool right_on_left = (h_badcap_left_circle < h_badcap_oppositeBound_right_circle) ;
    bool left_on_right = (h_badcap_right_circle < h_badcap_oppositeBound_left_circle) ;
    
    double segment_area_left_circle = circles[0].segment_area(left_circle_h_on_right);
    double segment_area_right_circle = circles[1].segment_area(right_circle_h_on_left);
    if(isnan(segment_area_left_circle) && (left_circle_h_on_right >= circles[0].get_radius() ||  left_circle_h_on_right <= -circles[0].get_radius())) {segment_area_left_circle = 0.0;}
    
    if(isnan(segment_area_right_circle) && (right_circle_h_on_left >= circles[1].get_radius() || right_circle_h_on_left <= -circles[1].get_radius())) {segment_area_right_circle = 0.0;}
    
    
    if(h_badcap_oppositeBound_left_circle > 0 && h_badcap_oppositeBound_right_circle > 0) //This shows that each of the two caps have intersection with the opposite boundary of the good patch. ONLY possible when dB_sign is -1.
    {
        if(right_on_left &&  left_on_right) //This shows that the opposite caps intersection is more prominent (boundary with the patch i.e. area) on a side than its original cap.
        {
            proj_area_left = circles[1].segment_area(oppositeBoundDist_right_circle);
            
            proj_area_right = circles[0].segment_area(oppositeBoundDist_left_circle);
            
        }
        else if((right_on_left && !left_on_right) || (!right_on_left && left_on_right))
        {
            printf("Symmetry is broken for bad patch caps in Spherocylinder right_on_left=%d\t left_on_right=%d\n", right_on_left, left_on_right); abort();
        }
        else
        {
             proj_area_left = circles[0].area() - segment_area_left_circle;
             proj_area_right = circles[1].area() - segment_area_right_circle;
        }
    }
    else
    {
         proj_area_left = circles[0].area() - segment_area_left_circle;
         proj_area_right = circles[1].area() - segment_area_right_circle;
    }
    
    
    
    
//    printf("seg_area_left = %10.10f\t seg_area_right = %10.10f\n ", segment_area_left_circle, segment_area_right_circle);
    
    //        printf("[segment_area_i_on_ngbleft segment_area_i_on_ngbright] = [%10.5f %10.5f]\n", segment_area_i_on_ngbleft, segment_area_i_on_ngbright);
    
    
    result[0] = proj_area_left;
    result[1] = spherocylinder_proj_area;
    result[2] = proj_area_right ;
    return result;
}



int Spherocylinder_cap_composite::isInside(std::vector<double>& point)
{
    //Write a test to check this i.e. if the point is inside any of the caps then it is inside the composite cluster
    for (size_t i=0; i<list_of_shapes.size(); i++)
    {
        if (list_of_shapes[i]->isInside(point) == 1)
        {
            return 1;
        }
        else{continue;}
    }
    return 0;
}

int Spherocylinder_cap_composite::nearSurface(std::vector<double> & point, double delta)
{
    for (size_t i=0; i<list_of_shapes.size(); i++)
    {
        if (list_of_shapes[i]->nearSurface(point, delta) == 1)
        {
            int outOfOthers = 1;
            for (size_t j=0; j<list_of_shapes.size(); j++)
            {
                if(i != j)
                {
                    if (list_of_shapes[j]->isInside(point) == 1)
                    {outOfOthers = 0;}
                }
            }
            if(outOfOthers==1)
            {return 1;}
            else{return 0;} //continue;
        }
    }
    return 0;  //If it gets out of the for loop without returning 1 within it then returns 0
}

std::vector<double> Spherocylinder_cap_composite::xy_spread()
{
    std::vector<double> thisSpread (4,0.0); //[Xmin, Xmax, Ymin, Ymax]
    thisSpread = list_of_shapes[0]->xy_spread();
    
    for(size_t i = 1; i<list_of_shapes.size(); i++)
    {
        std::vector<double> spread_ith_cap = list_of_shapes[i]->xy_spread();
        if(spread_ith_cap[0] <= thisSpread[0]){thisSpread[0] = spread_ith_cap[0];}
        if(spread_ith_cap[1] >= thisSpread[1]){thisSpread[1] = spread_ith_cap[1];}
        if(spread_ith_cap[2] <= thisSpread[2]){thisSpread[2] = spread_ith_cap[2];}
        if(spread_ith_cap[3] >= thisSpread[3]){thisSpread[3] = spread_ith_cap[3];}
    }
    return thisSpread;
}

std::vector<double> Spherocylinder_cap_composite::threeDim_spread()
{
    //To check the breach of 5 box planes, excluding the wall surface
    //[Xmin, Xmax, Ymin, Ymax, Zmax]
    //Start with central caps spread as the spread of cluster and then modify based on other caps
    std::vector<double> thisSpread = list_of_shapes[0]->threeDim_spread();
//    printf("Inside Spherocylinder_cap_composite three dim spread \n");
//    printf("list_of_shapes size = %d\n", (int)list_of_shapes.size());
    for(size_t i = 1; i<list_of_shapes.size(); i++)
    {
        std::vector<double> iCap3DSpread = list_of_shapes[i]->threeDim_spread();
        if(iCap3DSpread[0] <= thisSpread[0]) { thisSpread[0] = iCap3DSpread[0]; }
        if(iCap3DSpread[1] >= thisSpread[1]) { thisSpread[1] = iCap3DSpread[1]; }
        if(iCap3DSpread[2] <= thisSpread[2]) { thisSpread[2] = iCap3DSpread[2]; }
        if(iCap3DSpread[3] >= thisSpread[2]) { thisSpread[3] = iCap3DSpread[3]; }
        if(iCap3DSpread[4] >= thisSpread[4]) { thisSpread[4] = iCap3DSpread[4]; }
    }
    return thisSpread;
}
