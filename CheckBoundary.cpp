//
//  CheckBoundary.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/18/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "CheckBoundary.hpp"

CheckBoundary::CheckBoundary()
{
    
}

CheckBoundary::~CheckBoundary()
{
    
}

CheckBoundary::CheckBoundary(Surface* surf, MC* mc , DynamicBox* dyn_box, double density):
surface(surf), mc_engine(mc), dynamic_box(dyn_box), point_density(density)
{
    
}


void CheckBoundary::ManageBoxBreach(Shape* shape)
{
    cluster = shape;
    std::vector<int> box_breach = dynamic_box->CheckBoxBreach(cluster);
//    for(int i=0 ; i<(int)box_breach.size(); i++)
//    {
//        printf("box_breach [%d] = %d\n", i, box_breach[i]);
//    }
    
    std::vector<int> find_box_breach = FindBoxBreach(box_breach);
    if(!find_box_breach.empty()) //find_box_breach being empty suggests that box is not breached.
    {
        ExpandBox (find_box_breach);
        ManageBoxBreach(cluster);
        //dynamic_box->print_box();
    }
    
}


std::vector<int> CheckBoundary::FindBoxBreach(std::vector<int> box_breach)
{
    std::vector<int> result; //[breach in x-direction, y-direction, z-direction] <== Boolean (0,1)
    for(int i=0; i<(int)box_breach.size(); i++)
    {
        if(box_breach[i] == 1)
        {
            if(i==0)
            {   //x-direction breach
                if(box_breach[i+1] == 1)
                {
                    result.push_back(0); //Here 0 is the breach direction of x.
                }
                else{printf("box not breached symmetrically in x-direction\n"); abort();}
                
            }
            if(i==2)
            { //y-direction breach
                if(box_breach[i+1] == 1)
                {
                    result.push_back(1);
                }
                else{printf("box not breached symmetrically in y-direction\n"); abort();}
                
            }
            if(i==4)
            {//z-direction breach
                result.push_back(2);
                    
            }
            
        }
    }
    return result;
}

void CheckBoundary::ExpandBox (std::vector<int> find_box_breach)
{
    int direction_to_expand_in;
    std::vector<std::vector<double>> new_box;
    int extra_points_per_added_box;
//    printf("find_box_breach size = %d\n", (int)find_box_breach.size());
    for(size_t i=0 ; i<find_box_breach.size(); i++)
    {
//        printf("find_box_breach [%d] = %d\n, ", (int)i, find_box_breach[i]);
        direction_to_expand_in = find_box_breach[i];
        extra_points_per_added_box = dynamic_box->add_fix_box(direction_to_expand_in, point_density); //This also updates the box within the dynamic_box class
        new_box = dynamic_box->get_box();
        mc_engine->add_points(new_box, direction_to_expand_in, extra_points_per_added_box); //This also updates the box within the MC class
    }
}

bool CheckBoundary::CheckSpherocylinderBadPatch(double cyl_length, int dB_sign, double patch_width, std::vector<double> centre_left, std::vector<double> centre_right, double projected_radius)
{
    /*
     Currently this for an infinite bad patches surrounding a good patch. Therefore for dB_sign=1 i.e. the centre is on the bad patch,
     this function directly returns True, since there is no boundary to dictate the "breaking" of bad patch cap loop.
     */
    if(dB_sign == 1)
    {
        return true;
    }
    else
    {
        std::vector<double> caps_boundary_near_opposite_patch_boundary(2, -1.0); //[left, right]
        caps_boundary_near_opposite_patch_boundary[0] = projected_radius + centre_left[0] ;
        caps_boundary_near_opposite_patch_boundary[1] = centre_right[0] - projected_radius ;
        double dB_prime_left, dB_prime_right; //distance from opposite patch boundary;
        double chord_prime, chord_prime_left, chord_prime_right;
       
        bool left_patch_within_conditon = (caps_boundary_near_opposite_patch_boundary[0] <= 0.5 * patch_width);
        bool right_patch_within_conditoin = (caps_boundary_near_opposite_patch_boundary[1] >= -0.5 * patch_width);
        
        if(left_patch_within_conditon && right_patch_within_conditoin)
        {return true;}
        else if ((right_patch_within_conditoin && !left_patch_within_conditon) || (!right_patch_within_conditoin && left_patch_within_conditon) )
        {
            printf("Bad patches are not symmetrically crossing good patch bounds in CheckBoundary.cpp!\n");
            abort();
        }
        else
        {
            dB_prime_left = 0.5 * patch_width - centre_left[0];
            dB_prime_right = centre_right[0] - (-0.5 * patch_width);
            
            double inside_sq_left = projected_radius * projected_radius - dB_prime_left * dB_prime_left;
            double inside_sq_right = projected_radius * projected_radius - dB_prime_right * dB_prime_right;
            
            chord_prime_left = sqrt(inside_sq_left);
            chord_prime_right = sqrt(inside_sq_right);
            
            chord_prime = (chord_prime_left == chord_prime_right) || (abs(chord_prime_left - chord_prime_right) <= 1e-5) ? chord_prime_left : -100 ;
            
            if(chord_prime == -100)
            {
                printf("Chords on opposite boundary are unequal by more than 1e-5 diff = %10.10f\n", abs(chord_prime_left - chord_prime_right));
                abort();
                
            }
            else if(2 * chord_prime > cyl_length)
            {
                return false;
            }
            else {return true;}
        }
    }
    
}


