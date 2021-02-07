//
//  TestHourglass.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/4/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include "TestHourglass.hpp"
#include "MC.hpp"
#include "Spherical_cap.hpp"
#include "Cylinder.hpp"
#include "HourglassCluster.hpp"

int main(int argc, const char * argv[])
{
    /* Here the wall surface is in the x-y plane and the hourglass extends in the z-direction */
    
    printf("Number of arguments = %d\n", argc);
    double z_wall = 0.0;
    /* Geometry related variables */
    double R1, R2, Rc;
    double L_cyl, L_thick;
    double theta_w, theta_f;
    double z_cyl_centre_1, z_cyl_centre_2;
    std::vector<std::vector<double> > box(3);
    box[0].resize(2);
    box[1].resize(2);
    box[2].resize(2);
    
    /* Monte-Carlo related variables */
    int n_points; int n_points_array_size; //This is for the case when checking convergence for n = [10^6, 2*10^6, ...]
    int aSeed[3];
    double delta;
    
    R1  = atof(argv[1]);
    R2  = atof(argv[2]);
    Rc = atof(argv[3]);
    L_cyl = atof(argv[4]);
    L_thick = atof(argv[5]);
    theta_w = atof(argv[6]);
    theta_f = atof(argv[7]);
    z_cyl_centre_1 = atof(argv[8]);
    z_cyl_centre_2 = atof(argv[9]);
    
    
//    n_points = atoi(argv[10]);       //Total points in millions in the box for MC
//    n_points = n_points*1e06;
    n_points_array_size = atoi(argv[10]);        //Size of
    aSeed[0] = atoi(argv[11]);       //Seed in x-direction for MC
    aSeed[1] = atoi(argv[12]);       //Seed in y-direction for MC
    aSeed[2] = atoi(argv[13]);       //Seed in z-direction for MC
    delta    = atof(argv[14]);
    
    /* These two files can be used to print out either points (surface or volume) or to get data i.e. volume and SA for different values of n_points using the for loop in n_points at the end */
//    FILE* dataFile = fopen(argv[15], "w");
    FILE* pointsfile = fopen(argv[15], "w");


    printf("delta= %10.4f\n", delta);
    
    double box_xy_dim = std::max({R1, Rc, R2});
    printf("initial box_xy_dim = %10.5f\n", box_xy_dim);
    
    if(theta_w <= pi/2.0 && theta_f > pi/2.0)
    {
        box_xy_dim = std::max({R1*sin(theta_w), Rc, R2});
    }
    else if(theta_w <= pi/2.0 && theta_f <= pi/2.0)
    {
        box_xy_dim = std::max({R1*sin(theta_w), Rc, R2*sin(theta_f)});
    }
    else if (theta_w > pi/2.0 && theta_f <= pi/2.0)
    {
        box_xy_dim = std::max({R1, Rc, R2*sin(theta_f)});
    }
    else if (theta_w > pi/2.0 && theta_f > pi/2.0)
    {
        box_xy_dim = std::max({R1, Rc, R2});
    }
    
    printf("final box_xy_dim = %10.5f\n", box_xy_dim);
    
    double buffer = 0.05*std::max({R1, Rc, R2});
    printf("buffer = %10.5f\n", buffer);
    
    box[0][0] = -(box_xy_dim + buffer);    box[1][0] = -(box_xy_dim + buffer);    box[2][0] = z_wall;
    box[0][1] = box_xy_dim + buffer;     box[1][1] = box_xy_dim + buffer;     box[2][1] = L_thick;
    
    printf("Box = \n");
    printf("Box[0][0] = %10.5f\t Box[0][1] = %10.5f\n", box[0][0], box[0][1]);
    printf("Box[1][0] = %10.5f\t Box[1][1] = %10.5f\n", box[1][0], box[1][1]);
    printf("Box[2][0] = %10.5f\t Box[2][1] = %10.5f\n", box[2][0], box[2][1]);
    
    /* Setting up the shapes */
    /* Wall cap */
    std::vector<double> centre_wall(3);   //Centre of the good patch. Only x-y
    centre_wall[0] = 0.0; centre_wall[1] = 0.0; centre_wall[2] = z_wall;
    
    double projected_rw = R1*sin(theta_w); //projected radii of the circles
    Spherical_cap WallCap (centre_wall, projected_rw, theta_w, z_wall);
    
    /* Free cap */
    std::vector<double> centre_free(3);   //Centre of the good patch. Only x-y
    centre_free[0] = 0.0; centre_free[1] = 0.0; centre_free[2] = L_thick;
    
    double projected_rf = R2*sin(theta_f); //projected radii of the circles
    Spherical_cap FreeCap (centre_free, projected_rf, theta_f, L_thick, 1);
    
    /* Cylinder */
    std::vector<std::vector<double> > circle_centres (2);
    circle_centres[0].resize(3);            circle_centres[1].resize(3);
    circle_centres[0][0] = 0.0;             circle_centres[1][0] = 0.0;
    circle_centres[0][1] = 0.0;             circle_centres[1][1] = 0.0;
    circle_centres[0][2] = z_cyl_centre_1;  circle_centres[1][2] = z_cyl_centre_2 ;
    
    Cylinder cylinder (Rc, circle_centres);
    
    Shape* WallCap_shape = &WallCap;
    Shape* Cylinder_shape = &cylinder;
    Shape* FreeCap_shape = &FreeCap;
    std::vector<Shape*> list_of_shape_pointers = {WallCap_shape, Cylinder_shape, FreeCap_shape};
    
    Hourglass Hglas (list_of_shape_pointers);
    printf("list size = %d\n", (int)list_of_shape_pointers.size());
    Shape* shape_ptr = &Hglas;
    
    for(int i=0; i<n_points_array_size; i++)
    {
        n_points = (i+1)*1e06;
        /* Setting the MC engine */
        MC mc_engine (n_points, box, aSeed);
        printf("i=%d\t box_volume = %10.5f\n", i, mc_engine.box_volume());
        std::vector<double> measures = mc_engine.calc_volume_SA(shape_ptr, delta);
        mc_engine.print_volume_points(pointsfile);
        double volume = measures[0];
        double SA = measures[1];
        printf("Hourglass volume = %10.15f\t Surface Area = %10.15f\n", volume, SA);
        //fprintf(dataFile, "%d\t%10.10f\t%10.10f\n", i, volume, SA);
    }
    
    return 0;
}
