//
//  main_spherocylinder.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 3/17/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include "Spherocylinder_cap_composite.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "DynamicBox.hpp"
#include "MC.hpp"
#include "CheckBoundary.hpp"
#include "FreeEnergy.hpp"
#include "EvolveCluster.hpp"

using namespace std;


Patch setup_patch(double width_x, double width_y, double x_coordinate_centre, double theta, double z_surface)
{
    std::vector<double> centre(3);   //Centre of the good patch. Only x-y
    centre[0] = x_coordinate_centre;
    centre[1] = 0.0;
    centre[2] = z_surface;
    std::vector<double> dim(2,0.0);

    dim[0] = width_x; dim[1] = width_y; //Common dimensions for all Good patches
    
    Patch patch (theta, centre, dim);
    return patch ;
}

//void Calc_Volume_SA(Shape* shape, MC* mc_engine, FILE* output)
//{
//
//}


/* Have this function setup stripes based on num_patches and width. If the width of a certain type of patch is -1, then that means that patch is infinite*/
//Stripes setup_stripes (int num_patches, double pG_width_x, double pG_width_y, double pB_width_x, double pB_width_y, double theta_cg, double theta_cb, double x_coordinate_centre, double z_surface)
//{
//
//}


int main(int argc, const char * argv[]) {
    
    //Dynamically add pathces as the bad patch cap grows.
    
    //Density
    double Rho ;
    
    //Variables for spherocylinder
    double d_cyl_length, d_chord_length;
    double cyl_max; //Maximum length of the cylinder
    
    //Variables for caps
    double d_Rg, d_Rb;
    double Rg_max, Rb_max;
    
    //Variables for Patches
    double z_surface = 0.0; //Unless stated otherwise
    
    double pG_width_x, pG_width_y, pB_width_x, pB_width_y, theta_cg, theta_cb;
    int num_patches;       //This is the maximum number of patches that are allowed.
    
    // Monte-Carlo related variables
    double point_density; //Constant point density for all additional points added dynamically
    int n_points;       //Initial number of points in the first good patch box calcualted from the density
    int aSeed[3];
    double delta;
    
    //Dynamic box related variables
    double extension_length; //The fixed extension length in each direction
    std::string tag;
     /* Reading all variables from command line */
    pG_width_x  = atof(argv[1]);     //x width of the good patches (same for all good patches) (Angstroms)
    pG_width_y  = atof(argv[2]);     //y width of the good patches (same for all good patches) (Angstroms)
    pB_width_x  = atof(argv[3]);     //x width of the good patches (same for all good patches) (Angstroms)
    pB_width_y  = atof(argv[4]);     //y width of the good patches (same for all good patches) (Angstroms)
    theta_cg    = atof(argv[5]);       //Good patch contact angle
    theta_cb    = atof(argv[6]);       //Bad patch contact angle
    d_Rg        = atof(argv[7]);       //increments in good patch cap radius (Angstroms i.e. d_Rg=0.05 A)
    d_Rb        = atof(argv[8]);       //increments in bad patch cap radius  (Angstroms)
    d_cyl_length= atof(argv[9]);
    d_chord_length= atof(argv[10]);
    Rg_max      = atof(argv[11]);       //maximum limit of good patch cap.    (Angstroms)
    Rb_max      = atof(argv[12]);
    cyl_max     = atof(argv[13]);        //maximum limit of bad patch cap      (Angstroms)
    point_density = atof(argv[14]);     // Constant point density for MC. (points/Angstrom^3)
    aSeed[0]    = atoi(argv[15]);       //Seed in x-direction for MC
    aSeed[1]    = atoi(argv[16]);       //Seed in y-direction for MC
    aSeed[2]    = atoi(argv[17]);       //Seed in z-direction for MC
    delta       = atof(argv[18]);      //buffer region width=(2*\delta) for surface points (Angstroms)
    Rho         = atof(argv[19]);   //Moles per m3
    num_patches = atoi(argv[20]);
    extension_length = atof(argv[21]);
    tag.assign(argv[22]);
    
    std::string output_points_file_name = tag + "_points.txt";
    std::string outputfile_name = tag + "_Volume_SA_data.txt" ;
    FILE* output_points_file = fopen(output_points_file_name.c_str(), "w");
    FILE* outputfile = fopen(outputfile_name.c_str(), "w");
    //Starting with a single good patch and two infinite bad patches beside it.
    
    Patch good_patch = setup_patch(pG_width_x, pG_width_y, 0.0, theta_cg, z_surface);
    
    
    std::vector<double> centre_good ({0.0, 0.0, z_surface});
    std::vector<Patch> list_of_patches {good_patch};
    std::vector<std::vector<double> > orientations {centre_good};
    
    Surface* surface_ptr;
    Stripes stripes (list_of_patches, orientations, z_surface);
    surface_ptr = &stripes;
    
    //Here starting with a much smaller box to check if dynamic box works.
    double box_height = (Rg_max >= Rb_max) ? (double)(Rg_max/6.0) : (double)(Rb_max/6.0) ;
    // ? 2.0*Rg_max : 2.0*Rb_max ;
    
    stripes.calc_box(box_height);
    double initial_box_volume = stripes.box_volume();
    n_points = (int) (initial_box_volume * point_density);
    printf("n_points = %d\n", n_points);
    for(size_t i =0; i<stripes.box.size(); i++)
    {
        for(size_t j =0; j<stripes.box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
        }
    }
    
    std::vector<double> mc_volume_SA(2,0.0);
    //Instantiating dynamic_box
    DynamicBox dynamic_box (stripes.box, extension_length);
    MC mc_engine (n_points, stripes.box, aSeed); //branch_comm);
    CheckBoundary check_boundary (surface_ptr, &mc_engine , &dynamic_box, point_density);
    Shape* Cluster_shape_ptr;
    
    std::vector<double> maximum_limits ({Rg_max, Rb_max , cyl_max});
    std::vector<double> increments({d_Rg, d_Rb , d_cyl_length});
    std::vector<double> patch_widths({pG_width_x, pG_width_y, pB_width_x, pB_width_y});
    
    EvolveCluster evolve_cluster (&check_boundary,  &mc_engine, mc_volume_SA, maximum_limits , increments, theta_cg,  theta_cb, patch_widths , z_surface, delta, Rho);
    
    double Rg, Rb, cyl_length;
    double N, Volume;
    double SA;
    double projected_SA;
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> box_breach; //array of (0,1) to check box surface breach
    
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1;
    int len_cyl = (int) ((cyl_max-0.0)/d_cyl_length) + 1;
    
    for(int i = 0 ; i< len_Rg; i++) // 1
    {
        Rg = 0.0 + i*d_Rg ;
        printf("Rg = %10.10f\n", Rg);
        cyl_length = 0.0; //  startingRg + i*d_Rg            //sphere's radius
        double projected_rg = Rg*sin(theta_cg);
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        Cluster_shape_ptr= &GoodCap;
        check_boundary.ManageBoxBreach(Cluster_shape_ptr);

        
        std::vector<double> shape_xy_spread = Cluster_shape_ptr->xy_spread();
        std::vector<double> cPatchBounds = good_patch.patch_boundaries();
//        printf("shape_xy_spread = %10.10f\t%10.10f\t%10.10f\t%10.10f\n", shape_xy_spread[0], shape_xy_spread[1], shape_xy_spread[2], shape_xy_spread[3]);
        if(cPatchBounds[0] < shape_xy_spread[0] && cPatchBounds[1] > shape_xy_spread[1]) //Only spherical cap within the boundary
        {
            if((cPatchBounds[0] > shape_xy_spread[0] - d_Rg*sin(theta_cg)) && (cPatchBounds[0] <= shape_xy_spread[0]) && (cPatchBounds[1] < shape_xy_spread[1] + d_Rg*sin(theta_cg)) && (cPatchBounds[1] >= shape_xy_spread[1])  )
            {
                printf("Inside spherocylinder patch \n");
                for(int j = 0 ; j< len_cyl; j++) //j=0 is just the spherical cap that touches the boundary
                {
                    cyl_length = 0.0 + j*d_cyl_length;
                    printf("cyl_length = %10.10f\n", cyl_length);
                    std::vector<double> normal_to_the_circle({0.0 , 1.0, 0.0});
                    
                    SpheroCylinder sphero_cylinder (GoodCap, cyl_length, normal_to_the_circle, z_surface);
                    
                    evolve_cluster.EvolveBadCapWithSpherocylinder (&sphero_cylinder, Cluster_shape_ptr, cyl_length, d_chord_length, Rg, outputfile, output_points_file);
                }
            }
            else
            {
                //Only spherical cap
                mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                Volume = mc_volume_SA[0]; //GoodCap.getVolume();
                SA = mc_volume_SA[1];
                projected_SA = GoodCap.projected_SA();
                N = (Volume * 1e-30) * Rho * avogadro ;
                fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, N, Volume, SA, projected_SA);
            }
        }
        
        else if (cPatchBounds[0] > shape_xy_spread[0] && cPatchBounds[1] < shape_xy_spread[1])
        {
            //Not allowing rg to grow beyond the good patch boundary
            break;
        }
    }
    
}

//Earlier directly spherocylinder now for each cyl_length, bad cap
//
//
//Cluster_shape_ptr= &sphero_cylinder;
//
//check_boundary.ManageBoxBreach(Cluster_shape_ptr);
//
//mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
//Volume = mc_volume_SA[0]; //GoodCap.getVolume();
//SA = mc_volume_SA[1];
//projected_SA = sphero_cylinder.AnalytProjSurfAreaWithPlaneIntersect();
//N = (Volume * 1e-30) * Rho * avogadro ;
//
//if(j == len_cyl-2)
//{
//    printf("inside printig surface points\n");
//    mc_engine.print_surf_points(output_points_file);
//
//}
//fprintf(outputfile, "%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\t%10.10f\n" , Rg, cyl_length, N, Volume, SA, projected_SA);
