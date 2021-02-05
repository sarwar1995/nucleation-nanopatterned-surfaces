//
//  main.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 7/1/20.
//  Copyright © 2020 Sarwar Hussain. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string.h>
//#include "Composite_cluster.hpp"
//#include "Surface/Surface.hpp"
#include "Patch.hpp"
#include "Stripes.hpp"
#include "Spherical_cap.hpp"
#include "Composite_cluster.hpp"
#include "MC.hpp"
#include <chrono>

//Define these but CONSIDER taking input
double Rho;         //Moles per m3
double Mu;
double Sigma;       //Liquid-crystal surface tension
double T;
#define avogadro 6.022140857e23;
#define kb 1.38064852e-23;


using namespace std;
double free_energy (double volume, double SA, std::vector<double>& projected_SA, double theta_good, double theta_bad)
{
    /* Here projected_SA has the projected surface area lying exclusively on the [0]: good, [1]:bad on left and [2]:bad on right patches */
    double G;
    double Vcomp = -1*volume* Rho * Mu;
    double Scomp = Sigma*SA;
    double SA_proj_good = projected_SA[0];
    double SA_proj_bad_1 = projected_SA[1];
    double SA_proj_bad_2 = projected_SA[2];
    double Sproj_comp = -Sigma*(cos(theta_good)*SA_proj_good + cos(theta_bad)*SA_proj_bad_1 + cos(theta_bad)*SA_proj_bad_2);
    G = Vcomp + Scomp + Sproj_comp;
    return G;
    
}

//double calcCompCluster (Spherical_cap& GoodCap, Stripes& stripes, MC& mc_engine, std::vector<double> c_bad_left, std::vector<double> c_bad_right, double projected_rb, double theta_cb, double z_surface)
//{
//
//}


//Needs two dynamically allocated arrays (vectors?) of indices or position of points that lie inside the cluster and near the surface.

int main(int argc, const char * argv[]) {
    
    /// Inputs from command line ///
  
    /* Surface related variables */
    double z_surface = 0.0; //Unless stated otherwise
    double p_width_x, p_width_y, theta_cg, theta_cb;
    
    /* Cluster related variables */
    double d_Rg, d_Rb;
    double Rg_max, Rb_max;
//    int Nmin, Nmax;
//    double dN;
    
    /* Monte-Carlo related variables */
    int n_points;
    int aSeed[3];
    double delta;
    
    /* CellList related variables */
    std::vector<double> R_cutoff(3, 0.0);
    
    /* I/O related variables */
//    FILE* outputfile;
    FILE* outputfile_spreadout;
    FILE* outputfile_not_spreadout;
    FILE* pointsFile;
    
    /* Reading all variables from command line */
    p_width_x  = atof(argv[1]);     //x width of the patches (same for all)
    p_width_y  = atof(argv[2]);     //y width of the patches (same for all)
    theta_cg = atof(argv[3]);       //Good patch contact angle
    theta_cb = atof(argv[4]);       //Bad patch contact angle
    d_Rg     = atof(argv[5]);       //increments in good patch cap radius
    d_Rb     = atof(argv[6]);       //increments in bad patch cap radius
    Rg_max   = atof(argv[7]);       //maximum limit of good patch cap
    Rb_max   = atof(argv[8]);       //maximum limit of bad patch cap
    n_points = atoi(argv[9]);       //Total points in millions in the box for MC
    n_points = n_points*1e06;
    aSeed[0] = atoi(argv[10]);       //Seed in x-direction for MC
    aSeed[1] = atoi(argv[11]);       //Seed in y-direction for MC
    aSeed[2] = atoi(argv[12]);       //Seed in z-direction for MC
    delta    = atof(argv[13]);      //buffer region width=(2*\delta) for surface points
    R_cutoff[0] = atof(argv[14]);   //Cell size in x for CellList
    R_cutoff[1] = atof(argv[15]);   //Cell size in y for CellList
    R_cutoff[2] = atof(argv[16]);   //Cell size in z for CellList
    //outputfile = fopen (argv[17],"w"); // Output file
    //pointsFile = fopen (argv[12],"w");
    std::string tag(argv[17]);
    double startingRg = atof(argv[18]);
    
    std::string outFileName_spreadout = tag + "-spreadout.txt";
    std::string outFileName_not_spreadout = tag + "-not_spreadout.txt";
    std::string outFileName_nearsurf_points = tag + "near_surf_points.txt";
    std::cout<<"outFileName_spreadout = "<<outFileName_spreadout<<std::endl;
    std::cout<<"outFileName_not_spreadout = "<<outFileName_not_spreadout<<std::endl;
    outputfile_spreadout = fopen (outFileName_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the bad patch
    outputfile_not_spreadout = fopen (outFileName_not_spreadout.c_str(),"w"); // Output file for the centres of bad patches that are on the good patch
    pointsFile = fopen(outFileName_nearsurf_points.c_str(), "w");
    
    printf("d_Rg, d_Rb = %10.5f %10.5f\n",d_Rg, d_Rb);
    printf("Rg_max, Rb_max = %10.5f %10.5f\n",Rg_max, Rb_max);
    printf("delta=%10.5f\n",delta);
    int len_Rg = (int) ((Rg_max-0.0)/d_Rg) + 1; //This 1 is added because for loop is i<len_Rg
    int len_Rb; //Calculated seperately for each Rg. Rb's are symmetric for the 2 patches
    
    /* Setting the good patch */
    std::vector<double> centre_good(3);   //Centre of the good patch. Only x-y
    centre_good[0] = 0.0; centre_good[1] = 0.0; centre_good[2] = z_surface;
    std::vector<double> dim(2,0.0);
    dim[0] = p_width_x; dim[1] = p_width_y; //Common dimension for all patches
    Patch good_patch (theta_cg, centre_good, dim);
    
    /* Setting the bad patches */
    std::vector<double> centre_bad_left(3); //(-a, 0, 0)
    std::vector<double> centre_bad_right(3);//(a, 0, 0)
    centre_bad_left[0] = -p_width_x; centre_bad_left[1] = 0.0; centre_bad_left[2] = z_surface;
    centre_bad_right[0] = p_width_x; centre_bad_right[1] = 0.0; centre_bad_right[2] = z_surface;
    Patch bad_patch_left (theta_cb, centre_bad_left, dim);
    Patch bad_patch_right (theta_cb, centre_bad_right, dim);
    
    /* Setting the surface */
    std::vector<Patch> list_of_patches;
    std::vector<std::vector<double> > orientations;
    orientations.push_back(centre_good);        //Good Patch centre
    orientations.push_back(centre_bad_left);    //Bad left Patch centre
    orientations.push_back(centre_bad_right);   //Bad right Patch centre
    list_of_patches.push_back(good_patch);
    list_of_patches.push_back(bad_patch_left);
    list_of_patches.push_back(bad_patch_right);
    
    Stripes stripes (list_of_patches, orientations, z_surface);
    
    /* Setting the height of the box to R_max */
    double box_height = (Rg_max >= Rb_max) ? 2.0*Rg_max : 2.0*Rb_max ;
    //(theta_cg < pi/2.0) ? Rg_max : 2.0*Rg_max ;
    
    /* Calculating box. Here box[i] = 2-dimensional and represents the boundaries of the box in each direction. */
    stripes.calc_box(box_height);
    for(size_t i =0; i<stripes.box.size(); i++)
    {
        for(size_t j =0; j<stripes.box[i].size(); j++)
        {
            printf("box[%d][%d] = %10.5f\t",(int)i,(int)j,stripes.box[i][j]);
        }
    }
    
    /* Setting the MC engine */
    MC mc_engine (n_points, stripes.box, aSeed);
    
    
    printf("box_volume = %10.5f\n", mc_engine.box_volume());
    printf("potency good=%10.10f\t potency bad=%10.10f\n",potency_factor(theta_cg), potency_factor(theta_cb));
    
    /* Generating the cellList */
    CellList cell_list (mc_engine.get_points(), R_cutoff, mc_engine.get_box());
    
    /* Setting up vectors for N_arrays */
//    std::vector<double> NArray_Gmin;
//    std::vector<double> NArray_confs;
//    std::vector<std::vector<double> > NArray_quant;
    
    /* Setting up required variables and output data structures*/
    double Volume;
    double SA;
    double Volume_MC;
    double SA_MC;
    double Rg, Rb, Rb_min;
    double d_B;
    std::vector<double> c_bad_left(3, 0.0), c_bad_right(3, 0.0); //Projected centres of bad spherical caps. x-values will depend on Rb values.
    c_bad_left[1] = 0.0; c_bad_right[1] = 0.0;
    c_bad_left[2] = z_surface; c_bad_right[2] = z_surface;
    
    
    std::vector<int> stripes_bounds; //array of (0,1) to check crossing of boundaries
    std::vector<int> stripes_box_breach; //array of (0,1) to check box surface breach
    
    std::vector<double> mc_volume_SA(2,0.0);
    
    printf("Before starting Rg loop\n");
    
    Shape* Cluster_shape_ptr;
    
    /* Starting the loop for Rg.*/
    for(int i = 0 ; i<1; i++) //len_Rg
    {
        printf("i=%d\t",i);
        Rg = startingRg + i*d_Rg ;  //  0.0 + i*d_Rg            //sphere's radius
        double projected_rg = Rg*sin(theta_cg); //projected radii of the circles
        printf("projected_rg = %10.10f\n",projected_rg);
        
        Spherical_cap GoodCap (centre_good, projected_rg, theta_cg, z_surface);
        
        Cluster_shape_ptr= &GoodCap;
        
        stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
        
        if(stripes_bounds[0] == 0) //Cap on central cluster is within bounds
        {
           auto t_begin = std::chrono::high_resolution_clock::now();
            Volume = GoodCap.getVolume();
            SA = GoodCap.getSA();
            auto t_end = std::chrono::high_resolution_clock::now();
               auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_begin);
            int zero = 0;
        fprintf(outputfile_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume, SA, duration.count(), stripes_bounds[0], zero);
            fprintf(outputfile_not_spreadout, "%10.10f\t%d\t%10.10f\t%10.10f\t%lld\t%d\t%d\n", Rg, zero, Volume, SA, duration.count(), stripes_bounds[0], zero);
        }
        else
        {
            
            /* Rb >= sqrt(Rg^2 - (a/2)^2 ) */
            Rb_min = sqrt(Rg*Rg - ((p_width_x*p_width_x)/4.0));
            len_Rb = (int) ((Rb_max-Rb_min)/d_Rb) + 1;
            
            //Two cases for each Rb value. +1: spread out, -1: not spread out
            std::vector<int> dB_list = {+1, -1};
            for(int k=0; k<dB_list.size(); k++)
            {
                for(int j = 0 ; j<len_Rb; j++)
                {
                    printf("j = %d\tdB_list = %d\n",j, dB_list[k]);
                    /* Have a file output the values of near surface points*/
                    Rb = Rb_min + j*d_Rb ;
                    
                    if(Rb > Rb_max){throw std::invalid_argument("Invalid Rb > Rb_max");}
                    double projected_rb = Rb*sin(theta_cb);
                    
                    
                    d_B = dB_list[k] * sqrt(Rb*Rb - Rg*Rg + ((p_width_x*p_width_x)/4.0)) ;
                    /* Symmetric caps on either side of central patch */
                    c_bad_left[0] = -(p_width_x/2.0) - d_B;
                    c_bad_right[0] = (p_width_x/2.0) + d_B;
                    
                    
                    /*Checking that centres satisfy cb1 < 0 and cb2 > 0 and */
                    if(c_bad_left[0] >= 0.0 || c_bad_right[0] <= 0.0)
                    {
                        printf("cb1 = %10.10f\tcb2=%10.10f\n",c_bad_left[0],c_bad_right[0]);
                        break;
                    }
                    
                    /*Check for  rb1 + cb1 < rg < rb1 - cb1 which should be satisfied by the previous condition */
                    if(!(Rg < Rb - c_bad_left[0] && Rg > Rb + c_bad_left[0]))
                    {
                        printf("j=%d\t Rb=%10.10f Rg=%10.10f dB=%10.10f dB_list[k]=%d\n",j, Rb, Rg, d_B, dB_list[k] );
                        throw std::invalid_argument("Invalid Cb1 and Cb2");
                    }
                    
                    
                    
                    Spherical_cap BadCap_left(c_bad_left, projected_rb, theta_cb, z_surface);
                    Spherical_cap BadCap_right(c_bad_right, projected_rb, theta_cb, z_surface);
                    std::vector<Spherical_cap> capsList = {GoodCap, BadCap_left, BadCap_right};
                    
                    Composite_cluster Comp_cluster (capsList, stripes);
                    Cluster_shape_ptr = &Comp_cluster;
                    stripes_bounds = stripes.monitor_cluster_spread(Cluster_shape_ptr);
                    stripes_box_breach = stripes.monitor_box_breach(Cluster_shape_ptr);
                    
                    for(size_t q=0; q<stripes_box_breach.size(); q++)
                    {
                        if(stripes_box_breach[q] == 1)
                        {
                            printf("Box was breached for Rg=%10.10f\tRb=%10.10f\t and db=%d\n",Rg, Rb,dB_list[k]);
                            abort();
                        }
                    }
                    
                    if(stripes_bounds[1]==1 && stripes_bounds[2]==1)
                    {
                        printf("Bad patch bounds are crossed at Rb=%10.10f\t and db=%d\n",Rb,dB_list[k]);
                        break;
                    }
                    
                    
                    auto t1 = std::chrono::high_resolution_clock::now();
                    try
                    {
                        
                        if(j==0) //Beginning of growing bad patches
                        {
                           
                            mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                            //This printing part should be after the calc_volume_SA call
                            if(dB_list[k] == 1)
                            {
                                fprintf(pointsFile,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list[k]);
                                mc_engine.print_surf_points(pointsFile);
                            }
                            
                        }
                        else
                        {
                            mc_volume_SA = mc_engine.calc_volume_SA(Cluster_shape_ptr, delta);
                            
                            //mc_volume_SA = mc_engine.update_volume_SA(Cluster_shape_ptr, &cell_list, delta);
                        }
                    }
                    catch(const std::logic_error& e)
                    {
                        if(Rg==0.0 && Rb==0.0)
                        {
                            Volume_MC=0.0;
                            SA_MC = 0.0;
                        }
                        else
                        {
                            std::cerr << e.what();
                            abort();
                        }
                    }
                    auto t2 = std::chrono::high_resolution_clock::now();
                    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
                    Volume_MC = mc_volume_SA[0];
                    SA_MC = mc_volume_SA[1];
                    if(dB_list[k] == 1)
                    {
                        fprintf(outputfile_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC, SA_MC, duration_ms.count(),stripes_bounds[0], d_B);}
                    else
                    {
                       fprintf(outputfile_not_spreadout,"%10.10f\t%10.10f\t%10.10f\t%10.10f\t%lld\t%d\t%10.10f\n", Rg, Rb, Volume_MC, SA_MC, duration_ms.count(),stripes_bounds[0], d_B);
                    }
                    //fprintf(outputfile, "%10.10f\t%10.10f\t%lld\t%d\n", Volume_MC, SA_MC, duration_ms.count(), stripes_bounds[0]);
                    }
            }
        }
        
    }
    //fclose(outputfile);
    fclose(outputfile_spreadout);
    fclose(outputfile_not_spreadout);
    return 0;
}




//    Rg_max   = atof(argv[5]);
//    Rb_max   = atof(argv[6]);
//    Nmin     = atof(argv[7]);
//    dN       = atof(argv[8]);
//    Nmax     = atof(argv[9]);
//
//    int len_Rg = (int) (Rg_max/d_Rg) + 1;
//    int len_N  = (int) ((Nmax-Nmin)/dN) + 1;
//

//    //Free energy calculations.
//    double G_composite, N_composite;
//    vector<double>               NArray_Gmin(len_N, 0.0);       //Gmin for each N
//    vector<double>               NArray_confs(len_N, 0.0);      // number of confs for each N
//    vector<vector<double> >      NArray_quant(len_N);           //[N, Gmin, Rg_min, Rb_min]
//
//    //Starting Rg loop. For good branch
//    for(int i = 0 ; i<len_Rg; i++)
//    {
//        double Rg = 0.0 + i*d_Rg ;              //sphere's radius
//        double Rb;
//        double rg = Rg*sin(theta_cg), rb;       //projected radii of the circles
//
//        if(rg < 3*(p_width/2.0))                //rg is contained within the two bad patches
//        {
//            std::vector<double> Cg (3,0.0);         //Centres of spheres. g=good, b1=bad 1 on left, b2=bad 2 on right
//            std::vector<double> Cb1 (3,0.0);
//            std::vector<double> Cb2 (3,0.0);
//            std::vector<double> origin (3,0.0);
//            std::vector<double> z_axis {0.0, 0.0, 1.0};
//
//            Cg[0] = 0.0 ; Cg[1] = 0.0 ; Cg[2] = 0.0 - Rg*cos(theta_cg) ;
//
//            Sphere sphere_G (Rg, Cg);
//            Plane  wall     (origin, z_axis);
//            std::vector<Sphere> Spheres(3);
//            Spheres[0] = sphere_G;
//
//            if(rg <= (p_width/2.0))                     //good does not intersesct with either bad patches
//            {
//                Rb = 0.0;
//                rb = 0.0;
//                Sphere sphere_B1 (Rb, Cb1);
//                Sphere sphere_B2 (Rb, Cb2);
//                Spheres[1] = sphere_B1;
//                Spheres[2] = sphere_B2;
//
//                Composite_cluster Cluster(Spheres, wall);
//                Cluster.calc_intersections();
//
//                std::vector<double> net_measures = Cluster.get_net_measures();
//                std::vector<double> proj_SAs     = Cluster.get_proj_SA();
//                double Vnet                = net_measures[0];
//                double Snet                = net_measures[1];
//                double V_comp              = -1.0 * Rho * Mu * Vnet;
//                double S_comp              = Sigma*(Snet - cos(theta_cg)*(proj_SAs[0]) - cos(theta_cb)*(proj_SAs[1] + proj_SAs[2]));
//                G_composite         = V_comp + S_comp;
//                N_composite         = Vnet * Rho * avogadro;
//                add_to_N (N_composite, G_composite, Rg, Rb, dN, Nmin, len_N, NArray_Gmin, NArray_confs, NArray_quant);
//            }
//
//            else                                    //Rg intersescts with both bad patches
//            {
//                double rb_min = sqrt(rg*rg - (p_width/2.0)*(p_width/2.0)) ; //min value of the projected radius
//                double Rb_min = rb_min/sin(theta_cb);
//                int len_Rb    = (int) ((Rb_max - Rb_min)/d_Rb);
//                for(int j=0; j<len_Rb; j++)
//                {
//                    Rb = Rb_min + j*d_Rb ;  //two possibilities for each centre, but Symmetric
//                    rb = Rb*sin(theta_cb);
//                    double db = sqrt(rb*rb - rg*rg + (p_width/2.0)*(p_width/2.0));
//
//                    //CASE 1: Both centres lie on the bad patches
//                    Cb1[0] = -(p_width/2.0) - db;   Cb1[1] = 0.0;   Cb1[2] = 0.0 - Rb*cos(theta_cb) ;
//                    Cb2[0] = (p_width/2.0) + db;    Cb1[1] = 0.0;   Cb2[2] = 0.0 - Rb*cos(theta_cb) ;
//                    if(Cb1[0] - rb > -(3*(p_width/2.0)) && Cb2[0] + rb < (3*(p_width/2.0)))
//                    {
//                        Sphere sphere_B1 (Rb, Cb1);
//                        Sphere sphere_B2 (Rb, Cb2);
//                        Spheres[1] = sphere_B1;
//                        Spheres[2] = sphere_B2;
//
//                        Composite_cluster Cluster(Spheres, wall);
//                        Cluster.calc_intersections();
//
//                        std::vector<double> net_measures = Cluster.get_net_measures();
//                        std::vector<double> proj_SAs     = Cluster.get_proj_SA();
//                        double Vnet        = net_measures[0];
//                        double Snet        = net_measures[1];
//                        double V_comp      = -1.0 * Rho * Mu * Vnet;
//                        double S_comp      = Sigma*(Snet - cos(theta_cg)*(proj_SAs[0]) - cos(theta_cb)*(proj_SAs[1] + proj_SAs[2]));
//                        G_composite         = V_comp + S_comp;
//                        N_composite         = Vnet * Rho * avogadro;
//                        add_to_N (N_composite, G_composite, Rg, Rb, dN, Nmin, len_N, NArray_Gmin, NArray_confs, NArray_quant);
//
//
//                        //CASE 2: Both centres lie on the good patch
//                        Cb1[0] = -(p_width/2.0) + db;       //other two components same as top
//                        Cb2[0] = (p_width/2.0) - db;
//
//                        if(Cb1[0] < 0.0 && 0.0 < Cb2[0])    //(cb1<cg<cb2): ensures that no crystal intersects with the other patch
//                        {
//                            Sphere sphere_B1 (Rb, Cb1);
//                            Sphere sphere_B2 (Rb, Cb2);
//                            Spheres[1] = sphere_B1;
//                            Spheres[2] = sphere_B2;
//
//                            Composite_cluster Cluster(Spheres, wall);
//                            std::vector<double> net_measures = Cluster.get_net_measures();
//                            std::vector<double> proj_SAs     = Cluster.get_proj_SA();
//                            double Vnet                = net_measures[0];
//                            double Snet                = net_measures[1];
//                            double V_comp              = -1.0 * Rho * Mu * Vnet;
//                            double S_comp              = Sigma*(Snet - cos(theta_cg)*(proj_SAs[0]) - cos(theta_cb)*(proj_SAs[1] + proj_SAs[2]));
//                            G_composite         = V_comp + S_comp;
//                            N_composite         = Vnet * Rho * avogadro;
//                            add_to_N (N_composite, G_composite, Rg, Rb, dN, Nmin, len_N, NArray_Gmin, NArray_confs, NArray_quant);
//
//                        }
//                        else
//                        {
//                            continue;
//                        }
//                    }
//                    else    //cb +- rb is bigger than 3a/2
//                    {
//                        break;
//                    }
//
//                }
//            }
//
//        }
//        else    //rg is bigger than 3a/2
//        {
//            break;
//        }
//    }
//
//    return 0;
//}
