//
//  ExtraCodeMain.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 2/12/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//


//std::string outFileName_nearsurf_not_spreadout_points_1 = tag + "near_surf_points_not_spreadout_1.txt";
//std::string outFileName_nearsurf_not_spreadout_points_2 = tag + "near_surf_points_not_spreadout_2.txt";
//std::string outFileName_surfNgbs_not_spreadout_points_1 = tag + "near_surfNgb_points_not_spreadout_1.txt";
//std::string outFileName_surfNgbs_not_spreadout_points_2 = tag + "near_surfNgb_points_not_spreadout_2.txt";
//
//std::string outFileName_nearsurf_spreadout_points_1 = tag + "near_surf_points_spreadout_1.txt";
//std::string outFileName_nearsurf_spreadout_points_2 = tag + "near_surf_points_spreadout_2.txt";
//std::string outFileName_surfNgbs_spreadout_points_1 = tag + "near_surfNgb_points_spreadout_1.txt";
//std::string outFileName_surfNgbs_spreadout_points_2 = tag + "near_surfNgb_points_spreadout_2.txt";

//                            std::vector<std::vector<double> > ngb_surface_points ;
//This printing part should be after the calc_volume_SA call
//                            if(j == 1)
//                            {
//                                ngb_surface_points = cell_list.calc_neighbors (mc_engine.surface_points);
//                                if(dB_list[k] == -1)
//                                {
//                                    fprintf(pointsFile_not_spreadout_1,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list[k]);
//                                    mc_engine.print_surf_points(pointsFile_not_spreadout_1);
//                                    print_points_tofile(ngb_surface_points, NgbpointsFile_not_spreadout_1);
//                                }
//
//                                if(dB_list[k] == 1)
//                                {
//                                    fprintf(pointsFile_spreadout_1,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list[k]);
//                                    mc_engine.print_surf_points(pointsFile_spreadout_1);
//                                    print_points_tofile(ngb_surface_points, NgbpointsFile_spreadout_1);
//                                }
//                            }
//
//                            if(j == 2)
//                            {
//                                ngb_surface_points = cell_list.calc_neighbors (mc_engine.surface_points);
//                                if(dB_list[k] == -1)
//                                {
//                                    fprintf(pointsFile_not_spreadout_2,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list[k]);
//                                    mc_engine.print_surf_points(pointsFile_not_spreadout_2);
//                                    print_points_tofile(ngb_surface_points, NgbpointsFile_not_spreadout_2);
//                                }
//                                if(dB_list[k] == 1)
//                                {
//                                    fprintf(pointsFile_spreadout_2,"Rg=%10.10f\t Rb=%10.10f\t db=%d\n", Rg, Rb, dB_list[k]);
//                                    mc_engine.print_surf_points(pointsFile_spreadout_2);
//                                    print_points_tofile(ngb_surface_points, NgbpointsFile_spreadout_2);
//                                }
//                            }


//pointsFile_not_spreadout_1 = fopen(outFileName_nearsurf_not_spreadout_points_1.c_str(), "w");
//pointsFile_not_spreadout_2 = fopen(outFileName_nearsurf_not_spreadout_points_2.c_str(), "w");
//NgbpointsFile_not_spreadout_1 = fopen(outFileName_surfNgbs_not_spreadout_points_1.c_str(), "w");
//NgbpointsFile_not_spreadout_2 = fopen(outFileName_surfNgbs_not_spreadout_points_2.c_str(), "w");
//
//pointsFile_spreadout_1 = fopen(outFileName_nearsurf_spreadout_points_1.c_str(), "w");
//pointsFile_spreadout_2 = fopen(outFileName_nearsurf_spreadout_points_2.c_str(), "w");
//NgbpointsFile_spreadout_1 = fopen(outFileName_surfNgbs_spreadout_points_1.c_str(), "w");
//NgbpointsFile_spreadout_2 = fopen(outFileName_surfNgbs_spreadout_points_2.c_str(), "w");


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
