//
//  test_single_object_parallelism.cpp
//  Clusters
//
//  Created by Sarwar Hussain on 5/25/21.
//  Copyright © 2021 Sarwar Hussain. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include "manager_spherical_caps.hpp"
#include <mpi.h>


using namespace std;

//class TestInterface
//{
//public:
//    TestInterface(){}
//    ~TestInterface(){}
//    virtual void increment_var(){};
//    virtual void print_var (){};
//};
//
//class TestClass:public TestInterface
//{
//
//public:
//    TestClass(){instance_variable= 0;}
//    TestClass(int my_variable){instance_variable = my_variable;}
//    ~TestClass(){}
//    void increment_var() {instance_variable++;}
//    void print_var () {printf("var=%d\n", instance_variable);}
//protected:
//    int instance_variable;
//};

int myRank, nProcs;
int main(int argc, char * argv[])
{
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    
    int branches_per_node, Levels;
    if(myRank == 0)
    {
        branches_per_node = atoi(argv[1]);
        Levels = atoi(argv[2]);
        printf("argc = %d branches_per_node = %d\n Levels=%d\n", argc, branches_per_node, Levels );
    }
    MPI_Bcast (&branches_per_node, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast (&Levels, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    
    ManagerSphericalCaps manager (branches_per_node,  Levels);

    manager.setup(argv, 2);
    int total_args_read = 2 + manager.get_read_args();
    int profiling = 0; //whether this is the profiling step or not.
    int constant_modifier_run=0; //whether to do the simulation with constant clstr centre modifiers {-1 for bad patches, +1 for good patches}. When this is on levels != n_patches, instead levels=n_unique_patches
    int evolve_finished;
    if(argc > total_args_read+1) //because argc also has the name of the program
    {
        if(myRank == 0 )
        {
            profiling  = atoi(argv[total_args_read + 1]);
        }
        MPI_Bcast (&profiling, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    }
    if(argc > total_args_read+2)
    {
        if(myRank == 0 )
        {
            constant_modifier_run  = atoi(argv[total_args_read + 2]);
        }
        MPI_Bcast (&constant_modifier_run, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    }
    
    //    evolve_finished=0;
    //        EvolveSphericalCap evolve_cap = manager.get_evolve_spherical_cap();
    //        double Rg, Rb, Rg_second, Rb_second;
    //        Rg=2.800000000000000;
    //        Rb=4.418144406900000;
    //        Rg_second=1.083205120600000;
    
    //    if(myRank==0)
    //    {
    //        Rg  = atof(argv[total_args_read + 3]);
    //        Rb  = atof(argv[total_args_read + 4]);
    //        Rg_second  = atof(argv[total_args_read + 5]);
    //        Rb_second  = atof(argv[total_args_read + 6]);
    //    }
    //    MPI_Bcast (&Rg, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    //    MPI_Bcast (&Rb, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    //    MPI_Bcast (&Rg_second, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    //    MPI_Bcast (&Rb_second, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);
    //
    //    std::vector<double> input_Radii {Rg, Rb, Rg_second, Rb_second};
    //    evolve_cap.artificial_evolve(input_Radii);
    //
    //    printf("process=%d is at the barrier\n", myRank);
    
    
    if(profiling == 0)
    {
        evolve_finished = manager.evolve(constant_modifier_run);
    }
    else
    {
        evolve_finished = manager.evolve_profiling(constant_modifier_run);
        printf("out of evolve_profiling in manager rank=%d\n", myRank);
    }
    
    
    
                    //int evolve_finished = manager.evolve();

                //    if(myRank == 0)
                //    {
                //        printf("Before dummy evolve\n");
                //        manager.dummy_evolve(2);
                ////        manager.print_surface_ptr();
                ////        manager.print_box();
                ////        manager.print_mc_and_check_boundary();
                //    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    if(evolve_finished == 1)
    {
//        manager.gather();
        if(myRank == 0){printf("Gather done\n");}
//        manager.open_output_file();
        //printf("Before manager printing to file\n");
        if(profiling == 0){manager.print_to_file();}
    }
    else
    {
        if(myRank == 0){printf("Global evolve hasn't returned 1\n");}
    }
    
    manager.free_MPI_comms();
    MPI_Finalize();
    return EXIT_SUCCESS;
}













//int BranchSize;
////    if(myRank == 0)
////    {
////        printf("nProcs = %d\n", nProcs);
//BranchSize = atoi(argv[1]);
//printf("BranchSize = %d\n", BranchSize);
////    }
////    MPI_Bcast (&BranchSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//TestInterface* test_interface;
//
//
////    if(myRank == 0)
////    {
//if(BranchSize == 2)
//{
//    printf("Inside branchsize =2\n");
//    //            TestClass* test_class_ptr =
//    test_interface = new TestClass (BranchSize);
//}
//else
//{
//    printf("Inside branchsize not 2\n");
//    //            TestClass* test_class_ptr =
//    test_interface = new TestClass;
//}
//test_interface->increment_var();
//test_interface->print_var ();
//delete test_interface;
////    }

