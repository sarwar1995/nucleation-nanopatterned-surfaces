# nucleation-nanopatterned-surfaces
This C++ codebase is parallelizable with 2^n cores. It is used to generate spherical caps on striped surfaces with alternating crystal-philic and crystal-phobic surfaces.

The important code files are as follows:
1. run_manager_spherical_caps.cpp is the driver code that sets up the `ManagerSphericalCaps` class
2. The `ManagerSphericalCaps` (manager_spherical_caps.hpp) class sets up the following objects:
  * ParallelProcess parallel_process
  *  `EvolveSphericalCap evolve_spherical_cap `
  *  `CheckBoundary check_boundary`
  * ` MC mc_engine`
  *  `DynamicBox dynamic_box`
  *  `Stripes stripes`
  *  `Surface* surface_ptr`
  *  `PeriodicIO periodic_io`
3. The `ParallelProcess` (parallel_process.hpp) takes care of the parallelization scheme and is used to provide information about the topology of the processors to other classes. This topology is the distribution of calculation tasks (composite geometrical cluster's configuration) to processors and not their physical topology within the HPC supercomputing cluster.
4. The `MC` (monte-carlo engine) class takes in a pointer to the `Shape` class and calculates the volume and surface area. There is a prallelized version of this implementation in  MC_parallel.hpp which is what is used in the manager.
5. The `EvolveSphericalCap` class (evolve_spherical_caps.hpp) is responsible for starting the simulation (growth of cluster) and printing the output periodically from all processors. 
It also uses the `CheckBoundary` class classes to check the validity of each cap i.e. whether a patch's boundary is breached by a cap growing on top of it, in which case new caps are added to keep the contact angle valid on the new patch (the patch onto which a cap spills after breaching the original patch on which it was growing)
6. The `DynamicBox` class is used to grow the box around the growing cluster, by adding extra points (with constant density) so that Volume calculations can be done efficiently by avoiding unnecessary checking of points that would be needed in case of a fixed large box for identifying points that are inside a very small cluster as the majority of points in a large box woule be outside the cluster. By keeping the box as tightly bounded around the cluster as possible, we avoid this.
7. The actual geometrical components are defined in Shape.hpp, Sphere.hpp, Spherical_cap.hpp and Composite_cluster.hpp. Shape.hpp defines the interface for
implementing isInside and nearSurf methods to identify where a point lies with respect to a shape, which are the main methods used by `MC` class to calculate Volume/Surface Area.
8. The `Surface` and `Stripes` classes hold the geometry of the patterned surface on which the cluster grows and they interface with the `CheckBoundary` class to provide information about a clusters footprint on the surface, in terms of boundary crossing.

There are several other files in this repo, like SpheroCylinder.hpp and SphericalCapDivided.hpp that were used to test a slight modification of this theory. Other files are mostly ad-hoc ways of conducting unit and integration testing.
