# DyonicKNGeodesics
Implementation of analytic solutions to generic bound and plunging dyonic geodesics in Kerr Newmann, in terms of Jacobi elliptic functions

Consists of two primary packages  within Kernel, DyonicKNGeoOrbit which solves for bound orbits and, DyonicKNGeoPlunge which solves for plunging geodesics.
Packages can be loaded as seen in GeodesicRunFile as a paclet using the commands 
    - SetDirectory@NotebookDirectory[]
    - PacletDirectoryLoad["."]
    - << DyonicKNGeodesicsdfd
    
Each file has has one function, DyonicKNGeoPlunge[...] and DyonicKNGeoOrbit[...] output with certain associations
    
Associations:     ["a" -> a,
                  "Parametrization"->"Mino", 
                  "ChargeParameters",
                  "Energy", 
                  "AngularMomentum", 
                  "CarterConstant", 
                  "ConstantsOfMotion",
                  "Periods",
                  "RadialRoots",
                  "PolarRoots",
                  "EllipticBasis",
                  "RadialRootCrossingTimeMino",
                  "Trajectory"]
                  
Examolpe Usage of of these function can be found in GeodesicRunFile.nb, additional plots found in arxiv submission can we generated in Paper_Plots.nb
