# GS-Mixing
ALICE Event Mixing for NTuples using Gale-Shapely Algorithm
To be used with a 13e_clusterv2_small.root file

Please ignore the majority of the functions in this code, they were a folly attempt of an eager graduate student to quickly do some physics. The functions to look at are PTable, StableMarriage(), and PairingTest() (with Centrality list and entries_after_cuts as helper functions).

If for some reason you would like to run the code, load a .root file of your choosing, and do TTree->MakeClass() for a skeleton .h and .C file. Port over only the declaration of functions in the .h file and the whole body of code in the .C file. Also be sure to pay attention to the various #include's.

The code works as intended, and can be used with a root file. However, it must be readapted because it is terribly slow in it's current implementation.

This implementation will likeley be abondoned, but I will use what I learned from this and work with a post doc's framework with parallelization. This implemntation includes 1:Many pairing of events based on multiplicity, ensures an event is not used 'n' number of times. Because I did not think of the simple solution of replicating events, I introduced a status called "unfullfilled" when an event is not able to be paired with a specified number of events.

Hopefully my second attempt, this time collaborating with a more experienced coder will be more successfull. As my first major project in C++, I at least learned quite a bit.
