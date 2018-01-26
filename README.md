# GS-Mixing
ALICE Event Mixing for NTuples using Gale-Shapely Algorithm
To be used with a 13e_clusterv2_small.root file

If for some reason you would like to run the code, load a .root file of your choosing, and do TTree->MakeClass() for a skeleton .h and .C file. Port over only the declaration of functions in the .h file and the whole body of code in the .C file. Also be sure to pay attention to the various #include's.

The code works as intended, and can be used with a root file. However, it must be readapted because it is terribly slow in it's current implementation.

This implementation will likeley be abondoned, but I will use what I learned from this and work with a post doc's framework that utilizes parallelization. This implemntation includes 1-to-Many pairing of events based on multiplicity (primary vertex data was not available at the time) and ensures an event is not used 'n' number of times. Because I did not think of the simple solution of replicating events, I introduced a status called "unfullfilled" when an event is not able to be paired with a specified number of events.

Hopefully my second attempt, this time collaborating with a more experienced coder, will be more successfull. As my first major project in C++, I at least learned quite a bit.
