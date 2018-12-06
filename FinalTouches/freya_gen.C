#include "/Users/jeffreyburggraf/PycharmProjects/FREYA/MyFREYA/ROOT_tree_writer.C" // machine specific directory
#include <cstdlib>
#include <iostream>

FissionEventWriter  writer = FissionEventWriter();

void example(){
    int isotope = 98252;         // zzaaa
    double  time = 0;           // initial time
    double nubar = 3.7;            // mean neutron multiplicity is user specified.
    double erg = 10;            // energy of incident particle
    int fiss_type = 0;          // 0: SF    1: neutron-induced    2: photon-induced
    double dir[3] = {1,0,0};   // direction of incident particle

	for (int i=0; i<10000; i++){
        fissionEvent  * fe = new fissionEvent(isotope, time, nubar, erg, fiss_type, dir);
        writer.WriteEvent(fe);  // Write event to TTree. FREYA.root will appear in the present working directory.

    }

	TBrowser *tb = new TBrowser();

}

