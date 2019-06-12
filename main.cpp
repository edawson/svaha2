#include <vector>
#include <stdio>
#include "tinyFA.hpp"
#include "sparsepp/spp.h"
#include "gfakluge.hpp"

int construct_contig(){

}

int construct(){

}

int make_breakpoints(){

}

void usage(){
    cout << "svaha2: linear-time, low-memory construction of variation graphs." << endl;
    cout << "Usage: svaha2 [options] -r <ref.fa> -v <variants.vcf> " << endl;
    cout << "options:" << endl;
    cout << "-m / --max-node-size : maximum length (in basepairs) of a node in the graph." << endl;
    cout << "-t / --threads       : number of OMP threads to use in construction." << endl;
    cout << "version 0.1" << endl;
}

int main(int argc, char** argv){

    char* ref_file = NULL;
    char* var_file = NULL;
    int threads = 1;
    int max_node_size = 128;

    int c;
    int optind = 2;

    if (argc <= 2){
        usage();
        exit(1);
    }

    if (ref_file == NULL){
        cerr << "ERROR: svaha2 requires a reference file." << endl;
        usage();
        exit(1);
    }
    
    

    return 0;
}
