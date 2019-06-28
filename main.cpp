#include <iostream>
#include <vector>
#include <getopt.h>
#include <omp.h>
#include <map>
#include <set>
#include "pliib.hpp"
#include "tinyfa.hpp"
#include "tinyVCF.hpp"
#include "spp.h"
#include "gfakluge.hpp"

using namespace std;

namespace svaha {

struct node{
    std::uint64_t id;
    char* contig = NULL;
    char** paths = NULL;
    std::uint8_t num_paths;
    char* seq = NULL;
    std::uint32_t seqlen;
    std::uint64_t contig_seq_start;
    std::uint64_t contig_seq_end;
    // std::vector<node*> prev;
    // std::vector<node*> next;
    std::string emit(){
        gfak::sequence_elem s;
        s.id = id;
        s.name = std::to_string(this->id);
        s.length = this->seqlen;
        s.sequence = string(this->seq);
        return s.to_string_2();
    };
};

struct edge{
    // These two booleans define if the edges
    // go in the default from_end/to_start way or not.
    // Inversions don't follow this pattern,
    // instead going from the start to the end, for example.
    edge(){

    };
    edge(svaha::node*& first,
        svaha::node*& second,
        bool from_end = true,
        bool to_start = true,
        bool is_reverse = false){
            from_node = first;
            to_node = second;
            this->from_end = from_end;
            this->to_start = to_start;
            this->is_reverse = is_reverse;
    };
    bool from_end = true;
    bool to_start = true;
    bool is_reverse = false;
    svaha::node* to_node;
    svaha::node* from_node;
    std::string emit(){
        gfak::edge_elem e;
        e.source_name = std::to_string(from_node->id);
        e.sink_name = std::to_string(to_node->id);
        e.source_orientation_forward = from_end;
        e.sink_orientation_forward = to_start;
        e.ends.set(0,1);
        e.ends.set(1,1);
        e.ends.set(2,0);
        e.ends.set(3,0);
        e.type = 1;
        return e.to_string_2();
    };
};


struct path{
    std::uint64_t id;
    std::string name;
    std::vector<edge*> edges;
    std::vector<node*> nodes;
};

struct pre_contig{
    char* seq;
    std::uint32_t seqlen;
    std::uint64_t c_node_id = 0;
    std::uint64_t c_edge_id = 0;
    //spp::sparse_hash_map<uint64_t, node*> bp_to_node;
    svaha::node* bp_to_node;
    spp::sparse_hash_map<uint64_t, TVCF::variant> bp_to_variant;
    spp::sparse_hash_map<uint64_t, TVCF::variant> interchrom_variants;
    std::vector<std::uint64_t> breakpoints;
    // std::vector<edge*> edges;
    // std::vector<node*> nodes;

    node* create_node(){
        node* n = new node();
        n->id = ++c_node_id;
        return n;
    };
    void bump_node_ids(std::uint64_t bumpsz){
        c_node_id = bumpsz;
    };
    void set_edge_id_start(){
        this->c_edge_id = this->c_node_id + 1;
    };
    void add_bp(std::uint64_t b){
        breakpoints.push_back(b);
    };
    void sort_bps(){
        std::sort(breakpoints.begin(), breakpoints.end());
    };
    void uniquify_bps(){
        std::set<std::uint64_t> unis (breakpoints.begin(), breakpoints.end());
        std::vector<std::uint64_t> bps (unis.begin(), unis.end());
        breakpoints = bps;
    };
};


struct graph{
    uint64_t min_node_id = 0;
    uint64_t max_node_id = 0;
    uint64_t min_edge_id = 0;
    uint64_t max_edge_id = 0;
    uint64_t curr_node_id = 0;
    uint64_t curr_edge_id = 0;
    node* create_node(){
        node* n = new node();
        n->id = ++curr_node_id;
        return n;
    };
    std::uint64_t edge_id(){
        return ++curr_edge_id;
    };
    spp::sparse_hash_map<string, pre_contig> name_to_contig;
    //spp::sparse_hash_map<string, vector<TVCF::variant>> name_to_variants;
    void re_id(){
        std::uint64_t prev_id = 0;
    };
};

void make_maxnode_breakpoints(uint64_t seqlen, uint32_t maxnodelen, std::vector<std::uint64_t>& ret){
    ret.reserve(seqlen / maxnodelen);
    for (size_t i = maxnodelen; i < seqlen; i += maxnodelen){
        ret.push_back(i);
    }
};


};

void usage(){
    cout << "svaha2: linear-time, low-memory construction of variation graphs." << endl;
    cout << "Usage: svaha2 [options] -r <ref.fa> -v <variants.vcf> " << endl;
    cout << "options:" << endl;
    cout << "-m / --max-node-size : maximum length (in basepairs) of a node in the graph." << endl;
    cout << "-r / --reference     : The reference genome to use for construction." << endl;
    cout << "-v / --vcf           : A VCF file to use for construction." << endl;
    cout << "-b / --bed           : a bed file to restrict chromosomes to." << endl;
    cout << "-T / --no-translocations  : ignore interchromosomal variants." << endl;
    //cout << "-t / --threads     : number of OMP threads to use in construction." << endl;
    cout << "-f / --flat          : use flat alternates (every allele is represented by at least one node)." << endl;
    cout << "-p / --paths         : output path information for variants." << endl;
    //cout << "-I / --insertions         : FASTA file of insertion variant sequences." << endl;
    cout << "version 0.1" << endl;
}

int main(int argc, char** argv){

    char* ref_file = NULL;
    char* var_file = NULL;
    char* bed_file = NULL;
    char* insertion_file = NULL;
    bool flat = false;
    bool do_translocations = true;
    int threads = 1;
    int max_node_size = 128;

    int c;
    int optind = 2;

    if (argc <= 2){
        usage();
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},
            {"bed", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"max-node-size", required_argument, 0, 'm'},
            {"flat", no_argument, 0, 'f'},
            {"no-translocations", no_argument, 0, 'T'},
            {"paths", no_argument, 0, 'p'},
            {"insertions", no_argument, 0, 'I'},
            {0,0,0,0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "hr:v:b:TI:t:m:pf", long_options, &option_index);
        if (c == -1){
            break;
        }
        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'r':
                pliib::strcopy(optarg, ref_file);
                break;
	    case 'b':
		pliib::strcopy(optarg, bed_file);
		break;
            case 'I':
                pliib::strcopy(optarg, insertion_file);
                break;
            case 'v':
                pliib::strcopy(optarg, var_file);
                break;
            case 'T':
                do_translocations = false;
                break;
            case 'm':
                max_node_size = atoi(optarg);
                break;
            case 'f':
                flat = true;
                break;
            case 'h':
                usage();
                exit(1);
            default:
                abort();
        }
    
    }

    if (ref_file == NULL){
        cerr << "ERROR: svaha2 requires a reference file." << endl;
        usage();
        exit(1);
    }


    TFA::tiny_faidx_t tf;
    if (!TFA::checkFAIndexFileExists(ref_file)){
        cerr << "No .fai index found for " << ref_file << "; creating index." << endl;
        TFA::createFAIndex(ref_file, tf);
        TFA::writeFAIndex(ref_file, tf);
        cerr << "FASTA index created." << endl;

    }
    else{
        cerr << "FASTA index found for " << ref_file << endl;
        parseFAIndex(ref_file, tf);
        cerr << "Parsed fasta index with " << tf.seq_to_entry.size() << " entries." << endl;
    }

    TFA::tiny_faidx_t insertion_tf;
    if (insertion_file != NULL){
        if (TFA::checkFAIndexFileExists(insertion_file)){
            cerr << "No .fai index found for " << insertion_file << "; creating index." << endl;
            TFA::createFAIndex(insertion_file, insertion_tf);
            TFA::writeFAIndex(insertion_file, insertion_tf);
            cerr << "FASTA index created." << endl;
        }
        else{
            cerr << "FASTA index found for " << insertion_file << endl;
            parseFAIndex(insertion_file, insertion_tf);
            cerr << "Parsed fasta index with " << tf.seq_to_entry.size() << " entries." << endl;
        }
    }

    // Step -1: Create a bedfile mask if needed.
    std::set<std::string> acceptable_chroms;
    if (bed_file != NULL){
   	    std::ifstream bfi(bed_file);
	    if (bfi.good()){
		    char* line = new char[100];
		    while (bfi.getline(line, 100)){
			    std::string s(line);
			    pliib::strip(s, '\n');
			    pliib::strip(s, ' ');
			    acceptable_chroms.insert(s);
		    }
	    }
	    else{
		    cerr << "Error: no bed file found [ " << bed_file << " ]." << endl;
		    exit(1);
	    }
    }

    // Step 0: initialize a graph to hold contig information
    svaha::graph sg;
    for (auto& x : tf.seq_to_entry){
        svaha::pre_contig p;
	    std::string name(x.first);
	    if (acceptable_chroms.size() == 0 || acceptable_chroms.find(name) != acceptable_chroms.end()){
            sg.name_to_contig[name] = p;
            //std::vector<TVCF::variant> bps;
            //sg.name_to_variants[name] = bps;
	    }
    }

    // Step one: add breakpoints for max-node-length
     

    std::uint64_t var_count = 0;
    
    if (var_file != NULL){
        std::ifstream vfi(var_file);
        if (vfi.good()){
            cerr << "Loading variants from " << var_file << endl;
            
            char* current_contig = NULL;
            //std::vector<std::uint64_t> breakpoints;
            int mxlsz = 2000;
            char* line = new char[mxlsz];
            TVCF::variant* var = new TVCF::variant();
            while (vfi.getline(line, mxlsz)){
                if (line[0] != '#'){
                    TVCF::parse_variant_line(line, var);
                    if (current_contig != NULL && current_contig != var->chrom){
                        // We've moved on to a new contig;
                        // we need to write / flush the local copy of breakpoints
                        //sg.name_to_contig.at(string(current_contig)).breakpoints = breakpoints;
                        //breakpoints.clear();
                    }
                    // Add the breakpoint triggered by POS
                    //cerr << "Creating breakpoint at [" << var->zero_based_position() <<
                    //    ", " << var->get_reference_end(0) << "] for variant type " << var->get_info("SVTYPE") << endl;

                    if (TFA::hasSequence(tf, var->chrom)){
                        if (acceptable_chroms.size() == 0 ||
                        acceptable_chroms.find(string(var->chrom)) != acceptable_chroms.end()){
                            
                            sg.name_to_contig.at(string(var->chrom)).breakpoints.push_back(var->zero_based_position());
                            string svtype = var->get_sv_type();
                            TVCF::variant v(*var);
                            //sg.name_to_variants.at(string(var->chrom)).push_back(v);
                            sg.name_to_contig.at(string(var->chrom)).bp_to_variant[v.zero_based_position()] = v;
                            if (svtype != "TRA"){
                                sg.name_to_contig.at(string(var->chrom)).breakpoints.push_back(var->get_reference_end(0));
                            }
                            else if (do_translocations){
                                std::string c2 = var->get_info("CHR2");
                                if (acceptable_chroms.size() == 0 || acceptable_chroms.find(c2) != acceptable_chroms.end()){
                                    sg.name_to_contig.at(c2).breakpoints.push_back(var->get_reference_end(0));
                                }
                            }

                        }
                        else{
                            continue;
                        }
                    }
                    else{
                        cerr << "ERROR: chrom not found: " << var->chrom  << " or " << var->get_info("CHR2") << "; are you using a VCF and FASTA from the same reference?" << endl;
                        cerr << "EXITING" << endl;
                        exit(9);
                    }
                    var_count += 1;
                    
                }
            }

            delete [] line;
            delete var;
            vfi.close();
        }
        else{
            cerr << "ERROR: could not open variant file: " << var_file << endl;
        }
    }

    cerr << "Parsed " << var_count << " variants." << endl;
    cerr << "Created map of contig to breakpoints from variants." << endl;

    for (auto c : sg.name_to_contig){
        cerr << "Contig " << c.first << " has " << c.second.breakpoints.size() << " breakpoints." << endl;
    }

    cerr << endl << endl;
    cerr << "Adding max-node-length breakpoints." << endl;

    for (auto& c : sg.name_to_contig){
        TFA::getSequenceLength(tf, c.first.c_str(), c.second.seqlen);
        svaha::make_maxnode_breakpoints( c.second.seqlen, max_node_size, c.second.breakpoints);
        c.second.breakpoints.push_back(c.second.seqlen);
        cerr << "Contig " << c.first << " has " << c.second.breakpoints.size() << " breakpoints." << endl;
    }
    
    // Sort breakpoints and
    // Uniquify breakpoints.
    for (auto& c : sg.name_to_contig){
        c.second.sort_bps();
        c.second.uniquify_bps();
    }
    cerr << "Sorted and uniquified breakpoints" << endl;

    std::uint64_t greatest_prev_id = 0;

    for (auto& c : sg.name_to_contig){

        // Get reference sequence
        if (! TFA::hasSequence(tf, c.first.c_str())){
            cerr << "No sequence found for " << c.first << endl;
            exit(1);
        }
        TFA::getSequence(tf, c.first.c_str(), c.second.seq);
        std::size_t numbp = c.second.breakpoints.size();
        
        //std::vector<TVCF::variant> contig_vars = sg.name_to_variants.at(c.first);
        std::vector<std::uint64_t> bps = c.second.breakpoints;
        c.second.bp_to_node = new svaha::node [c.second.seqlen];
        bool snptrip = false;
        
        std::uint64_t pos = 0;
        std::uint64_t last_ref_node_pos = 0;
        svaha::node* prev_ref_node;

        //std::uint64_t total_seq = 0;

        for (size_t i = 0; i < bps.size() - 1; ++i){

            svaha::node* n = sg.create_node();
            pliib::strcopy(c.first.c_str(), n->contig);
            std::uint64_t brk = bps[i];
	    //cerr << c.second.seqlen << " " << pos << " " << brk - pos << endl;
            pliib::strcopy(c.second.seq + pos,  brk - pos, n->seq);
            n->seqlen = brk - pos;

            // Emit the node, caching it if we need it later for a variant.
            cout << n->emit() << endl;
            //cout << n->id << " " << n->contig  << " " << n->seqlen << endl;
            if (c.second.bp_to_variant.find(pos) != c.second.bp_to_variant.end()){
                c.second.bp_to_node[pos] = *n;
                c.second.bp_to_node[brk - 1] = *n;
            }
            // total_seq += brk - pos;
            // cout << total_seq << " " << c.second.seqlen <<endl;




            //cout << "pos: " << pos << " brk: " << brk << endl;


            

            TVCF::variant v;
            string vtype;
            std::uint64_t vpos;
            std::uint64_t vend;
            if (c.second.bp_to_variant.find(brk) != c.second.bp_to_variant.end()){
                v = c.second.bp_to_variant.at(brk);
                vtype = v.get_info("SVTYPE");
            }
            else if (c.second.interchrom_variants.find(brk) != c.second.interchrom_variants.end()){
                v = c.second.interchrom_variants.at(brk);
                vtype = v.get_info("SVTYPE");
            }

            if (vtype == "DEL"){
                // Get the start and end of the variant
            }
            else if (vtype == "INS"){
                svaha::node* ins_node = sg.create_node();
            }
            else if (vtype == "INV"){
                // Get start and end of variant,
                // then get the corresponding nodes.
            }

                    // IFF the variant is an INS type:
        //      get its sequence from its seq field OR an external FASTA
        
        // IFF the variant is a DEL or INV type:
        //      add two edges to the nodes defining the variant,
        //      2+ nodes in the case of FLAT alts
        //      and 1+i nodes in the case of graphical alts
                    
        // IFF the variant is an INV type:
        //      wire up its end to the previous END 
        //      and wire up its start to the next start.
        // IFF the variant is a TRA / breakend
        //      emit the right edges to link from chrom to chr2

        // Wire up ref nodes on the backbone if we got 'em
        // if (pos != 0 && 
        //     c.second.bp_to_node[pos] != NULL &&
        //     c.second.bp_to_node[last_ref_node_pos] != NULL &&
        //     c.second.bp_to_node[pos]->contig != NULL &&
        //     c.second.bp_to_node[last_ref_node_pos]->contig != NULL){
        //     svaha::edge e(c.second.bp_to_node[last_ref_node_pos], c.second.bp_to_node[pos]);
        //     cout << e.emit() << endl;
        // }

        if (pos != 0 && 
          n->contig != NULL && 
          prev_ref_node->contig != NULL ){
              svaha::edge e(prev_ref_node, n);
              cout << e.emit() << endl;
          }
        //  if (pos != 0){
        //      cout << c.second.bp_to_node[pos]->id << endl;
        //      cout << c.second.bp_to_node[last_ref_node_pos]->id << endl;
        //  }
        last_ref_node_pos = pos;
        prev_ref_node = n;
        pos = brk;
    }
    bps.clear();
    c.second.breakpoints.clear();
    delete [] c.second.seq;

    }

    // Handle interchromosomal variants
    for (auto& c : sg.name_to_contig){
        for (auto& pv : c.second.interchrom_variants){
            try{
                svaha::node* first;
                first = ( sg.name_to_contig.at(std::string(pv.second.chrom)).bp_to_node + pv.second.zero_based_position() );
                svaha::node* second;
                second = (sg.name_to_contig.at(std::string(pv.second.get_info("CHR2"))).bp_to_node + std::stoull(pv.second.get_info("END")));
            }
            catch(const std::out_of_range& oor){

            }

        }
    }

    return 0;
}
