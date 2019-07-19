#include <iostream>
#include <vector>
#include <getopt.h>
#include <omp.h>
#include <map>
#include <set>
#include "pliib.hpp"
#include "tinyFA.hpp"
#include "tinyVCF.hpp"
#include "spp.h"
#include "gfakluge.hpp"

using namespace std;

namespace svaha {


    struct node{
        std::uint64_t id;
        char* contig = nullptr;
        char** paths = nullptr;
        std::uint8_t num_paths = 0;
        std::uint8_t pathcap = 0;
        char* seq = nullptr;
        std::uint32_t seqlen;
        std::uint64_t contig_seq_start;
        std::uint64_t contig_seq_end;
        // std::vector<node*> prev;
        // std::vector<node*> next;


        ~node(){
            pliib::strdelete(contig);
            pliib::strdelete(seq);
        }
        std::string emit(){
            gfak::sequence_elem s;
            s.id = id;
            s.name = std::to_string(this->id);
            s.length = this->seqlen;
            s.sequence = string(this->seq);
            return s.to_string_2();
        };

        void create_paths(std::uint8_t np = 3){
            char** xpaths = new char*[np];
            for (std::uint8_t i = 0; i < num_paths; ++i){
                pliib::strcopy(paths[i], xpaths[i]);
            }
            pathcap = np;
        };
        void add_path(const char* pname){
            if (num_paths + 1 > pathcap){
                create_paths();
            }
            pliib::strcopy(pname, paths[num_paths - 1]);
        };
    };

    struct edge{
        // These two booleans define if the edges
        // go in the default from_end/to_start way or not.
        // Inversions don't follow this pattern,
        // instead going from the start to the end, for example.
        bool from_end = true;
        bool to_start = true;
        bool is_reverse = false;
        svaha::node* to_node;
        svaha::node* from_node;
        edge(){
            from_node = nullptr;
            to_node = nullptr;
        };
        ~edge(){

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
            return e.to_string_1();
        };
    };

    struct path_occ_t{
        svaha::node* node = nullptr;
        std::uint32_t rank = 0;
        char* pname = nullptr;
        bool isforward = true;
        inline string to_string(){
            if (node != nullptr && pname != nullptr){
                std::ostringstream st;
                st << "W" << "\t" <<
                pname << "\t" <<
                rank << "\t" <<
                node->id << "\t" <<
                (isforward ? "+" : "-");
                return st.str();
            }
            return "";
        }
    };


    struct path{
        std::uint64_t id;
        std::string name;
        //std::vector<edge*> edges;
        std::vector<path_occ_t> nodes;
    };

    struct pre_contig{
        char* seq = nullptr;
        std::uint32_t seqlen;
        std::uint64_t c_node_id = 0;
        std::uint64_t c_edge_id = 0;
        //spp::sparse_hash_map<uint64_t, node*> bp_to_node;
        svaha::node** bp_to_node;
        //spp::sparse_hash_map<uint64_t, TVCF::variant*> bp_to_variant;
        spp::sparse_hash_map<uint64_t, TVCF::minimal_allele_t*> bp_to_allele;
        spp::sparse_hash_map<uint64_t, node*> bp_to_inserted_node;
        //spp::sparse_hash_map<uint64_t, TVCF::variant*> interchrom_variants;
        spp::sparse_hash_map<uint64_t, TVCF::minimal_allele_t*> bp_to_interchrom;
        std::vector<std::uint64_t> breakpoints;

        ~pre_contig(){
            pliib::strdelete(seq);
            if (bp_to_node != nullptr){
                //delete [] bp_to_node;
            }
        };
        void clear_seq(){
            pliib::strdelete(seq);
            this->seq = nullptr;
        };
        // std::vector<edge*> edges;
        // std::vector<node*> nodes;

        // node* create_node(){
        //     node* n = new node();
        //     n->id = ++c_node_id;
        //     return n;
        // };
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
        spp::sparse_hash_map<string, path> name_to_path;
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
    cout << "-I / --insertions         : FASTA file of insertion variant sequences." << endl;
    cout << "version 0.1" << endl;
}

int main(int argc, char** argv){

    char* ref_file = nullptr;
    char* var_file = nullptr;
    char* bed_file = nullptr;
    char* insertion_file = nullptr;
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

    if (ref_file == nullptr){
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
    if (insertion_file != nullptr){
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
    if (bed_file != nullptr){
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
    else{
        for (auto x : tf.seq_to_entry){
            acceptable_chroms.insert(std::string(x.first));
        }
    }

    // Step 0: initialize a graph to hold contig information
    svaha::graph sg;
    for (auto& x : tf.seq_to_entry){
        svaha::pre_contig p;
        std::string name(x.first);
        if (acceptable_chroms.find(name) != acceptable_chroms.end()){
            sg.name_to_contig[name] = p;
            //std::vector<TVCF::variant> bps;
            //sg.name_to_variants[name] = bps;
        }
    }


    // Step 1: add variants from the breakpoint file (i.e. a VCF)
    std::uint64_t var_count = 0;

    if (var_file != nullptr){
        std::ifstream vfi(var_file);
        if (vfi.good()){
            cerr << "Loading variants from " << var_file << endl;

            char* current_contig = nullptr;
            //std::vector<std::uint64_t> breakpoints;
            int mxlsz = 2000;
            char* line = new char[mxlsz];

            while (vfi.getline(line, mxlsz)){
                if (line[0] != '#'){
                    TVCF::variant* var = new TVCF::variant();
                    TVCF::parse_variant_line(line, var);

                    //if (current_contig != NULL && current_contig != var->chrom){
                    // We've moved on to a new contig;
                    // we need to write / flush the local copy of breakpoints
                    //sg.name_to_contig.at(string(current_contig)).breakpoints = breakpoints;
                    //breakpoints.clear();
                    //}
                    // Add the breakpoint triggered by POS
                    //cerr << "Creating breakpoint at [" << var->zero_based_position() <<
                    //    ", " << var->get_reference_end(0) << "] for variant type " << var->get_info("SVTYPE") << endl;

                    if (TFA::hasSequence(tf, var->chrom)){
                        if (acceptable_chroms.find(string(var->chrom)) != acceptable_chroms.end()){
                            std::uint64_t on_chrom_position = var->zero_based_position() + 1;
                            string svtype = var->get_sv_type();
                            TVCF::minimal_allele_t* var_allele = new TVCF::minimal_allele_t();
                            pliib::strcopy(var->chrom, var_allele->chrom);
                            pliib::strcopy(var->get_sv_type().c_str(), var_allele->type);
                            std::uint64_t svend = var->get_sv_end() + 1;
                            if (svend == 0 || on_chrom_position == 0){
                                pliib::strdelete(var_allele->chrom);
                                pliib::strdelete(var_allele->type);
                                delete var;
                                continue;
                            }

                            // Correct for flat alleles
                            if (flat && (svtype == "DEL" || svtype == "INS")){
                                on_chrom_position = on_chrom_position - 1;
                                sg.name_to_contig.at(string(var->chrom)).breakpoints.push_back(on_chrom_position + 1);
                            }
                            sg.name_to_contig.at(string(var->chrom)).breakpoints.push_back(on_chrom_position);
                            var_allele->pos = on_chrom_position;
                            // TODO: going to have to refactor this to handle SNVs
                                
                            if (svtype == "INS"){
                                svend = on_chrom_position + 1;
                            }
                            // Zero-orient that breakpoint
                            var_allele->end = svend -1 ;

                            string second_chr = var->get_info("CHR2");
                            if (second_chr == ""){
                                second_chr = string(var_allele->chrom);
                            }
                            pliib::strcopy(second_chr.c_str(), var_allele->chrom_2);
                            sg.name_to_contig.at(string(var_allele->chrom_2)).breakpoints.push_back(var_allele->end);


                            //TVCF::variant v(*var);
                            //sg.name_to_variants.at(string(var->chrom)).push_back(v);
                            //sg.name_to_contig.at(string(var->chrom)).bp_to_variant[on_chrom_position] = var;
                            //if (svtype != "TRA"){
                            // else if (do_translocations && svtype == "TRA"){
                            //     std::string c2 = var->get_info("CHR2");
                            //     if (acceptable_chroms.find(c2) != acceptable_chroms.end()){
                            //         sg.name_to_contig.at(c2).breakpoints.push_back(var->get_reference_end(0));
                            //     }
                            // }
                            sg.name_to_contig.at(string(var_allele->chrom)).bp_to_allele[var_allele->pos] = var_allele;
                            // if (var_allele->chrom != var_allele->chrom_2)
                            //     sg.name_to_contig.at(string(var_allele->chrom_2)).bp_to_allele[var_allele->end] = var_allele;

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
                    delete var;
                }
            }

            delete [] line;
            //delete var;
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

//    for (auto c : sg.name_to_contig){
//        for (auto v : c.second.breakpoints){
//            cerr << v << endl;
//        }
//    }


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
        c.second.breakpoints.push_back(c.second.seqlen);
    }
    cerr << "Sorted and uniquified breakpoints" << endl;

    std::uint64_t greatest_prev_id = 0;

    for (auto& c : sg.name_to_contig){

	    cerr << "Building reference backbone for " << c.first << "." << endl;

        // Get reference sequence
        if (! TFA::hasSequence(tf, c.first.c_str())){
            cerr << "No sequence found for " << c.first << endl;
            exit(1);
        }
        TFA::getSequence(tf, c.first.c_str(), c.second.seq);
        std::size_t numbp = c.second.breakpoints.size();

        //std::vector<TVCF::variant> contig_vars = sg.name_to_variants.at(c.first);
        std::vector<std::uint64_t> bps = c.second.breakpoints;
        c.second.bp_to_node = new svaha::node* [c.second.seqlen];
        bool insertion_state = false;

        std::uint64_t pos = 0;
        std::uint64_t last_ref_node_pos = 0;
        svaha::node* prev_ref_node = nullptr;

        //std::uint64_t total_seq = 0;
        

        for (size_t i = 0; i < bps.size() - 1; ++i){


            std::uint64_t brk = bps[i];
            // if (brk == 0){
            //     continue;
            // }
            svaha::node* n = sg.create_node();
            pliib::strcopy(c.first.c_str(), n->contig);
            //cerr << c.second.seqlen << " " << pos << " " << brk - pos << endl;
            pliib::strcopy(c.second.seq + pos,  brk - pos, n->seq);
            n->seqlen = brk - pos;

            // Emit the node, caching it if we need it later for a variant.
            cout << n->emit() << endl;
            //pliib::strdelete(n->seq);
            //cout << n->id << " " << n->contig  << " " << n->seqlen << endl;

            TVCF::minimal_allele_t* allele;
            string vtype;
            std::uint64_t vpos;
            std::uint64_t vend;
            // bool position_is_variant = c.second.bp_to_variant.find(brk) != c.second.bp_to_variant.end() ||
            //         (c.second.interchrom_variants.find(brk) != c.second.interchrom_variants.end() && do_translocations);
            bool position_is_variant = c.second.bp_to_allele.find(brk) != c.second.bp_to_allele.end() ||
                    (c.second.bp_to_interchrom.find(brk) != c.second.bp_to_interchrom.end() && do_translocations);

            if (position_is_variant){
                //cerr << prev_ref_node->id << "->" << n->id << endl;
                allele = c.second.bp_to_allele.at(brk);
                vtype = string(allele->type);
                cerr << allele->to_string() << endl;
                // The below line is bad and idk why :'(
                //if (pos != allele->pos - 1) continue;
                
                // Right now, just check if a node exists at our insertion position,
                // but we should instead use a vector at each position and check the variant hashes.
                if (vtype == "DEL" && flat && c.second.bp_to_inserted_node[allele->pos] == nullptr){
                    svaha::node* ins_node = sg.create_node();
                    c.second.bp_to_inserted_node[allele->pos] = ins_node;
                    cerr << "Created ins node at " <<
                    allele->pos << " : id=" << 
                    ins_node->id << 
                    " for DEL." << endl; 
                    // ins_node->seq = n->seq;
                    // ins_node->seqlen = strlen(n->seq);
                    //cout << ins_node->emit() << endl;
                }
                else if (vtype == "INS"){
                    svaha::node* ins_node = sg.create_node();
                    c.second.bp_to_inserted_node[allele->pos] = ins_node;
                }
                else if (vtype == "INV" && flat && c.second.bp_to_inserted_node[allele->pos] == nullptr){
                    svaha::node* ins_node = sg.create_node();
                    c.second.bp_to_inserted_node[allele->pos] = ins_node;
                    cerr << "Created ins node at " << allele->pos <<
                     " : id=" << ins_node->id <<
                     "for INV." << endl; 
                }
                
                //cerr << "N: " << n->id << " brk-1: " << brk - 1 << " pos: " << pos << endl;
                // c.second.bp_to_node[allele->pos] = n;
                // c.second.bp_to_node[allele->pos - 1] = prev_ref_node;
                // c.second.bp_to_node[allele->end] = n;
            }
            else{
                //cerr << "N: " << n->id << " brk-1: " << brk - 1 << " pos: " << pos << endl;
                // c.second.bp_to_node[brk - 1] = n;
                // c.second.bp_to_node[pos+1] = n;
            }
            cerr << "N: " << n->id << " brk-1: " << brk - 1 << " pos: " << pos << endl;
            c.second.bp_to_node[brk - 1] = n;
            c.second.bp_to_node[pos] = n;

            if (pos > 0 && 
		        prev_ref_node != nullptr &&
                n->contig != NULL && 
                prev_ref_node->contig != NULL ){
                svaha::edge e(prev_ref_node, n);
                cout << e.emit() << endl;
                //delete prev_ref_node;
                prev_ref_node = n;
            }
            else if(pos == 0){
                prev_ref_node = n;
            }
            //  if (pos != 0){
            //      cout << c.second.bp_to_node[pos]->id << endl;
            //      cout << c.second.bp_to_node[last_ref_node_pos]->id << endl;
            //  }
            //last_ref_node_pos = pos;

            pos = brk;
        }
        bps.clear();
        c.second.breakpoints.clear();
        c.second.clear_seq();

    }

    spp::sparse_hash_map<std::string, std::uint8_t> processed_alleles;
    // Handle all variants
    for (auto& c : sg.name_to_contig){
        for (auto& pv : c.second.bp_to_allele){

            

            TVCF::minimal_allele_t* allele = pv.second;
            // Alleles are stored on both their pos and end,
            // so we need to check to see if we've processed them so
            // we don't do it twice.
            if (processed_alleles.find(allele->make_id()) != processed_alleles.end()){
                continue;
            }
            // Get the start and end of the variant
            string vtype = string(allele->type);
            std::uint64_t zero_based_vpos = allele->pos - 1;
            std::uint64_t zero_based_vend = allele->end - 1;

            cerr << "Variant allele: " << allele->to_string() << " id: " << allele->make_id() << endl;

            // Get the nodes right before, right after, and right at the start of the variant.
            // For flat variants, S holds the single-base ref sequence.
            svaha::node* s = c.second.bp_to_node[zero_based_vpos + 1];
            svaha::node* last_var_node = c.second.bp_to_node[zero_based_vend];
            svaha::node* pre = c.second.bp_to_node[zero_based_vpos];
            svaha::node* post = c.second.bp_to_node[zero_based_vend + 1];
            
            //#define DEBUG
            #ifdef DEBUG
            cerr << "Processing allele at " << zero_based_vpos << " TO " << zero_based_vend << endl;
                cerr << "Pre: " << pre->id <<
                " variant_node: " << s->id <<
                " Post: " <<  post->id << endl;
            #endif
            // Rules for each variant:
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

            
            if (vtype == "DEL"){    
                if (!flat){
                    // Generate an edge from the node before the variant to the node after the variant.
                    svaha::edge e(pre, post);
                    cout << e.emit() << endl;
                }
                else{
                    // Retrieve the inserted flat-allele alt node, which doesn't yet have a sequence.
                    svaha::node* insy = c.second.bp_to_inserted_node[zero_based_vpos + 1];
                    #ifdef DEBUG
                    cerr << "node " << insy->id << " composes our ref path." << endl;
                    cerr << insy->id << "\'s sequence will be: " << s->seq << endl;
                    #endif
                    // Fill insy's sequence with the sequence of node s.
                    pliib::strcopy(s->seq, insy->seq);
                    insy->seqlen = s->seqlen;

                    // Write the ref node out
                    cout << insy->emit() << endl;
                    // Link that node to our pre ref node
                    svaha::edge e(pre, insy);
                    cout << e.emit() << endl;
                    // Link it to our post ref node
                    svaha::edge f(insy, post);
                    cout << f.emit() << endl;
                }
            }
            else if (vtype == "INS" || vtype == ""){


                if (flat){
                    // Retrieve the inserted inverted node
                    svaha::node* insy = c.second.bp_to_inserted_node[zero_based_vpos + 1];

                }
            }
            else if (vtype == "INV"){
                #ifdef DEBUG
                cerr << "Last var node: " << last_var_node->id << endl;
                #endif
                
                if (!flat){
                    // Generate a tail->tail edge from pre to the last node of the variant run.
                    svaha::edge e(pre, last_var_node, true, false);
                    cout << e.emit() << endl;

                    // Generate a head->head edge from the first variant node to the post node.
                    svaha::edge f(s, post, false, true);
                    cout << f.emit() << endl;
                }
                else{
                    // Grab our inserted inversion node
                    svaha::node* invy = c.second.bp_to_inserted_node[zero_based_vpos + 1];
                    // Fill node sequence with the reverse of the ref node.

                }

            }
            else if (vtype == "TRA"){
                // pre holds the relevant node on this chromosome/contig
                // Get the node on the other contig
                svaha::node* other_contig_node = sg.name_to_contig.at(string(allele->chrom_2)).bp_to_node[zero_based_vend + 1];
                svaha::edge e(pre, other_contig_node);
                cout << e.emit() << endl;

                
            }

            processed_alleles[allele->make_id()] = 1;
        }
    }

    pliib::strdelete(ref_file);
    pliib::strdelete(var_file);
    pliib::strdelete(bed_file);
    pliib::strdelete(insertion_file);

    return 0;
}
