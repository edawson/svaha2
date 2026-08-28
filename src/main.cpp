#include "svaha.hpp"
#include <getopt.h>
#include <fstream>
#include <sstream>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <future>
#include <thread>
#include <locale>
#include <iomanip>
#include <functional>

#include <atomic>

using namespace svaha;

/**
 * @brief Phase 0: Parse variants from VCF or MAF and identify breakpoints
 */
void phase_0_variant_parsing(const std::string& var_path, 
                            const std::string& ref_path,
                            std::map<std::string, BreakpointVector>& chrom_to_bps,
                            std::vector<Variant>& variants,
                            const std::optional<Region>& region = std::nullopt,
                            const std::string& explicit_format = "",
                            bool relax_chrom = false) {
    
    std::cerr << "Phase 0: Parsing variants from " << var_path << "..." << std::endl;
    size_t start_count = variants.size();
    FastaIndex fai(ref_path + ".fai");
    auto normalize_chrom = [&](const std::string& c) -> std::string {
        if (fai.get(c)) return c;
        if (relax_chrom) {
            if (c.size() > 3 && c.substr(0, 3) == "chr") {
                std::string no_chr = c.substr(3);
                if (fai.get(no_chr)) return no_chr;
            }
            std::string with_chr = "chr" + c;
            if (fai.get(with_chr)) return with_chr;
        }
        return c;
    };

    auto process_variant = [&](const Variant& v_in) {
        Variant v = v_in;
        std::string orig_chrom = v.chrom;
        v.chrom = normalize_chrom(v.chrom);
        if (!v.chrom_2.empty()) v.chrom_2 = normalize_chrom(v.chrom_2);

        if (region && !region->overlaps(v.chrom, v.pos, v.end)) return;
        
        if (!fai.get(v.chrom)) {
            std::cerr << "Error: Chromosome '" << orig_chrom << "' not found in reference." << std::endl;
            if (!relax_chrom) {
                std::cerr << "Use --relax-chrom to allow matching between 'chr1' and '1'." << std::endl;
                exit(1);
            }
            return;
        }
        
        variants.push_back(v);
        chrom_to_bps[v.chrom].add(v.pos);
        chrom_to_bps[v.chrom].add(v.end);
        if (!v.chrom_2.empty() && v.chrom != v.chrom_2) {
            if (fai.get(v.chrom_2)) {
                chrom_to_bps[v.chrom_2].add(v.pos_2);
                chrom_to_bps[v.chrom_2].add(v.pos_2 + 1); // Anchor for BND
            } else if (!relax_chrom) {
                std::cerr << "Error: Chromosome '" << v.chrom_2 << "' not found in reference." << std::endl;
                exit(1);
            }
        }
    };

    std::string format = explicit_format;
    if (format.empty()) {
        // Detect format
        std::ifstream test_ifs(var_path);
        std::string first_line;
        std::getline(test_ifs, first_line);
        test_ifs.close();

        if (first_line.find("Sample_Id") != std::string::npos && (first_line.find("Connection_Type") != std::string::npos || first_line.find("SV_Status") != std::string::npos)) {
            format = "sv";
        } else if (first_line.find("Hugo_Symbol") != std::string::npos || first_line.find("#genome_nexus") != std::string::npos) {
            format = "cbio_maf";
        } else if (var_path.length() >= 4 && var_path.substr(var_path.find_last_of(".") + 1) == "maf" && !first_line.empty() && first_line[0] == '#') {
            format = "maf";
        } else {
            format = "vcf";
        }
    }

    if (format == "sv" || format == "bedpe") {
        SvParser::parse(var_path, process_variant);
    } else if (format == "cbio_maf") {
        CbioMafParser::parse(var_path, process_variant);
    } else if (format == "maf") {
        MafParser::parse(var_path, [&](const MafEntry& ref, const MafEntry& alt) {
            std::string clean_alt = "";
            for (char c : alt.text) if (c != '-') clean_alt += c;
            Variant v = {ref.src, ref.start, ref.start + ref.size, "", 
                        PackedSequence::pack(clean_alt), (uint32_t)clean_alt.length(), 
                        "MAF_BLOCK", ref.src, ref.start + ref.size};
            process_variant(v);
        });
    } else {
        // Assume VCF
        std::ifstream vcf_file(var_path);
        std::string line;
        while (std::getline(vcf_file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::stringstream ss(line);
            std::string chrom, pos_str, id, ref, alts, qual, filter, info;
            std::getline(ss, chrom, '\t');
            std::getline(ss, pos_str, '\t');
            
            uint64_t pos = std::stoull(pos_str) - 1; // 0-based
            std::getline(ss, id, '\t');
            std::getline(ss, ref, '\t');
            std::getline(ss, alts, '\t');
            std::getline(ss, qual, '\t');
            std::getline(ss, filter, '\t');
            std::getline(ss, info, '\t');

            std::stringstream ss_alts(alts);
            std::string alt;
            while (std::getline(ss_alts, alt, ',')) {
                uint64_t end = pos + ref.length();
                std::string type = "";
                std::string chrom_2 = chrom;
                uint64_t pos_2 = end;

                if (info.find("SVTYPE=DEL") != std::string::npos) type = "DEL";
                else if (info.find("SVTYPE=INS") != std::string::npos) type = "INS";
                else if (info.find("SVTYPE=INV") != std::string::npos) type = "INV";
                else if (info.find("SVTYPE=BND") != std::string::npos || info.find("SVTYPE=TRA") != std::string::npos) type = "BND";
                else if (ref.length() == 1 && alt.length() == 1) type = "SNP";
                else type = "INDEL";

                if (info.find("END=") != std::string::npos) {
                    size_t start_idx = info.find("END=") + 4;
                    size_t end_idx = info.find(';', start_idx);
                    end = std::stoull(info.substr(start_idx, end_idx - start_idx));
                }

                if (type == "BND") {
                    if (info.find("CHR2=") != std::string::npos) {
                        size_t start_idx = info.find("CHR2=") + 5;
                        size_t end_idx = info.find(';', start_idx);
                        chrom_2 = info.substr(start_idx, end_idx - start_idx);
                    }
                    if (info.find("POS2=") != std::string::npos) {
                        size_t start_idx = info.find("POS2=") + 5;
                        size_t end_idx = info.find(';', start_idx);
                        pos_2 = std::stoull(info.substr(start_idx, end_idx - start_idx));
                    }
                }

                // Handle anchor bases for small variants
                uint64_t v_pos = pos;
                uint64_t v_end = end;
                std::string v_alt = alt;
                if ((type == "SNP" || type == "INDEL") && !ref.empty() && !alt.empty()) {
                    if (ref[0] == alt[0]) {
                        v_pos++;
                        v_alt = alt.substr(1);
                        if (ref.length() == 1) v_end = v_pos; // Insertion
                        else v_end = v_pos + ref.length() - 1;
                    }
                }

                Variant v = {chrom, v_pos, v_end, ref, PackedSequence::pack(v_alt), (uint32_t)v_alt.length(), type, chrom_2, pos_2};
                process_variant(v);
            }
        }
    }
    std::cerr << "  Done. Loaded " << variants.size() - start_count << " variants from this file." << std::endl;
}

/**
 * @brief Phase 1: Build Topology (Nodes and Edges)
 */
void phase_1_topology(const std::map<std::string, BreakpointVector>& chrom_to_bps,
                      const std::vector<Variant>& variants,
                      std::vector<Node>& nodes,
                      std::vector<Edge>& edges,
                      std::map<std::pair<std::string, uint64_t>, uint64_t>& pos_to_node_id,
                      std::map<std::pair<std::string, uint64_t>, uint64_t>& end_to_node_id) {
    std::cerr << "Phase 1: Building graph topology..." << std::endl;
    uint64_t node_id = 1;

    for (auto const& [chrom, bps] : chrom_to_bps) {
        for (size_t i = 0; i < bps.size() - 1; ++i) {
            Node n = {node_id++, chrom, bps[i], 0, (uint32_t)(bps[i+1] - bps[i])};
            nodes.push_back(n);
            pos_to_node_id[{chrom, bps[i]}] = n.id;
            end_to_node_id[{chrom, bps[i+1]}] = n.id;
            if (i > 0) edges.push_back({nodes[nodes.size()-2].id, n.id, true, true});
        }
    }
    size_t backbone_nodes = nodes.size();

    for (auto const& v : variants) {
        // ... (existing variant node creation)
        if (v.type == "SNP" || v.type == "INDEL" || v.type == "INS" || v.type == "MAF_BLOCK" || v.type == "DNP" || v.type == "TNP" || v.type == "ONP") {
            Node alt = {node_id++, v.chrom, v.pos, 1, v.alt_len};
            nodes.push_back(alt);
            uint64_t pre_id = end_to_node_id[{v.chrom, v.pos}];
            uint64_t post_id = pos_to_node_id[{v.chrom, v.end}];
            if (pre_id) edges.push_back({pre_id, alt.id, true, true});
            if (post_id) edges.push_back({alt.id, post_id, true, true});
        } else if (v.type == "DEL" || v.type == "DELETION") {
            uint64_t pre_id = end_to_node_id[{v.chrom, v.pos}];
            uint64_t post_id = pos_to_node_id[{v.chrom, v.end}];
            if (pre_id && post_id) edges.push_back({pre_id, post_id, true, true});
        } else if (v.type == "INV" || v.type == "INVERSION") {
            uint64_t pre_id = end_to_node_id[{v.chrom, v.pos}];
            uint64_t post_id = pos_to_node_id[{v.chrom, v.end}];
            uint64_t start_node_id = pos_to_node_id[{v.chrom, v.pos}];
            uint64_t end_node_id = end_to_node_id[{v.chrom, v.end}];
            
            // If from cBio SV parser, we might have specific orientations
            if (!v.from_forward || !v.to_forward) {
                if (pre_id && post_id) edges.push_back({pre_id, post_id, v.from_forward, v.to_forward});
            } else {
                // Standard VCF inversion logic
                if (pre_id && end_node_id) edges.push_back({pre_id, end_node_id, true, false});
                if (start_node_id && post_id) edges.push_back({start_node_id, post_id, false, true});
            }
        } else if (v.type == "BND" || v.type == "TRANSLOCATION" || v.type == "FUSION") {
            uint64_t pre_id = end_to_node_id[{v.chrom, v.pos}];
            uint64_t post_id = pos_to_node_id[{v.chrom_2, v.pos_2}];
            if (pre_id && post_id) edges.push_back({pre_id, post_id, v.from_forward, v.to_forward});
        } else if (v.type == "DUPLICATION") {
            uint64_t start_id = pos_to_node_id[{v.chrom, v.pos}];
            uint64_t end_id = end_to_node_id[{v.chrom, v.end}];
            if (start_id && end_id) edges.push_back({end_id, start_id, true, true});
        }
    }
    std::cerr << "  Done. Created " << backbone_nodes << " backbone nodes and " << nodes.size() - backbone_nodes << " alternate nodes." << std::endl;
    std::cerr << "  Total edges: " << edges.size() << std::endl;
}

/**
 * @brief Helper to fetch sequence from mmap'ed ref with FAI
 */
std::string fetch_ref_seq(const MemoryMappedFile& ref_mmap, const FastaIndex& fai, 
                         const std::string& chrom, uint64_t start, uint64_t length) {
    auto entry = fai.get(chrom);
    if (!entry) return std::string(length, 'N');
    std::string seq;
    seq.reserve(length);
    uint64_t remaining = length;
    uint64_t cur_pos = start;
    while (remaining > 0) {
        uint64_t line_offset = entry->offset + (cur_pos / entry->line_bases) * entry->line_width + (cur_pos % entry->line_bases);
        uint64_t bases_in_line = entry->line_bases - (cur_pos % entry->line_bases);
        uint64_t to_read = std::min(remaining, bases_in_line);
        seq += ref_mmap.get_substring(line_offset, to_read);
        remaining -= to_read;
        cur_pos += to_read;
    }
    return seq;
}

/**
 * @brief Phase 2: Materialize GFA with Lazy Sequence Loading
 */
void phase_2_materialize(const std::string& ref_path,
                         const std::vector<Node>& nodes,
                         const std::vector<Edge>& edges,
                         const std::vector<Variant>& variants) {
    std::cerr << "Phase 2: Materializing GFA..." << std::endl;
    MemoryMappedFile ref_mmap(ref_path);
    FastaIndex fai(ref_path + ".fai");
    
    // Index alt alleles for fast lookup
    std::map<std::pair<std::string, uint64_t>, const Variant*> alt_map;
    for (const auto& v : variants) {
        alt_map[{v.chrom, v.pos}] = &v;
    }

    std::cout << "H\tVN:Z:1.1" << std::endl;
    size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
    size_t chunk_size = (nodes.size() + num_threads - 1) / num_threads;
    std::vector<std::future<std::vector<std::string>>> futures;
    std::atomic<size_t> nodes_processed{0};
    
    for (size_t t = 0; t < num_threads; ++t) {
        size_t start_idx = t * chunk_size;
        size_t end_idx = std::min(start_idx + chunk_size, nodes.size());
        if (start_idx >= end_idx) break;
        futures.push_back(std::async(std::launch::async, [&, start_idx, end_idx]() {
            std::vector<std::string> chunk_seqs;
            for (size_t i = start_idx; i < end_idx; ++i) {
                const auto& n = nodes[i];
                std::string seq = "*";
                if (n.SR == 0) {
                    if (ref_mmap.is_open()) seq = fetch_ref_seq(ref_mmap, fai, n.SN, n.SO, n.length);
                } else {
                    auto it = alt_map.find({n.SN, n.SO});
                    if (it != alt_map.end()) {
                        seq = PackedSequence::unpack(it->second->packed_alt, it->second->alt_len);
                    }
                }
                chunk_seqs.push_back(seq);
                
                size_t processed = ++nodes_processed;
                if (processed % 1000000 == 0) {
                    std::cerr << "  Processed " << processed << " / " << nodes.size() << " nodes..." << std::endl;
                }
            }
            return chunk_seqs;
        }));
    }
    size_t node_idx = 0;
    for (auto& f : futures) {
        for (const auto& seq : f.get()) std::cout << nodes[node_idx++].to_gfa_s(seq) << std::endl;
    }
    std::cerr << "  Writing edges..." << std::endl;
    for (const auto& e : edges) std::cout << e.to_gfa_l() << std::endl;
    std::cerr << "Phase 2 Complete." << std::endl;
}

void print_usage() {
    std::cout << "svaha: high-performance variation graph construction\n"
              << "Usage: svaha <command> [options]\n\n"
              << "Commands:\n"
              << "  build    Construct a GFA from reference and variants\n"
              << "           Options:\n"
              << "             -r <file>  Reference FASTA (required)\n"
              << "             -v <file>  VCF variants\n"
              << "             --maf <file>  Mutation Annotation Format (cBioPortal)\n"
              << "             --sv <file>   Structural Variants (cBioPortal data_sv.txt)\n"
              << "             --bedpe <file> BEDPE structural variants\n"
              << "             -m <int>   Max node size (default 32)\n"
              << "             -R <reg>   Region (chrom:start-end)\n"
              << "             --relax-chrom Allow 'chr1' to match '1' and vice-versa\n"
              << "  stats    Show statistics for a GFA file\n"
              << "  view     Visualize a portion of the graph\n";
}

int main(int argc, char** argv) {
    if (argc < 2) { print_usage(); return 1; }
    std::string command = argv[1];
    if (command == "build") {
        std::vector<std::pair<std::string, std::string>> var_inputs;
        std::string ref_path;
        std::optional<Region> region;
        int max_node_size = 32;
        bool relax_chrom = false;
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if ((arg == "-v" || arg == "--vcf") && i + 1 < argc) var_inputs.push_back({argv[++i], "vcf"});
            else if (arg == "--maf" && i + 1 < argc) var_inputs.push_back({argv[++i], "cbio_maf"});
            else if (arg == "--sv" && i + 1 < argc) var_inputs.push_back({argv[++i], "sv"});
            else if (arg == "--bedpe" && i + 1 < argc) var_inputs.push_back({argv[++i], "bedpe"});
            else if (arg == "-r" && i + 1 < argc) ref_path = argv[++i];
            else if (arg == "-m" && i + 1 < argc) max_node_size = std::stoi(argv[++i]);
            else if (arg == "-R" && i + 1 < argc) region = Region::parse(argv[++i]);
            else if (arg == "--relax-chrom") relax_chrom = true;
        }
        if (var_inputs.empty() || ref_path.empty()) {
            std::cerr << "Error: Variant input and -r (Reference) are required for build." << std::endl;
            return 1;
        }

        std::map<std::string, BreakpointVector> chrom_to_bps;
        std::vector<Variant> variants;
        for (const auto& input : var_inputs) {
            phase_0_variant_parsing(input.first, ref_path, chrom_to_bps, variants, region, input.second, relax_chrom);
        }

        // Global breakpoints
        std::cerr << "Finalizing breakpoints and discretizing reference..." << std::endl;
        FastaIndex fai(ref_path + ".fai");
        std::vector<std::string> target_chroms;
        if (region) target_chroms.push_back(region->chrom);
        else for (auto const& [c, bps] : chrom_to_bps) target_chroms.push_back(c);

        for (const auto& c : target_chroms) {
            auto entry = fai.get(c);
            uint64_t c_start = region ? region->start : 0;
            uint64_t c_end = region ? region->end : (entry ? entry->length : 0);
            chrom_to_bps[c].add(c_start);
            chrom_to_bps[c].add(c_end);
            for (uint64_t b = c_start + max_node_size; b < c_end; b += max_node_size) chrom_to_bps[c].add(b);
        }

        for (auto& entry : chrom_to_bps) entry.second.finalize();

        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::map<std::pair<std::string, uint64_t>, uint64_t> pos_to_node_id, end_to_node_id;
        phase_1_topology(chrom_to_bps, variants, nodes, edges, pos_to_node_id, end_to_node_id);
        phase_2_materialize(ref_path, nodes, edges, variants);
    } else if (command == "stats") {
        std::string gfa_path;
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-h" || arg == "--help") { std::cout << "Usage: svaha stats <file.gfa>\n"; return 0; }
            gfa_path = arg;
        }
        if (gfa_path.empty()) { std::cerr << "Error: GFA file required for stats." << std::endl; return 1; }
        std::ifstream gfa_file(gfa_path);
        if (!gfa_file.is_open()) { std::cerr << "Error: Could not open GFA file: " << gfa_path << std::endl; return 1; }
        uint64_t node_count = 0, edge_count = 0, total_len = 0;
        std::string line;
        while (std::getline(gfa_file, line)) {
            if (line.empty()) continue;
            if (line[0] == 'S') {
                node_count++;
                std::stringstream ss(line);
                std::string type, id, seq;
                ss >> type >> id >> seq;
                if (seq != "*") total_len += seq.length();
            } else if (line[0] == 'L') edge_count++;
        }
        std::cout << "GFA Statistics for: " << gfa_path << "\nNodes: " << node_count << "\nEdges: " << edge_count << "\nTotal Sequence Length: " << total_len << " bp" << std::endl;
    } else if (command == "view") {
        std::cout << "View command (Implementation pending...)" << std::endl;
    } else { print_usage(); return 1; }
    return 0;
}
