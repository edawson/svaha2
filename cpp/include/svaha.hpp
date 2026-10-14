#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <map>
#include <algorithm>
#include <memory>
#include <sstream>
#include <optional>
#include <fstream>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <mutex>
#include <condition_variable>
#include <locale>
#include <iomanip>
#include <functional>

namespace svaha {

/**
 * @brief Memory-mapped file wrapper for lazy sequence loading
 */
class MemoryMappedFile {
public:
    MemoryMappedFile(const std::string& path) {
        fd = open(path.c_str(), O_RDONLY);
        if (fd == -1) return;

        struct stat sb;
        if (fstat(fd, &sb) == -1) {
            close(fd);
            fd = -1;
            return;
        }
        file_size = sb.st_size;

        data = static_cast<char*>(mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            close(fd);
            fd = -1;
            data = nullptr;
        }
    }

    ~MemoryMappedFile() {
        if (data) munmap(data, file_size);
        if (fd != -1) close(fd);
    }

    bool is_open() const { return data != nullptr; }
    const char* get_data() const { return data; }
    size_t size() const { return file_size; }

    std::string get_substring(size_t start, size_t length) const {
        if (!data || start >= file_size) return "";
        size_t actual_length = std::min(length, file_size - start);
        return std::string(data + start, actual_length);
    }

private:
    int fd = -1;
    char* data = nullptr;
    size_t file_size = 0;
};

/**
 * @brief Simple FASTA index (.fai) parser
 */
struct FaiEntry {
    std::string name;
    uint64_t length;
    uint64_t offset;
    uint32_t line_bases;
    uint32_t line_width;
};

class FastaIndex {
public:
    FastaIndex(const std::string& fai_path) {
        std::ifstream ifs(fai_path);
        if (!ifs.is_open()) return;
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            FaiEntry entry;
            std::getline(ss, entry.name, '\t');
            ss >> entry.length >> entry.offset >> entry.line_bases >> entry.line_width;
            entries[entry.name] = entry;
        }
    }

    std::optional<FaiEntry> get(const std::string& name) const {
        auto it = entries.find(name);
        if (it != entries.end()) return it->second;
        return std::nullopt;
    }

private:
    std::map<std::string, FaiEntry> entries;
};

/**
 * @brief 2-bit sequence packer
 */
class PackedSequence {
public:
    static std::vector<uint64_t> pack(const std::string& seq) {
        if (seq == "-" || seq.empty()) return {};
        size_t n_blocks = (seq.length() + 31) / 32;
        std::vector<uint64_t> packed(n_blocks, 0);
        for (size_t i = 0; i < seq.length(); ++i) {
            uint64_t val = 0;
            switch (toupper(seq[i])) {
                case 'A': val = 0; break;
                case 'C': val = 1; break;
                case 'G': val = 2; break;
                case 'T': val = 3; break;
                default: val = 0; // Treat N as A for 2-bit
            }
            packed[i / 32] |= (val << (2 * (i % 32)));
        }
        return packed;
    }

    static std::string unpack(const std::vector<uint64_t>& packed, size_t length) {
        if (length == 0) return "";
        std::string seq(length, ' ');
        for (size_t i = 0; i < length; ++i) {
            uint64_t val = (packed[i / 32] >> (2 * (i % 32))) & 0x3;
            switch (val) {
                case 0: seq[i] = 'A'; break;
                case 1: seq[i] = 'C'; break;
                case 2: seq[i] = 'G'; break;
                case 3: seq[i] = 'T'; break;
            }
        }
        return seq;
    }
};

/**
 * @brief rGFA Node representation with Stable Coordinates
 */
struct Node {
    uint64_t id;
    std::string SN;         // Stable Name (e.g., "chr1")
    uint64_t SO;            // Stable Offset (0-based)
    uint32_t SR;            // Stable Rank (0=ref, 1+=alt)
    uint32_t length;
    
    std::string to_gfa_s(const std::string& sequence = "*") const {
        std::stringstream ss;
        ss << "S\t" << id << "\t" << sequence 
           << "\tSN:Z:" << SN << "\tSO:i:" << SO << "\tSR:i:" << SR;
        return ss.str();
    }
};

/**
 * @brief Edge representation using IDs
 */
struct Edge {
    uint64_t from;
    uint64_t to;
    bool from_forward = true;
    bool to_forward = true;

    std::string to_gfa_l() const {
        std::stringstream ss;
        ss << "L\t" << from << "\t" << (from_forward ? "+" : "-") 
           << "\t" << to << "\t" << (to_forward ? "+" : "-") << "\t0M";
        return ss.str();
    }
};

/**
 * @brief Simple Bit-Vector for breakpoints (Phase 0/1)
 * In a full implementation, we'd use sdsl-lite's bit_vector for rank/select.
 */
class BreakpointVector {
public:
    void add(uint64_t pos) {
        breakpoints.push_back(pos);
    }
    
    void finalize() {
        std::sort(breakpoints.begin(), breakpoints.end());
        auto last = std::unique(breakpoints.begin(), breakpoints.end());
        breakpoints.erase(last, breakpoints.end());
    }

    const std::vector<uint64_t>& get() const { return breakpoints; }
    
    size_t size() const { return breakpoints.size(); }
    uint64_t operator[](size_t i) const { return breakpoints[i]; }

private:
    std::vector<uint64_t> breakpoints;
};

/**
 * @brief Genomic Region representation
 */
struct Region {
    std::string chrom;
    uint64_t start; // 0-based
    uint64_t end;   // 0-based

    static Region parse(const std::string& s) {
        Region r;
        size_t colon = s.find(':');
        if (colon == std::string::npos) {
            r.chrom = s;
            r.start = 0;
            r.end = UINT64_MAX;
        } else {
            r.chrom = s.substr(0, colon);
            size_t dash = s.find('-', colon);
            r.start = std::stoull(s.substr(colon + 1, dash - colon - 1));
            r.end = std::stoull(s.substr(dash + 1));
        }
        return r;
    }

    bool contains(const std::string& c, uint64_t p) const {
        return c == chrom && p >= start && p < end;
    }

    bool overlaps(const std::string& c, uint64_t s, uint64_t e) const {
        return c == chrom && s < end && e > start;
    }
};

/**
 * @brief Variant representation for Phase 0
 */
struct Variant {
    std::string chrom;
    uint64_t pos;           // 0-based start
    uint64_t end;           // 0-based end
    std::string ref;
    std::vector<uint64_t> packed_alt; // 2-bit packed alt sequence
    uint32_t alt_len;
    std::string type;       // SNP, INS, DEL, INV, etc.
    std::string chrom_2;    // For translocations
    uint64_t pos_2;         // For translocations
    bool from_forward = true;
    bool to_forward = true;
};

/**
 * @brief cBioPortal MAF Parser
 */
class CbioMafParser {
public:
    static void parse(const std::string& path, std::function<void(const Variant&)> on_variant) {
        std::ifstream ifs(path);
        if (!ifs.is_open()) return;

        std::string line;
        std::map<std::string, int> header;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string token;
            int col = 0;
            if (header.empty()) {
                while (std::getline(ss, token, '\t')) header[token] = col++;
                continue;
            }

            std::vector<std::string> row;
            while (std::getline(ss, token, '\t')) row.push_back(token);
            if (row.size() < header.size()) continue;

            try {
                Variant v;
                v.chrom = row[header["Chromosome"]];
                v.pos = std::stoull(row[header["Start_Position"]]) - 1;
                v.end = std::stoull(row[header["End_Position"]]);
                v.ref = row[header["Reference_Allele"]];
                std::string alt = row[header["Tumor_Seq_Allele2"]];
                v.type = row[header["Variant_Type"]];
                
                if (alt == "-") alt = "";
                v.packed_alt = PackedSequence::pack(alt);
                v.alt_len = alt.length();
                
                if (v.type == "SNP" && v.end == v.pos) v.end = v.pos + 1;
                on_variant(v);
            } catch (...) { continue; }
        }
    }
};

/**
 * @brief cBioPortal Structural Variant (SV) Parser
 */
class SvParser {
public:
    static void parse(const std::string& path, std::function<void(const Variant&)> on_variant) {
        std::ifstream ifs(path);
        if (!ifs.is_open()) return;

        std::string line;
        std::map<std::string, int> header;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string token;
            int col = 0;
            if (header.empty()) {
                while (std::getline(ss, token, '\t')) {
                    header[token] = col++;
                }
                continue;
            }

            std::vector<std::string> row;
            while (std::getline(ss, token, '\t')) row.push_back(token);
            if (row.size() < header.size()) continue;

            try {
                Variant v;
                v.chrom = row[header["Site1_Chromosome"]];
                v.pos = std::stoull(row[header["Site1_Position"]]) - 1;
                v.chrom_2 = row[header["Site2_Chromosome"]];
                v.pos_2 = std::stoull(row[header["Site2_Position"]]) - 1;
                v.type = row[header["Class"]]; 
                
                std::string conn = row[header["Connection_Type"]];
                if (conn == "3to5") { v.from_forward = true; v.to_forward = true; }
                else if (conn == "5to5") { v.from_forward = false; v.to_forward = true; }
                else if (conn == "3to3") { v.from_forward = true; v.to_forward = false; }
                else if (conn == "5to3") { v.from_forward = false; v.to_forward = false; }

                if (v.chrom == v.chrom_2) {
                    v.end = v.pos_2;
                    if (v.end < v.pos) std::swap(v.pos, v.end);
                } else {
                    v.end = v.pos + 1;
                }
                on_variant(v);
            } catch (...) { continue; }
        }
    }
};

/**
 * @brief MAF alignment block entry
 */
struct MafEntry {
    std::string src;
    uint64_t start;
    uint64_t size;
    char strand;
    uint64_t srcSize;
    std::string text;
};

/**
 * @brief High-performance MAF parser
 */
class MafParser {
public:
    static void parse(const std::string& maf_path, 
                     std::function<void(const MafEntry&, const MafEntry&)> on_block) {
        std::ifstream ifs(maf_path);
        if (!ifs.is_open()) return;

        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#' || line[0] == 'a') continue;
            if (line[0] == 's') {
                MafEntry ref = parse_s_line(line);
                if (std::getline(ifs, line) && line[0] == 's') {
                    MafEntry alt = parse_s_line(line);
                    on_block(ref, alt);
                }
            }
        }
    }

private:
    static MafEntry parse_s_line(const std::string& line) {
        std::stringstream ss(line);
        std::string s, src, start, size, strand, srcSize, text;
        ss >> s >> src >> start >> size >> strand >> srcSize >> text;
        return {src, std::stoull(start), std::stoull(size), strand[0], std::stoull(srcSize), text};
    }
};

} // namespace svaha
