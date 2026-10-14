#include "svaha.hpp"
#include <cassert>
#include <iostream>

using namespace svaha;

void test_region() {
    Region r = Region::parse("chr1:100-200");
    assert(r.chrom == "chr1");
    assert(r.start == 100);
    assert(r.end == 200);
    assert(r.contains("chr1", 150));
    assert(!r.contains("chr1", 250));
    assert(r.overlaps("chr1", 50, 150));
    assert(!r.overlaps("chr1", 250, 300));
    std::cout << "test_region passed" << std::endl;
}

void test_packed_sequence() {
    std::string seq = "ACGTACGT";
    auto packed = PackedSequence::pack(seq);
    std::string unpacked = PackedSequence::unpack(packed, seq.length());
    assert(seq == unpacked);
    std::cout << "test_packed_sequence passed" << std::endl;
}

void test_breakpoint_vector() {
    BreakpointVector bv;
    bv.add(100);
    bv.add(50);
    bv.add(100);
    bv.finalize();
    assert(bv.size() == 2);
    assert(bv[0] == 50);
    assert(bv[1] == 100);
    std::cout << "test_breakpoint_vector passed" << std::endl;
}

int main() {
    test_region();
    test_packed_sequence();
    test_breakpoint_vector();
    return 0;
}
