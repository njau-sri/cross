#pragma once

#include <string>
#include <vector>

using allele_t = unsigned char;

struct VcfEntry
{
    std::string chr;
    std::string id;
    std::vector<std::string> as;
    std::vector<allele_t> gt;
    int pos = -1;
    int ploidy = -1;
};

struct Genotype
{
    std::vector<std::string> ind;
    std::vector<std::string> loc;
    std::vector<std::string> chr;
    std::vector<int> pos;
    std::vector< std::vector<allele_t> > dat;
    std::vector< std::vector<std::string> > allele;
    int ploidy = -1;
};

int parse_vcf_header(const std::string &s, std::vector<std::string> &v);

int parse_vcf_entry(const std::string &s, VcfEntry &e);

int read_vcf(const std::string &filename, Genotype &gt);

int write_vcf(const Genotype &gt, const std::string &filename, bool force_diploid = true);
