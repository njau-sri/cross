#include <memory>
#include <random>
#include <fstream>
#include <iostream>
#include "appcross.h"
#include "cmdline.h"
#include "strsplit.h"
#include "util.h"
#include "hapmapio.h"
#include "vcfio.h"
#include "plinkio.h"
#include "rtmio.h"

namespace {

std::mt19937 RNG;

int read_effect(const string &filename, vector<string> &trait, vector<string> &locus, vector<string> &allele, vector<double> &effect)
{
    std::ifstream ifs(filename);

    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;
    string curr;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if ( vs.empty() )
            continue;

        if ( vs.size() == 1 && vs[0].find('>') == 0 ) {
            curr = vs[0].substr(1);
            continue;
        }

        if (vs.size() != 3) {
            std::cerr << "ERROR: column count doesn't match at line " << ln
                      << " (" << vs.size() << "!=3): " << filename << "\n";
            return 1;
        }

        if ( ! curr.empty() ) {
            trait.push_back(curr);
            locus.push_back(vs[0]);
            allele.push_back(vs[1].substr(0,vs[1].find(':')));
            effect.push_back(number<double>(vs[2]));
        }
    }

    return 0;
}

int read_map(const string &filename, vector<string> &loc, vector<string> &chr, vector<double> &pos)
{
    std::ifstream ifs(filename);

    if ( ! ifs ) {
        std::cerr << "ERROR: can't open file: " << filename << "\n";
        return 1;
    }

    auto delim = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };

    size_t ln = 0;

    for (string line; std::getline(ifs,line); ) {
        ++ln;

        vector<string> vs;
        strsplit(delim, line.begin(), line.end(), vs);
        if ( vs.empty() )
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at line " << ln << ": " << filename << "\n";
            return 1;
        }

        bool ok = false;
        auto x = number<double>(vs[2], &ok);
        if ( ! ok || x < 0.0 )
            return 1;

        loc.push_back(vs[0]);
        chr.push_back(vs[1]);
        pos.push_back(x);
    }

    return 0;
}

bool is_homozygous(const vector<allele_t> &v)
{
    auto n = v.size() / 2;

    for (size_t i = 0; i < n; ++i) {
        if (v[i*2] != v[i*2+1])
            return false;
    }

    return true;
}

double percentile(int p, const vector<double> &v)
{
    if (p <= 1)
        return v.front();

    if (p >= 100)
        return v.back();

    double x = v.size() * p / 100.0;
    auto j = static_cast<size_t>(x);

    if (x - j == 0.0)
        return (v[j-1] + v[j]) / 2;

    return v[j];
}

template<typename T>
double mean(const vector<T> &v)
{
    return std::accumulate(v.begin(), v.end(), double(0)) / v.size();
}

template<typename T>
double sd(const vector<T> &v)
{
    auto m = mean(v);
    double ss = 0;
    for (auto e : v)
        ss += (e-m)*(e-m);
    auto n = v.size();
    return std::sqrt(ss/(n-1));
}

} // namespace

int AppCross::run(int argc, char *argv[])
{
    auto cmd = std::make_shared<CmdLine>("cross [options]");

    cmd->add("--vcf", "VCF file", "");
    cmd->add("--ped", "PLINK PED/MAP file prefix", "");
    cmd->add("--hmp", "HapMap file",  "");
    cmd->add("--geno", "genotype data file", "");
    cmd->add("--map", "linkage map data file", "");
    cmd->add("--effect", "allele effect file", "");
    cmd->add("--pheno", "phenotype data file", "");
    cmd->add("--out", "output file", "appcross.out");
    cmd->add("--type", "cross type", "single");
    cmd->add("--size", "sample size", "2000");
    cmd->add("--perc", "percentiles of sample",  "1,25,50,75,100");

    cmd->add("--indep", "force independent assortment model");

    if (argc < 2) {
        cmd->help();
        return 1;
    }

    cmd->parse(argc, argv);

    m_par.vcf = cmd->get("--vcf");
    m_par.ped = cmd->get("--ped");
    m_par.hmp = cmd->get("--hmp");
    m_par.geno = cmd->get("--geno");
    m_par.map = cmd->get("--map");
    m_par.effect = cmd->get("--effect");
    m_par.pheno = cmd->get("--pheno");
    m_par.out = cmd->get("--out");
    m_par.type = cmd->get("--type");
    m_par.perc = cmd->get("--perc");

    m_par.size = number<int>(cmd->get("--size"));
    if (m_par.size < 1) {
        std::cerr << "WARNING: invalid argument: --size " << cmd->get("--size") << "\n";
        m_par.size = 2000;
    }

    m_par.indep = cmd->has("--indep");

    cmd.reset();

    std::transform(m_par.type.begin(), m_par.type.end(), m_par.type.begin(), ::toupper);

    int info = perform();

    return info;
}

int AppCross::perform()
{
    parse_percentile();

    load_phenotype();

    load_genotype();

    load_linakge_map();

    load_allele_effect();

    save_allele_matrix();

    bool linkage = check_linkage_model();

    if ( linkage && m_par.type != "SINGLE") {
        linkage = false;
        std::cerr << "WARNING: linkage model is not available for " << m_par.type << " cross\n";
    }

    if ( linkage )
        make_genetic_map();

    predict_parent();

    write_phenotype(m_pred, m_par.out + ".parent.txt");

    std::cerr << "INFO: recombination simulation model: " << (linkage ? "linkage" : "independence") << "\n";

    std::cerr << "INFO: performing " << m_par.type << " cross prediction...\n";

    if (m_par.type == "SINGLE") {
        if ( linkage )
            cross_single_linkage();
        else
            cross_single_indep();
    }
    else if (m_par.type == "THREEWAY") {
        if ( linkage )
            cross_threeway_linkage();
        else
            cross_threeway_indep();
    }
    else if (m_par.type == "DOUBLE") {
        if ( linkage )
            cross_double_linkage();
        else
            cross_double_indep();
    }
    else {
        std::cerr << "ERROR: invalid type of cross: " << m_par.type << "\n";
        return 1;
    }

    std::cerr << "INFO: CROSS completed successfully\n";

    return 1;
}

void AppCross::parse_percentile()
{
    m_perc.clear();

    auto delim = [](char c) { return c == ',' || c == ' ' || c == '\t' || c == '\r'; };

    vector<string> vs;
    strsplit(delim, m_par.perc.begin(), m_par.perc.end(), vs);
    for (auto &e : vs) {
        auto x = number<int>(e);
        if (x > 0 && x <= 100)
            m_perc.push_back(x);
        else
            std::cerr << "WARNING: invalid percent value: " << e << "\n";
    }

    if ( m_perc.empty() )
        m_perc = {1, 25, 50, 75, 100};
}

void AppCross::load_phenotype()
{
    if ( m_par.pheno.empty() )
        return;

    std::cerr << "INFO: reading phenotype file...\n";

    int ret = read_phenotype(m_par.pheno, m_pt);

    if (ret != 0) {
        m_pt.phe.clear();
        m_pt.ind.clear();
    }

    std::cerr << "INFO: " << m_pt.ind.size() << " observations and " << m_pt.phe.size() << " phenotypes were observed\n";
}

void AppCross::load_genotype()
{
    if ( m_par.vcf.empty() && m_par.ped.empty() && m_par.hmp.empty() && m_par.geno.empty() )
        return;

    std::cerr << "INFO: reading genotype file...\n";

    int info = 0;

    if ( ! m_par.vcf.empty() )
        info = read_vcf(m_par.vcf, m_gt);
    else if ( ! m_par.ped.empty() ) {
        info = read_plink(m_par.ped, m_gt);
        auto m = m_gt.loc.size();
        for (size_t j = 0; j < m; ++j)
            m_gt.dist[j] *= 100.0;
    }
    else if ( ! m_par.hmp.empty() )
        info = read_hapmap(m_par.hmp, m_gt);
    else if ( ! m_par.geno.empty() )
        info = read_genotype(m_par.geno, m_gt);

    if (info != 0) {
        m_gt.loc.clear();
        m_gt.ind.clear();
        m_gt.dat.clear();
    }
    else
        update_chrlen(m_gt.chr, m_gt.dist);

    std::cerr << "INFO: " << m_gt.ind.size() << " individuals and " << m_gt.loc.size() << " loci were observed\n";
}

void AppCross::load_linakge_map()
{
    if ( m_par.map.empty() )
        return;

    std::cerr << "INFO: reading linkage map file...\n";

    vector<string> loc, chr;
    vector<double> pos;

    int info = read_map(m_par.map, loc, chr, pos);
    if (info != 0)
        loc.clear();

    std::cerr << "INFO: " << loc.size() << " loci were observed\n";

    if ( loc.empty() )
        return;

    update_chrlen(chr, pos);

    if (loc == m_gt.loc && chr == m_gt.chr) {
        m_gt.dist = pos;
        return;
    }

    auto m = m_gt.loc.size();

    for (size_t j = 0; j < m; ++j) {
        auto itr = std::find(loc.begin(), loc.end(), m_gt.loc[j]);
        if (itr == loc.end())
            continue;

        auto wh = std::distance(loc.begin(), itr);
        if (chr[wh] == m_gt.chr[j])
            m_gt.dist[j] = pos[wh];
    }
}

void AppCross::load_allele_effect()
{
    if ( m_par.effect.empty() )
        return;

    std::cerr << "INFO: reading allele effect file...\n";
    vector<string> ae_trait, ae_locus, ae_allele;
    vector<double> ae_effect;
    read_effect(m_par.effect, ae_trait, ae_locus, ae_allele, ae_effect);
    std::cerr << "INFO: " << ae_locus.size() << " records were observed\n";

    auto qtl = unique(ae_locus);
    m_trait = stable_unique(ae_trait);

    auto t = m_trait.size();
    auto m = m_gt.loc.size();

    vector<size_t> idx;
    for (size_t i = 0; i < m; ++i) {
        if ( std::binary_search(qtl.begin(), qtl.end(), m_gt.loc[i]) )
            idx.push_back(i);
    }

    subset(m_gt.loc,idx).swap(m_gt.loc);
    subset(m_gt.chr,idx).swap(m_gt.chr);
    subset(m_gt.pos,idx).swap(m_gt.pos);
    subset(m_gt.dist,idx).swap(m_gt.dist);
    subset(m_gt.dat,idx).swap(m_gt.dat);
    subset(m_gt.allele,idx).swap(m_gt.allele);

    m = m_gt.loc.size();

    m_effect.clear();
    for (size_t j = 0; j < m; ++j) {
        auto q = m_gt.allele[j].size();
        m_effect.emplace_back(q*t, 0.0);
    }

    auto n = ae_locus.size();
    for (size_t i = 0; i < n; ++i) {
        auto itr = std::find(m_gt.loc.begin(), m_gt.loc.end(), ae_locus[i]);
        if (itr == m_gt.loc.end()) {
            std::cerr << "WARNING: can't find genotype for locus: " << ae_locus[i] << "\n";
            continue;
        }

        auto j = std::distance(m_gt.loc.begin(), itr);
        const auto &allele = m_gt.allele[j];

        auto itr2 = std::find(allele.begin(), allele.end(), ae_allele[i]);
        if (itr2 == allele.end()) {
            std::cerr << "WARNING: can't find genotype for allele: " << ae_locus[i] << " -> " << ae_allele[i] << "\n";
            continue;
        }

        auto k = std::distance(allele.begin(), itr2);
        auto l = index(m_trait, ae_trait[i]);
        auto q = allele.size();

        m_effect[j][q*l+k] = ae_effect[i];
    }
}

void AppCross::save_allele_matrix() const
{
    std::cerr << "INFO: saving QTL-allele matrix...\n";

    std::ofstream ofs(m_par.out + ".qam.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".qam\n";
        return;
    }

    int p = m_gt.ploidy;
    bool diploid = p == 2;

    auto t = m_trait.size();
    auto n = m_gt.ind.size();
    auto m = m_gt.loc.size();

    ofs << "Locus\tChromosome\tPosition\tDistance\tAllele";
    for (size_t i = 0; i < t; ++i)
        ofs << "\t" << m_trait[i];
    for (size_t i = 0; i < n; ++i)
        ofs << "\t" << m_gt.ind[i];
    ofs << "\n";

    string line;
    for (size_t j = 0; j < m; ++j) {
        int q = m_gt.allele[j].size();
        for (int k = 0; k < q; ++k) {
            ofs << m_gt.loc[j] << "\t" << m_gt.chr[j] << "\t" << m_gt.pos[j] << "\t"
                << m_gt.dist[j] << "\t" << m_gt.allele[j][k];
            for (size_t i = 0; i < t; ++i)
                ofs << "\t" << m_effect[j][q*i+k];
            line.resize(0);
            for (size_t i = 0; i < n; ++i) {
                line.push_back('\t');
                auto a = m_gt.dat[j][i*p];
                auto b = diploid ? m_gt.dat[j][i*p+1] : a;
                line.push_back((a == k+1 || b == k+1) ? '1' : '0');
            }
            ofs << line << "\n";
        }
    }

    std::cerr << "INFO: QTL-allele matrix was successfully saved\n";
}

void AppCross::update_chrlen(const vector<string> &chr, const vector<double> &pos)
{
    auto n = chr.size();
    for (auto &e : unique(chr)) {
        double len = 0.0;
        for (size_t i = 0; i < n; ++i) {
            if (chr[i] == e && pos[i] > len)
                len = pos[i];
        }
        auto itr = m_chrlen.find(e);
        if (itr == m_chrlen.end())
            m_chrlen[e] = len;
        else if (len > itr->second)
            itr->second = len;
    }
}

bool AppCross::check_linkage_model() const
{
    if ( m_par.indep )
        return false;

    auto uchr = unique(m_gt.chr);

    for (auto e : m_chrlen) {
        if ( std::binary_search(uchr.begin(), uchr.end(), e.first) && e.second <= 0.0 ) {
            std::cerr << "WARNING: incorrect chromosome length: " << e.first << ": " << e.second << "\n";
            return false;
        }
    }

    for (auto e : m_gt.dist) {
        if ( ! std::isfinite(e) || e < 0.0 ) {
            std::cerr << "WARNING: invalid genetic distance: " << e << "\n";
            return false;
        }
    }

    for (auto &e : uchr) {
        vector<double> pos;
        auto m = m_gt.loc.size();
        for (size_t i = 0; i < m; ++i) {
            if (m_gt.chr[i] == e)
                pos.push_back(m_gt.dist[i]);
        }
        if ( ! std::is_sorted(pos.begin(), pos.end()) ) {
            std::cerr << "WARNING: genetic distance is not sorted\n";
            return false;
        }
        if (unique(pos).size() != pos.size()) {
            std::cerr << "WARNING: duplicate genetic distance detected\n";
            return false;
        }
    }

    return true;
}

void AppCross::make_genetic_map()
{
    auto m = m_gt.loc.size();

    m_gmap.chr.clear();
    m_gmap.len.clear();
    m_gmap.idx.clear();
    m_gmap.pos.clear();

    for (auto &chr : unique(m_gt.chr)) {
        vector<size_t> idx;
        vector<double> pos;
        for (size_t i = 0; i < m; ++i) {
            if (m_gt.chr[i] == chr) {
                idx.push_back(i);
                pos.push_back(m_gt.dist[i]);
            }
        }
        m_gmap.chr.push_back(chr);
        m_gmap.len.push_back(m_chrlen[chr]);
        m_gmap.idx.push_back(idx);
        m_gmap.pos.push_back(pos);
    }
}

void AppCross::predict_parent()
{
    m_pred.ind = m_gt.ind;
    m_pred.phe = m_trait;
    calc_genotypic_value(m_gt.dat, m_pred.dat);

    Phenotype pt = m_pred;

    int m = pt.phe.size();
    int n = pt.ind.size();

    for (int j = 0; j < m; ++j) {
        auto itr = std::find(m_pt.phe.begin(), m_pt.phe.end(), pt.phe[j]);
        if (itr == m_pt.phe.end()) {
            std::cerr << "WARNING: can't find phenotype data for trait: " << pt.phe[j] << "\n";
            continue;
        }
        auto k = std::distance(m_pt.phe.begin(), itr);
        for (int i = 0; i < n; ++i) {
            auto itr2 = std::find(m_pt.ind.begin(), m_pt.ind.end(), pt.ind[i]);
            if (itr2 != m_pt.ind.end()) {
                auto l = std::distance(m_pt.ind.begin(), itr2);
                pt.dat[j][i] = m_pt.dat[k][l];
            }
            else
                std::cerr << "WARNING: can't find phenotype data for individual: " << pt.ind[i] << "\n";
        }
    }

    m_pt.phe.swap(pt.phe);
    m_pt.ind.swap(pt.ind);
    m_pt.dat.swap(pt.dat);
}

void AppCross::cross_single_indep()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            simulate_indep(popsize, i, j, geno);

            calc_genotypic_value(geno, ys);

            adjust_prediction(i, j, ys);

            ofs << m_gt.ind[i] << "\t" << m_gt.ind[j];

            for (auto &v : ys) {
                ofs << "\t" << mean(v) << "\t" << sd(v);
                std::sort(v.begin(), v.end());
                for (auto e : m_perc)
                    ofs << "\t" << percentile(e,v);
            }

            ofs << "\n";
        }
    }
}

void AppCross::cross_single_linkage()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            simulate_linkage(popsize, i, j, geno);

            calc_genotypic_value(geno, ys);

            adjust_prediction(i, j, ys);

            ofs << m_gt.ind[i] << "\t" << m_gt.ind[j];

            for (auto &v : ys) {
                ofs << "\t" << mean(v) << "\t" << sd(v);
                std::sort(v.begin(), v.end());
                for (auto e : m_perc)
                    ofs << "\t" << percentile(e,v);
            }

            ofs << "\n";
        }
    }
}

void AppCross::cross_threeway_indep()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2\tP3";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                simulate_indep(popsize, i, j, k, geno);

                calc_genotypic_value(geno, ys);

                adjust_prediction(i, j, k, ys);

                ofs << m_gt.ind[i] << "\t" << m_gt.ind[j] << "\t" << m_gt.ind[k];

                for (auto &v : ys) {
                    ofs << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : m_perc)
                        ofs << "\t" << percentile(e,v);
                }

                ofs << "\n";
            }
        }
    }
}

void AppCross::cross_threeway_linkage()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2\tP3";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                simulate_linkage(popsize, i, j, k, geno);

                calc_genotypic_value(geno, ys);

                adjust_prediction(i, j, k, ys);

                ofs << m_gt.ind[i] << "\t" << m_gt.ind[j] << "\t" << m_gt.ind[k];

                for (auto &v : ys) {
                    ofs << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : m_perc)
                        ofs << "\t" << percentile(e,v);
                }

                ofs << "\n";
            }
        }
    }
}

void AppCross::cross_double_indep()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2\tP3\tP4";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                for (int l = k + 1; l < n; ++l) {
                    simulate_indep(popsize, i, j, k, l, geno);

                    calc_genotypic_value(geno, ys);

                    adjust_prediction(i, j, k, l, ys);

                    ofs << m_gt.ind[i] << "\t" << m_gt.ind[j] << "\t" << m_gt.ind[k] << "\t" << m_gt.ind[l];

                    for (auto &v : ys) {
                        ofs << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : m_perc)
                            ofs << "\t" << percentile(e,v);
                    }

                    ofs << "\n";
                }
            }
        }
    }
}

void AppCross::cross_double_linkage()
{
    int n = m_gt.ind.size();
    int t = m_trait.size();
    int popsize = m_par.size;

    vector< vector<allele_t> > geno;
    vector< vector<double> > ys;

    std::ofstream ofs(m_par.out + ".cross.txt");

    if ( ! ofs ) {
        std::cerr << "ERROR: can't open file: " << m_par.out << ".cross\n";
        return;
    }

    ofs << "P1\tP2\tP3\tP4";
    for (int i = 0; i < t; ++i) {
        auto s = m_trait[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : m_perc)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                for (int l = k + 1; l < n; ++l) {
                    simulate_linkage(popsize, i, j, k, l, geno);

                    calc_genotypic_value(geno, ys);

                    adjust_prediction(i, j, k, l, ys);

                    ofs << m_gt.ind[i] << "\t" << m_gt.ind[j] << "\t" << m_gt.ind[k] << "\t" << m_gt.ind[l];

                    for (auto &v : ys) {
                        ofs << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : m_perc)
                            ofs << "\t" << percentile(e,v);
                    }

                    ofs << "\n";
                }
            }
        }
    }
}

void AppCross::adjust_prediction(int p1, int p2, vector< vector<double> > &ys) const
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = m_pt.dat[j][p1];
        auto y2 = m_pt.dat[j][p2];
        auto g1 = m_pred.dat[j][p1];
        auto g2 = m_pred.dat[j][p2];
        auto a = (y1 - g1 + y2 - g2) / 2;
        for (auto &e : ys[j])
            e += a;
    }
}

void AppCross::adjust_prediction(int p1, int p2, int p3, vector< vector<double> > &ys) const
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = m_pt.dat[j][p1];
        auto y2 = m_pt.dat[j][p2];
        auto y3 = m_pt.dat[j][p3];
        auto g1 = m_pred.dat[j][p1];
        auto g2 = m_pred.dat[j][p2];
        auto g3 = m_pred.dat[j][p3];
        auto a = (y1 - g1 + y2 - g2 + y3 - g3) / 3;
        for (auto &e : ys[j])
            e += a;
    }
}

void AppCross::adjust_prediction(int p1, int p2, int p3, int p4, vector< vector<double> > &ys) const
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = m_pt.dat[j][p1];
        auto y2 = m_pt.dat[j][p2];
        auto y3 = m_pt.dat[j][p3];
        auto y4 = m_pt.dat[j][p4];
        auto g1 = m_pred.dat[j][p1];
        auto g2 = m_pred.dat[j][p2];
        auto g3 = m_pred.dat[j][p3];
        auto g4 = m_pred.dat[j][p4];
        auto a = (y1 - g1 + y2 - g2 + y3 - g3 + y4 - g4) / 4;
        for (auto &e : ys[j])
            e += a;
    }
}

void AppCross::simulate_indep(int n, int p1, int p2, vector< vector<allele_t> > &geno) const
{
    int p = m_gt.ploidy;
    auto m = m_gt.loc.size();

    int idx[4];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    if (p == 2) {
        idx[1] += 1;
        idx[3] += 1;
    }

    geno.resize(m);

    allele_t gam[4];
    std::uniform_int_distribution<> dis(0,3);

    for (size_t j = 0; j < m; ++j) {
        for (int i = 0; i < 4; ++i)
            gam[i] = m_gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homozygous = true;
        for (int i = 1; i < 4 && homozygous; ++i)
            homozygous = gam[i] == gam[0];

        if ( homozygous ) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[dis(RNG)];
        }
    }
}

void AppCross::simulate_linkage(int n, int p1, int p2, vector< vector<allele_t> > &geno) const
{
    int p = m_gt.ploidy;
    auto m = m_gt.loc.size();
    int r = m_gmap.chr.size();

    int idx[4];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    if (p == 2) {
        idx[1] += 1;
        idx[3] += 1;
    }

    geno.resize(m);
    for (auto &v : geno)
        v.resize(n);

    vector<allele_t> dad, mom, kid;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < r; ++j) {
            dad.clear();
            mom.clear();
            for (auto k : m_gmap.idx[j]) {
                dad.push_back(m_gt.dat[k][idx[0]]);
                dad.push_back(m_gt.dat[k][idx[1]]);
                mom.push_back(m_gt.dat[k][idx[2]]);
                mom.push_back(m_gt.dat[k][idx[3]]);
            }
            for (;;) {
                docross(m_gmap.len[j], m_gmap.pos[j], dad, mom, kid);
                if ( is_homozygous(kid) )
                    break;
                dad = mom = kid;
            }
            int t = 0;
            for (auto k : m_gmap.idx[j]) {
                geno[k][i] = kid[t*2];
                ++t;
            }
        }
    }
}

void AppCross::simulate_indep(int n, int p1, int p2, int p3, vector< vector<allele_t> > &geno) const
{
    int p = m_gt.ploidy;
    auto m = m_gt.loc.size();

    int idx[6];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    idx[4] = idx[5] = p*p3;
    if (p == 2) {
        ++idx[1];
        ++idx[3];
        ++idx[5];
    }

    geno.resize(m);

    allele_t gam[6];
    std::uniform_int_distribution<> dis(0,5);

    for (size_t j = 0; j < m; ++j) {
        for (int i = 0; i < 6; ++i)
            gam[i] = m_gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homozygous = true;
        for (int i = 1; i < 6 && homozygous; ++i)
            homozygous = gam[i] == gam[0];

        if ( homozygous ) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[dis(RNG)];
        }
    }
}

void AppCross::simulate_linkage(int n, int p1, int p2, int p3, vector< vector<allele_t> > &geno) const
{
    // TODO: implement
}

void AppCross::simulate_indep(int n, int p1, int p2, int p3, int p4, vector< vector<allele_t> > &geno) const
{
    int p = m_gt.ploidy;
    auto m = m_gt.loc.size();

    int idx[8];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    idx[4] = idx[5] = p*p3;
    idx[6] = idx[7] = p*p4;
    if (p == 2) {
        ++idx[1];
        ++idx[3];
        ++idx[5];
        ++idx[7];
    }

    geno.resize(m);

    allele_t gam[8];
    std::uniform_int_distribution<> dis(0,7);

    for (size_t j = 0; j < m; ++j) {
        for (int i = 0; i < 8; ++i)
            gam[i] = m_gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homozygous = true;
        for (int i = 1; i < 8 && homozygous; ++i)
            homozygous = gam[i] == gam[0];

        if ( homozygous ) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[dis(RNG)];
        }
    }
}

void AppCross::simulate_linkage(int n, int p1, int p2, int p3, int p4, vector< vector<allele_t> > &geno) const
{
    // TODO: implement
}

void AppCross::calc_genotypic_value(const vector< vector<allele_t> > &geno, vector< vector<double> > &pheno) const
{
    auto t = m_trait.size();
    auto m = m_gt.loc.size();
    auto n = geno.empty() ? 0 : geno[0].size();

    pheno.clear();
    for (size_t i = 0; i < t; ++i)
        pheno.emplace_back(n, 0.0);

    for (size_t j = 0; j < m; ++j) {
        auto q = m_gt.allele[j].size();
        for (size_t i = 0; i < n; ++i) {
            auto a = geno[j][i] - 1;
            for (size_t k = 0; k < t; ++k)
                pheno[k][i] += m_effect[j][k*q+a];
        }
    }
}

// Haldane, J.B.S. The combination of linkage values, and the calculation of
//   distance between the loci of linked factors. J Genet 8, 299-309 (1919).
void AppCross::crossover(double len, vector<double> &xo) const
{
    int n = std::poisson_distribution<>(len/100)(RNG);
    std::uniform_real_distribution<> dis(0.0, len);

    xo.resize(n);
    for (int i = 0; i < n; ++i)
        xo[i] = dis(RNG);

    std::sort(xo.begin(), xo.end());
}

void AppCross::meiosis(const vector<double> &xo, const vector<double> &pos, vector<bool> &gam) const
{
    auto n = pos.size();

    bool strand = std::uniform_real_distribution<>()(RNG) < 0.5;
    gam.assign(n, strand);

    if ( xo.empty() )
        return;

    size_t j = 0;

    for (auto e : xo) {
        for (size_t i = j; i < n; ++i) {
            if (pos[i] < e) {
                gam[i] = strand;
                ++j;
            }
            else
                break;
        }
        strand = !strand;
    }

    for (size_t i = j; i < n; ++i)
        gam[i] = strand;
}

void AppCross::docross(double len, const vector<double> &pos, const vector<allele_t> &dad,
                       const vector<allele_t> &mom, vector<allele_t> &kid) const
{
    auto n = pos.size();
    kid.resize(dad.size());

    vector<bool> gam;
    vector<double> xo;

    crossover(len, xo);
    meiosis(xo, pos, gam);

    for (size_t i = 0; i < n; ++i) {
        int k = gam[i] ? 1 : 0;
        kid[i*2] = dad[i*2+k];
    }

    crossover(len, xo);
    meiosis(xo, pos, gam);

    for (size_t i = 0; i < n; ++i) {
        int k = gam[i] ? 1 : 0;
        kid[i*2+1] = mom[i*2+k];
    }
}
