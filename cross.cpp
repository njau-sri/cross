#include <cmath>
#include <limits>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "vcfio.h"
#include "cmdline.h"
#include "util.h"
#include "split.h"

namespace {

static const double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Par
{
    std::string vcf;
    std::string map;
    std::string pheno;
    std::string effect;
    std::string wt;
    std::string out;
    std::string pct;
    int type = 2;
    int size = 2000;
};

struct Phenotype
{
    std::vector<std::string> ind;
    std::vector<std::string> phe;
    std::vector<std::string> env;
    std::vector<std::string> blk;
    std::vector< std::vector<double> > dat;
};

struct LinkageMap
{
    std::vector<std::string> chr;
    std::vector<double> len;
    std::vector< std::vector<int> > idx;
    std::vector< std::vector<double> > pos;
};

Par par;
std::mt19937 rng;

int parse_command_line(int argc, char *argv[])
{
    CmdLine cmd("cross [options]");

    cmd.add("--vcf", "VCF genotype file", "");
    cmd.add("--map", "linkage map (cM) file", "");
    cmd.add("--effect", "allele effect file", "");
    cmd.add("--pheno", "phenotype data file", "");
    cmd.add("--wt", "trait weights file", "");
    cmd.add("--out", "output file prefix", "cross.out");
    cmd.add("--type", "cross type 2/3/4", "2");
    cmd.add("--size", "sample size", "2000");
    cmd.add("--pct", "sample percentiles", "0,25,50,75,100");

    if (argc < 2) {
        cmd.help();
        return 1;
    }

    cmd.parse(argc, argv);

    par.vcf = cmd.get("--vcf");
    par.map = cmd.get("--map");
    par.effect = cmd.get("--effect");
    par.pheno = cmd.get("--pheno");
    par.out = cmd.get("--out");
    par.wt = cmd.get("--wt");
    par.pct = cmd.get("--pct");

    par.type = std::stoi(cmd.get("--type"));
    if (par.type < 2 || par.type > 4) {
        std::cerr << "WARNING: invalid cross type: " << cmd.get("--type") << "\n";
        par.type = 2;
    }

    par.size = std::stoi(cmd.get("--size"));
    if (par.size < 1) {
        std::cerr << "WARNING: invalid sample size: " << cmd.get("--size") << "\n";
        par.size = 2000;
    }

    if (par.type != 2 && !par.map.empty()) {
        std::cerr << "ERROR: linkage model is not available for cross type: " << par.type << "\n";
        return 1;
    }

    return 0;
}

void parse_percent(const std::string &pcts, std::vector<int> &pcti)
{
    std::vector<std::string> vs;
    split(pcts, ", \t", vs);

    std::vector<int> z;

    for (auto &e : vs) {
        int x = std::stoi(e);
        if (x < 0 || x > 100) {
            std::cerr << "WARNING: invalid percent values: " << pcts << "\n";
            return;
        }
        z.push_back(x);
    }

    pcti.swap(z);
}

int read_pheno(const std::string &filename, Phenotype &pt)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    int ln = 0;
    std::vector<std::string> colnames;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);
        split(line, " \t", colnames);
        if (!colnames.empty())
            break;
    }

    std::vector<int> jphe;
    int jenv = -1, jblk = -1;
    int ncols = colnames.size();

    for (int j = 1; j < ncols; ++j) {
        if (colnames[j] == "_ENV_")
            jenv = j;
        else if (colnames[j] == "_BLK_")
            jblk = j;
        else
            jphe.push_back(j);
    }

    for (auto j : jphe)
        pt.phe.push_back(colnames[j]);

    std::vector< std::vector<double> > dat;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if (vs.empty())
            continue;

        if (vs.size() != colnames.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << ": "
                      << vs.size() << "!=" << colnames.size() << "\n";
            return 1;
        }

        pt.ind.push_back(vs[0]);

        if (jenv > 0)
            pt.env.push_back(vs[jenv]);

        if (jblk > 0)
            pt.blk.push_back(vs[jblk]);

        std::vector<double> v;

        for (auto j : jphe) {
            if (vs[j] == "?" || vs[j] == "NA" || vs[j] == ".")
                v.push_back(kNaN);
            else
                v.push_back(std::stod(vs[j]));
        }

        dat.push_back(v);
    }

    int m = pt.phe.size();
    int n = pt.ind.size();

    for (int j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i)
            v[i] = dat[i][j];
        pt.dat.push_back(v);
    }

    return 0;
}

double calc_pheno_mean(const std::string &wh, const std::vector<std::string> &ind, const std::vector<double> &dat)
{
    int n = 0;
    double y = 0.0;

    int m = ind.size();
    for (int i = 0; i < m; ++i) {
        if (ind[i] == wh && std::isfinite(dat[i])) {
            ++n;
            y += dat[i];
        }
    }

    return n == 0 ? kNaN : y / n;
}

int write_pheno(const Phenotype &pt, const std::string &filename)
{
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "ERROR: can't open file for writing: " << filename << "\n";
        return 1;
    }

    int m = pt.phe.size();
    int n = pt.ind.size();

    ofs << "Indiv";
    for (int j = 0; j < m; ++j)
        ofs << "\t" << pt.phe[j];
    ofs << "\n";

    for (int i = 0; i < n; ++i) {
        ofs << pt.ind[i];
        for (int j = 0; j < m; ++j)
            ofs << "\t" << pt.dat[j][i];
        ofs << "\n";
    }

    return 0;
}

int read_effect(const std::string &filename, std::vector<std::string> &trait, std::vector<std::string> &locus,
                std::vector<std::string> &allele, std::vector<double> &effect)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    int ln = 0;
    std::string curr;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if (vs.empty())
            continue;

        if (vs.size() == 1 && vs[0].find('>') == 0) {
            curr = vs[0].substr(1);
            continue;
        }

        if (vs.size() != 3)
            continue;

        if (vs[0].find("_ENV_") != std::string::npos)
            continue;

        if (vs[0].find("_BLK_") != std::string::npos)
            continue;

        if (!curr.empty()) {
            std::vector<std::string> as;
            split(vs[1], ":/|", as);
            if (as.size() == 1 || (as.size() == 2 && as[0] == as[1])) {
                trait.push_back(curr);
                locus.push_back(vs[0]);
                allele.push_back(as[0]);
                effect.push_back(std::stod(vs[2]));
            }
        }
    }

    return 0;
}

int make_qtl_allele_matrix(Genotype &gt, std::vector<std::string> &phe, std::vector< std::vector<double> > &eff)
{
    std::vector<std::string> fphe, floc, fall;
    std::vector<double> feff;

    std::cerr << "INFO: reading allele effect file...\n";
    int info = read_effect(par.effect, fphe, floc, fall, feff);
    if (info != 0)
        return 1;
    std::cerr << "INFO: " << floc.size() << " allele records were observed\n";

    auto qtl = unique(floc);
    phe = stable_unique(fphe);

    int p = phe.size();
    int m = gt.loc.size();

    std::vector<int> idx;
    for (int i = 0; i < m; ++i) {
        if (std::binary_search(qtl.begin(), qtl.end(), gt.loc[i]))
            idx.push_back(i);
    }

    subset(gt.loc, idx).swap(gt.loc);
    subset(gt.chr, idx).swap(gt.chr);
    subset(gt.pos, idx).swap(gt.pos);
    subset(gt.dat, idx).swap(gt.dat);
    subset(gt.allele, idx).swap(gt.allele);

    m = gt.loc.size();
    for (int j = 0; j < m; ++j) {
        int q = gt.allele[j].size();
        eff.emplace_back(p*q, 0);
    }

    int n = floc.size();
    for (int i = 0; i < n; ++i) {
        auto j = index(gt.loc, floc[i]);
        if (j == gt.loc.size()) {
            std::cerr << "ERROR: can't find genotype for locus: " << floc[i] << "\n";
            return 1;
        }
        const auto &allele = gt.allele[j];

        auto k = index(allele, fall[i]);
        if (k == allele.size())
            continue;

        auto t = index(phe, fphe[i]);
        auto q = allele.size();
        eff[j][t*q + k] = feff[i];
    }

    // !!! NOTICE !!!
    // The specified genotype may be different from the genotype used for
    // effect estimation, therefore, allele may not have an effect.

    return 0;
}

int read_map(const std::string &filename, std::vector<std::string> &loc, std::vector<std::string> &chr, std::vector<double> &pos)
{
    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    int ln = 0;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if (vs.empty())
            continue;

        if (vs.size() < 3) {
            std::cerr << "ERROR: expected at least 3 columns at line " << ln << "\n";
            return 1;
        }

        loc.push_back(vs[0]);
        chr.push_back(vs[1]);
        pos.push_back(std::stod(vs[2]));
    }

    return 0;
}

int make_linkage_map(const Genotype &gt, LinkageMap &lm)
{
    std::vector<std::string> floc, fchr;
    std::vector<double> fpos;

    std::cerr << "INFO: reading linkage map file...\n";
    int info = read_map(par.map, floc, fchr, fpos);
    if (info != 0)
        return 1;
    std::cerr << "INFO: " << floc.size() << " map records were observed\n";

    std::map<std::string, double> chrlen;

    for (auto &e : unique(fchr)) {
        double len = 0.0;
        int n = fchr.size();
        for (int i = 0; i < n; ++i) {
            if (fchr[i] == e && fpos[i] > len)
                len = fpos[i];
        }
        auto itr = chrlen.find(e);
        if (itr == chrlen.end())
            chrlen[e] = len;
        else if (len > itr->second)
            itr->second = len;
    }

    for (auto &chr : unique(gt.chr)) {
        if (chrlen.find(chr) == chrlen.end()) {
            std::cerr << "ERROR: can't find linkage map for chromosome: " << chr << "\n";
            return 1;
        }

        if (chrlen[chr] <= 0.0) {
            std::cerr << "ERROR: incorrect chromosome length: " << chr << " " << chrlen[chr] << "\n";
            return 1;
        }

        std::vector<int> idx;
        std::vector<double> pos;

        int m = gt.loc.size();
        for (int i = 0; i < m; ++i) {
            if (gt.chr[i] != chr)
                continue;

            auto j = index(floc, gt.loc[i]);

            if (j == floc.size()) {
                std::cerr << "ERROR: can't find linkage map position for locus: " << gt.loc[i] << "\n";
                return 1;
            }

            if (fchr[j] != chr) {
                std::cerr << "ERROR: chromosome doesn't match for locus: "
                          << gt.loc[i] << " " << chr << " " << fchr[j] << "\n";
                return 1;
            }

            if (!std::isfinite(fpos[j]) || fpos[j] < 0.0) {
                std::cerr << "ERROR: invalid linkage map position: " << fpos[j] << "\n";
                return 1;
            }

            idx.push_back(i);
            pos.push_back(fpos[j]);
        }

        if (!std::is_sorted(pos.begin(), pos.end())) {
            auto ord = order(pos);
            subset(idx, ord).swap(idx);
            subset(pos, ord).swap(pos);
        }

        if (unique(pos).size() != pos.size()) {
            std::cerr << "ERROR: duplicate positions in linkage map are not allowed\n";
            return 1;
        }

        lm.chr.push_back(chr);
        lm.len.push_back(chrlen[chr]);
        lm.idx.push_back(idx);
        lm.pos.push_back(pos);
    }

    return 0;
}

int read_weights(const std::string &filename, Phenotype &wt)
{
    // phe -> composite trait names
    // ind -> component trait names
    // dat -> weights

    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
        return 1;
    }

    int ln = 0;
    std::vector<std::string> colnames;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);
        split(line, " \t", colnames);
        if (!colnames.empty())
            break;
    }

    wt.phe.insert(wt.phe.end(), colnames.begin() + 1, colnames.end());

    std::vector< std::vector<double> > dat;

    for (std::string line; std::getline(ifs, line); ) {
        ++ln;
        line.erase(line.find_last_not_of("\r\n") + 1);

        std::vector<std::string> vs;
        split(line, " \t", vs);
        if (vs.empty())
            continue;

        if (vs.size() != colnames.size()) {
            std::cerr << "ERROR: column count doesn't match at line " << ln << ": "
                      << vs.size() << "!=" << colnames.size() << "\n";
            return 1;
        }

        wt.ind.push_back(vs[0]);

        std::vector<double> v;
        for (auto itr = vs.begin() + 1; itr != vs.end(); ++itr)
            v.push_back(std::stod(*itr));

        dat.push_back(v);
    }

    int m = wt.phe.size();
    int n = wt.ind.size();

    for (int j = 0; j < m; ++j) {
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i)
            v[i] = dat[i][j];
        wt.dat.push_back(v);
    }

    return 0;
}

int format_weights(const std::vector<std::string> &phe, Phenotype &wt)
{
    int t = 0;

    for (auto e : wt.ind) {
        if (std::find(phe.begin(), phe.end(), e) == phe.end())
            std::cerr << "WARNING: component trait is not included in prediction: " << e << "\n";
        else
            ++t;
    }

    if (t == 0) {
        std::cerr << "WARNING: no valid component trait found\n";
        wt.phe.clear();
    }

    int n = phe.size();
    int m = wt.phe.size();

    std::vector< std::vector<double> > dat(m, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        auto k = index(wt.ind, phe[i]);
        if (k != wt.ind.size()) {
            for (int j = 0; j < m; ++j)
                dat[j][i] = wt.dat[j][k];
        }
    }

    wt.ind = phe;
    wt.dat.swap(dat);

    return 0;
}

void calc_add(const Genotype &gt, const std::vector< std::vector<double> > &eff,
              const std::vector< std::vector<allele_t> > &geno, std::vector< std::vector<double> > &pheno)
{
    int m = gt.loc.size();
    int n = geno.empty() ? 0 : geno[0].size();

    int p = eff.empty() || gt.allele.empty() ? 0 : eff[0].size() / gt.allele[0].size();

    pheno.clear();
    for (int i = 0; i < p; ++i)
        pheno.emplace_back(n, kNaN);

    for (int j = 0; j < m; ++j) {
        int q = gt.allele[j].size();
        for (int i = 0; i < n; ++i) {
            int a = geno[j][i] - 1;
            for (int t = 0; t < p; ++t) {
                auto& y = pheno[t][i];
                auto x = eff[j][t*q+a];
                if ( std::isnan(y) )
                    y = x;
                else
                    y += x;
            }
        }
    }
}

// t = [0,100], 0 = min, 100 = max
double quantile(int t, const std::vector<double> &v)
{
    if (t <= 0)
        return v.front();

    if (t >= 100)
        return v.back();

    double np = v.size() * t / 100.0;
    double j = 0;
    double g = std::modf(np, &j);

    if (g == 0.0)
        return (v[j-1] + v[j]) / 2;

    return v[j];
}

double mean(const std::vector<double> &x)
{
    return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

double sd(const std::vector<double> &x)
{
    auto m = mean(x);
    double ss = 0;
    for (auto e : x)
        ss += (e - m)*(e - m);
    return std::sqrt(ss / (x.size() - 1));
}

bool is_pure_line(const std::vector<allele_t> &v)
{
    int n = v.size() / 2;

    for (int i = 0; i < n; ++i) {
        if (v[i * 2] != v[i * 2 + 1])
            return false;
    }

    return true;
}

// Haldane, J.B.S. The combination of linkage values, and the calculation of
//   distance between the loci of linked factors. J Genet 8, 299-309 (1919).
void crossover(double len, std::vector<double> &xo)
{
    int n = std::poisson_distribution<>(len / 100)(rng);
    std::uniform_real_distribution<> unif(0.0, len);

    xo.resize(n);
    for (int i = 0; i < n; ++i)
        xo[i] = unif(rng);

    std::sort(xo.begin(), xo.end());
}

void meiosis(const std::vector<double> &xo, const std::vector<double> &pos, std::vector<bool> &gam)
{
    int n = pos.size();

    bool strand = std::uniform_real_distribution<>()(rng) < 0.5;
    gam.assign(n, strand);

    if (xo.empty())
        return;

    int j = 0;

    for (auto e : xo) {
        for (int i = j; i < n; ++i) {
            if (pos[i] < e) {
                gam[i] = strand;
                ++j;
            }
            else
                break;
        }
        strand = !strand;
    }

    for (int i = j; i < n; ++i)
        gam[i] = strand;
}

struct Recombinant
{
    std::vector<bool> gam;
    std::vector<double> xo;

    void docross(double len, const std::vector<double> &pos,
                 const std::vector<allele_t> &dad, const std::vector<allele_t> &mom,
                 std::vector<allele_t> &kid)
    {
        int n = pos.size();
        kid.resize(dad.size());

        crossover(len, xo);
        meiosis(xo, pos, gam);

        for (int i = 0; i < n; ++i) {
            int k = gam[i] ? 1 : 0;
            kid[i * 2] = dad[i * 2 + k];
        }

        crossover(len, xo);
        meiosis(xo, pos, gam);

        for (int i = 0; i < n; ++i) {
            int k = gam[i] ? 1 : 0;
            kid[i * 2 + 1] = mom[i * 2 + k];
        }
    }
};

void simulate_indep(int n, int p1, int p2, const Genotype &gt, std::vector< std::vector<allele_t> > &geno)
{
    int p = gt.ploidy;
    int m = gt.loc.size();

    int idx[4];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    if (p == 2) {
        ++idx[1];
        ++idx[3];
    }

    geno.resize(m);

    allele_t gam[4];
    std::uniform_int_distribution<> unif(0, 3);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < 4; ++i)
            gam[i] = gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homo = true;
        for (int i = 1; i < 4 && homo; ++i)
            homo = gam[i] == gam[0];

        if (homo) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[unif(rng)];
        }
    }
}

void simulate_linkage(int n, int p1, int p2, const Genotype &gt, const LinkageMap &lm, std::vector< std::vector<allele_t> > &geno)
{
    int p = gt.ploidy;
    auto m = gt.loc.size();
    int nchr = lm.chr.size();

    int idx[4];
    idx[0] = idx[1] = p*p1;
    idx[2] = idx[3] = p*p2;
    if (p == 2) {
        ++idx[1];
        ++idx[3];
    }

    geno.resize(m);
    for (auto &v : geno)
        v.resize(n);

    Recombinant rec;
    std::vector<allele_t> dad, mom, kid;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nchr; ++j) {
            dad.clear();
            mom.clear();
            for (auto k : lm.idx[j]) {
                dad.push_back(gt.dat[k][idx[0]]);
                dad.push_back(gt.dat[k][idx[1]]);
                mom.push_back(gt.dat[k][idx[2]]);
                mom.push_back(gt.dat[k][idx[3]]);
            }
            for (;;) {
                rec.docross(lm.len[j], lm.pos[j], dad, mom, kid);
                if (is_pure_line(kid))
                    break;
                dad = mom = kid;
            }
            int t = 0;
            for (auto k : lm.idx[j]) {
                geno[k][i] = kid[t * 2];
                ++t;
            }
        }
    }
}

void simulate_indep(int n, int p1, int p2, int p3, const Genotype &gt, std::vector< std::vector<allele_t> > &geno)
{
    int p = gt.ploidy;
    int m = gt.loc.size();

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
    std::uniform_int_distribution<> unif(0, 5);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < 6; ++i)
            gam[i] = gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homo = true;
        for (int i = 1; i < 6 && homo; ++i)
            homo = gam[i] == gam[0];

        if (homo) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[unif(rng)];
        }
    }
}

void simulate_indep(int n, int p1, int p2, int p3, int p4, const Genotype &gt, std::vector< std::vector<allele_t> > &geno)
{
    int p = gt.ploidy;
    int m = gt.loc.size();

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
    std::uniform_int_distribution<> unif(0, 7);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < 8; ++i)
            gam[i] = gt.dat[j][idx[i]];

        auto &v = geno[j];
        v.resize(n);

        bool homo = true;
        for (int i = 1; i < 8 && homo; ++i)
            homo = gam[i] == gam[0];

        if (homo) {
            std::fill(v.begin(), v.end(), gam[0]);
        }
        else {
            for (int i = 0; i < n; ++i)
                v[i] = gam[unif(rng)];
        }
    }
}

void adjust_pred(int p1, int p2, const Phenotype &pv, const Phenotype &gv, std::vector< std::vector<double> > &ys)
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = pv.dat[j][p1];
        auto y2 = pv.dat[j][p2];
        auto g1 = gv.dat[j][p1];
        auto g2 = gv.dat[j][p2];
        auto a = (y1 - g1 + y2 - g2) / 2;
        std::for_each(ys[j].begin(), ys[j].end(), [a](double &x) { x += a; });
    }
}

void adjust_pred(int p1, int p2, int p3, const Phenotype &pv, const Phenotype &gv, std::vector< std::vector<double> > &ys)
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = pv.dat[j][p1];
        auto y2 = pv.dat[j][p2];
        auto y3 = pv.dat[j][p3];
        auto g1 = gv.dat[j][p1];
        auto g2 = gv.dat[j][p2];
        auto g3 = gv.dat[j][p3];
        auto a = (y1 - g1 + y2 - g2 + y3 - g3) / 3;
        std::for_each(ys[j].begin(), ys[j].end(), [a](double &x) { x += a; });
    }
}

void adjust_pred(int p1, int p2, int p3, int p4, const Phenotype &pv, const Phenotype &gv, std::vector< std::vector<double> > &ys)
{
    int m = ys.size();
    for (int j = 0; j < m; ++j) {
        auto y1 = pv.dat[j][p1];
        auto y2 = pv.dat[j][p2];
        auto y3 = pv.dat[j][p3];
        auto y4 = pv.dat[j][p4];
        auto g1 = gv.dat[j][p1];
        auto g2 = gv.dat[j][p2];
        auto g3 = gv.dat[j][p3];
        auto g4 = gv.dat[j][p4];
        auto a = (y1 - g1 + y2 - g2 + y3 - g3 + y4 - g4) / 4;
        std::for_each(ys[j].begin(), ys[j].end(), [a](double &x) { x += a; });
    }
}

void calc_composite(const std::vector< std::vector<double> > &ys, const std::vector<double> &w, std::vector<double> &z)
{
    int m = ys.size();
    int n = ys.empty() ? 0 : ys[0].size();
    z.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j)
            z[i] += w[j] * ys[j][i];
    }
}

void cross_pred_2(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm,
                  const Phenotype &pv, const Phenotype &gv, const Phenotype &wt,
                  const std::vector< std::vector<double> > &eff)
{
    int n = gt.ind.size();
    int t = gv.phe.size();

    std::vector< std::vector<allele_t> > geno;
    std::vector< std::vector<double> > ys, ys2, ys3;

    std::ofstream ofs(par.out + ".pred");
    if (!ofs) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
        return;
    }

    ofs << "P1\tP2";
    for (int i = 0; i < t; ++i) {
        auto s = gv.phe[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : pct)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    std::ofstream ofs2;
    if (!pv.ind.empty()) {
        ofs2.open(par.out + ".pred.adj");
        if (!ofs2) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
            return;
        }
        ofs2 << "P1\tP2";
        for (int i = 0; i < t; ++i) {
            auto s = gv.phe[i];
            ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs2 << "\t" << s << ".P" << e;
        }
        ofs2 << "\n";
    }

    std::ofstream ofs3;
    if (!wt.phe.empty()) {
        ofs3.open(par.out + ".compo");
        if (!ofs3) {
            std::cerr << "ERROR: can't open file: " << par.out << ".compo\n";
            return;
        }
        ofs3 << "P1\tP2";
        for (int i = 0; i < t; ++i) {
            auto s = wt.phe[i];
            ofs3 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs3 << "\t" << s << ".P" << e;
        }
        ofs3 << "\n";
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (lm.chr.empty())
                simulate_indep(par.size, i, j, gt, geno);
            else
                simulate_linkage(par.size, i, j, gt, lm, geno);

            calc_add(gt, eff, geno, ys);

            if (!pv.ind.empty()) {
                ys2 = ys;
                adjust_pred(i, j, pv, gv, ys2);
            }

            if (!wt.phe.empty()) {
                ys3.clear();
                std::vector<double> z;
                for (auto &w : wt.dat) {
                    if (ys2.empty())
                        calc_composite(ys, w, z);
                    else
                        calc_composite(ys2, w, z);
                    ys3.push_back(z);
                }
            }

            ofs << gt.ind[i] << "\t" << gt.ind[j];
            for (auto &v : ys) {
                ofs << "\t" << mean(v) << "\t" << sd(v);
                std::sort(v.begin(), v.end());
                for (auto e : pct)
                    ofs << "\t" << quantile(e, v);
            }
            ofs << "\n";

            if (ofs2.is_open()) {
                ofs2 << gt.ind[i] << "\t" << gt.ind[j];
                for (auto &v : ys2) {
                    ofs2 << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : pct)
                        ofs2 << "\t" << quantile(e, v);
                }
                ofs2 << "\n";
            }

            if (ofs3.is_open()) {
                ofs3 << gt.ind[i] << "\t" << gt.ind[j];
                for (auto &v : ys3) {
                    ofs3 << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : pct)
                        ofs3 << "\t" << quantile(e, v);
                }
                ofs3 << "\n";
            }
        }
    }
}

void cross_pred_3(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm,
                  const Phenotype &pv, const Phenotype &gv, const Phenotype &wt,
                  const std::vector< std::vector<double> > &eff)
{
    (void) lm;

    int n = gt.ind.size();
    int t = gv.phe.size();

    std::vector< std::vector<allele_t> > geno;
    std::vector< std::vector<double> > ys, ys2, ys3;

    std::ofstream ofs(par.out + ".pred");
    if (!ofs) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
        return;
    }

    ofs << "P1\tP2\tP3";
    for (int i = 0; i < t; ++i) {
        auto s = gv.phe[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : pct)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    std::ofstream ofs2;
    if (!pv.ind.empty()) {
        ofs2.open(par.out + ".pred.adj");
        if (!ofs2) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
            return;
        }
        ofs2 << "P1\tP2\tP3";
        for (int i = 0; i < t; ++i) {
            auto s = gv.phe[i];
            ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs2 << "\t" << s << ".P" << e;
        }
        ofs2 << "\n";
    }

    std::ofstream ofs3;
    if (!wt.phe.empty()) {
        ofs3.open(par.out + ".compo");
        if (!ofs3) {
            std::cerr << "ERROR: can't open file: " << par.out << ".compo\n";
            return;
        }
        ofs3 << "P1\tP2\tP3";
        for (int i = 0; i < t; ++i) {
            auto s = wt.phe[i];
            ofs3 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs3 << "\t" << s << ".P" << e;
        }
        ofs3 << "\n";
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                simulate_indep(par.size, i, j, k, gt, geno);

                calc_add(gt, eff, geno, ys);

                if (!pv.ind.empty()) {
                    ys2 = ys;
                    adjust_pred(i, j, k, pv, gv, ys2);
                }

                if (!wt.phe.empty()) {
                    ys3.clear();
                    std::vector<double> z;
                    for (auto &w : wt.dat) {
                        if (ys2.empty())
                            calc_composite(ys, w, z);
                        else
                            calc_composite(ys2, w, z);
                        ys3.push_back(z);
                    }
                }

                ofs << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k];
                for (auto &v : ys) {
                    ofs << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : pct)
                        ofs << "\t" << quantile(e, v);
                }
                ofs << "\n";

                if (ofs2.is_open()) {
                    ofs2 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k];
                    for (auto &v : ys2) {
                        ofs2 << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : pct)
                            ofs2 << "\t" << quantile(e, v);
                    }
                    ofs2 << "\n";
                }

                if (ofs3.is_open()) {
                    ofs3 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k];
                    for (auto &v : ys3) {
                        ofs3 << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : pct)
                            ofs3 << "\t" << quantile(e, v);
                    }
                    ofs3 << "\n";
                }
            }
        }
    }
}

void cross_pred_4(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm,
                  const Phenotype &pv, const Phenotype &gv, const Phenotype &wt,
                  const std::vector< std::vector<double> > &eff)
{
    (void) lm;

    int n = gt.ind.size();
    int t = gv.phe.size();

    std::vector< std::vector<allele_t> > geno;
    std::vector< std::vector<double> > ys, ys2, ys3;

    std::ofstream ofs(par.out + ".pred");
    if (!ofs) {
        std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
        return;
    }

    ofs << "P1\tP2\tP3\tP4";
    for (int i = 0; i < t; ++i) {
        auto s = gv.phe[i];
        ofs << "\t" << s << ".MEAN\t" << s << ".SD";
        for (auto e : pct)
            ofs << "\t" << s << ".P" << e;
    }
    ofs << "\n";

    std::ofstream ofs2;
    if (!pv.ind.empty()) {
        ofs2.open(par.out + ".pred.adj");
        if (!ofs2) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
            return;
        }
        ofs2 << "P1\tP2\tP3\tP4";
        for (int i = 0; i < t; ++i) {
            auto s = gv.phe[i];
            ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs2 << "\t" << s << ".P" << e;
        }
        ofs2 << "\n";
    }

    std::ofstream ofs3;
    if (!wt.phe.empty()) {
        ofs3.open(par.out + ".compo");
        if (!ofs3) {
            std::cerr << "ERROR: can't open file: " << par.out << ".compo\n";
            return;
        }
        ofs3 << "P1\tP2\tP3\tP4";
        for (int i = 0; i < t; ++i) {
            auto s = wt.phe[i];
            ofs3 << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs3 << "\t" << s << ".P" << e;
        }
        ofs3 << "\n";
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                for (int l = k + 1; l < n; ++l) {
                    simulate_indep(par.size, i, j, k, l, gt, geno);

                    calc_add(gt, eff, geno, ys);

                    if (!pv.ind.empty()) {
                        ys2 = ys;
                        adjust_pred(i, j, k, l, pv, gv, ys2);
                    }

                    if (ofs3.is_open()) {
                        ys3.clear();
                        std::vector<double> z;
                        for (auto &w : wt.dat) {
                            if (ys2.empty())
                                calc_composite(ys, w, z);
                            else
                                calc_composite(ys2, w, z);
                            ys3.push_back(z);
                        }
                    }

                    ofs << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k] << "\t" << gt.ind[l];
                    for (auto &v : ys) {
                        ofs << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : pct)
                            ofs << "\t" << quantile(e, v);
                    }
                    ofs << "\n";

                    if (ofs2.is_open()) {
                        ofs2 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k] << "\t" << gt.ind[l];
                        for (auto &v : ys2) {
                            ofs2 << "\t" << mean(v) << "\t" << sd(v);
                            std::sort(v.begin(), v.end());
                            for (auto e : pct)
                                ofs2 << "\t" << quantile(e, v);
                        }
                        ofs2 << "\n";
                    }

                    if (ofs3.is_open()) {
                        ofs3 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k] << "\t" << gt.ind[l];
                        for (auto &v : ys3) {
                            ofs3 << "\t" << mean(v) << "\t" << sd(v);
                            std::sort(v.begin(), v.end());
                            for (auto e : pct)
                                ofs3 << "\t" << quantile(e, v);
                        }
                        ofs3 << "\n";
                    }
                }
            }
        }
    }
}

void write_qtl_allele_matrix(const Genotype &gt, const std::vector<std::string> &traits, const std::vector< std::vector<double> > &effects)
{
    std::cerr << "INFO: saving QTL-allele matrix...\n";

    std::ofstream ofs(par.out + ".qam");
    if (!ofs) {
        std::cerr << "ERROR: can't open file for writing: " << par.out << ".qam\n";
        return;
    }

    int p = gt.ploidy;
    bool diploid = p == 2;

    int t = traits.size();
    int n = gt.ind.size();
    int m = gt.loc.size();

    ofs << "Locus\tChromosome\tPosition\tAllele";
    for (int i = 0; i < t; ++i)
        ofs << "\t" << traits[i];
    for (int i = 0; i < n; ++i)
        ofs << "\t" << gt.ind[i];
    ofs << "\n";

    std::string line;
    for (int j = 0; j < m; ++j) {
        int q = gt.allele[j].size();
        for (int k = 0; k < q; ++k) {
            ofs << gt.loc[j] << "\t" << gt.chr[j] << "\t" << gt.pos[j] << "\t" << gt.allele[j][k];
            for (int i = 0; i < t; ++i)
                ofs << "\t" << effects[j][q*i + k];
            line.resize(0);
            for (int i = 0; i < n; ++i) {
                line.push_back('\t');
                auto a = gt.dat[j][i*p];
                auto b = diploid ? gt.dat[j][i*p + 1] : a;
                line.push_back((a == k + 1 || b == k + 1) ? '1' : '0');
            }
            ofs << line << "\n";
        }
    }

    std::cerr << "INFO: QTL-allele matrix was successfully saved\n";
}

} // namespace

int cross(int argc, char *argv[])
{
    std::cerr << "CROSS (Built on " __DATE__ " " __TIME__ ")\n";

    int info = parse_command_line(argc, argv);
    if (info != 0)
        return 1;

    std::vector<int> pct = {0, 25, 50, 75, 100};
    parse_percent(par.pct, pct);

    Genotype gt;
    std::cerr << "INFO: reading genotype file...\n";
    info = read_vcf(par.vcf, gt);
    if (info != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::vector<std::string> phe;
    std::vector< std::vector<double> > eff;
    info = make_qtl_allele_matrix(gt, phe, eff);
    if (info != 0)
        return 1;

    Phenotype pt;
    if (!par.pheno.empty()) {
        std::cerr << "INFO: reading phenotype file...\n";
        info = read_pheno(par.pheno, pt);
        if (info != 0)
            return 1;
        std::cerr << "INFO: " << pt.ind.size() << " individuals, " << pt.phe.size() << " traits\n";
    }

    if (!pt.ind.empty()) {
        std::vector< std::vector<double> > dat;
        for (auto &curr_phe : phe) {
            auto j = index(pt.phe, curr_phe);
            if (j == pt.phe.size()) {
                std::cerr << "WARNING: can't find phenotype data for trait: " << curr_phe << "\n";
                dat.push_back(std::vector<double>(gt.ind.size(), kNaN));
                continue;
            }

            std::vector<double> v;
            for (auto &curr_ind : gt.ind) {
                double y = calc_pheno_mean(curr_ind, pt.ind, pt.dat[j]);
                if (std::isnan(y))
                    std::cerr << "WARNING: can't find phenotype data for individual: " << curr_ind << "\n";
                v.push_back(y);
            }

            dat.push_back(v);
        }

        pt.phe = phe;
        pt.ind = gt.ind;
        pt.dat = dat;
        pt.env.clear();
        pt.blk.clear();
    }

    Phenotype wt;
    if (!par.wt.empty()) {
        std::cerr << "INFO: reading trait weights file...\n";
        info = read_weights(par.wt, wt);
        if (info != 0)
            return 1;
        std::cerr << "INFO: " << wt.ind.size() << " traits, " << wt.phe.size() << " composite traits\n";
        format_weights(phe, wt);
    }

    LinkageMap lm;
    if (!par.map.empty()) {
        info = make_linkage_map(gt, lm);
        if (info != 0)
            return 1;
    }

    Phenotype gv;
    gv.phe = phe;
    gv.ind = gt.ind;
    if (gt.ploidy != 1) {
        auto dat = gt.dat;
        int n = gt.ind.size();
        for (auto &v : dat) {
            for (int i = 0; i < n; ++i)
                v[i] = v[i * gt.ploidy];
            v.resize(n);
        }
        calc_add(gt, eff, dat, gv.dat);
    }
    else
        calc_add(gt, eff, gt.dat, gv.dat);

    if (par.type == 3)
        cross_pred_3(pct, gt, lm, pt, gv, wt, eff);
    else if (par.type == 4)
        cross_pred_4(pct, gt, lm, pt, gv, wt, eff);
    else
        cross_pred_2(pct, gt, lm, pt, gv, wt, eff);

    write_qtl_allele_matrix(gt, phe, eff);

    write_pheno(gv, par.out + ".parent");

    std::cerr << "INFO: CROSS completed successfully\n";

    return 0;
}
