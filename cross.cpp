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

using std::size_t;

/// !!! TODO: parallel programming !!!

namespace {

    static const double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Par
    {
        std::string vcf;
        std::string map;
        std::string pheno;
        std::string effect;
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

    void parse_percentile(const std::string &pcts, std::vector<int> &pcti)
    {
        std::vector<std::string> vs;
        split(pcts, ", \t", vs);

        std::vector<int> z;

        for (auto &e : vs) {
            int x = std::stoi(e);
            if (x <= 0 || x > 100) {
                std::cerr << "WARNING: invalid percentile values: " << pcts << "\n";
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

        size_t ln = 0;
        std::vector<std::string> colnames;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);
            split(line, " \t", colnames);
            if (!colnames.empty())
                break;
        }

        std::vector<size_t> jphe;
        size_t jenv = 0, jblk = 0;

        for (size_t j = 1; j < colnames.size(); ++j) {
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
                std::cerr << "ERROR: column count doesn't match at line " << ln << " (" << vs.size() << "!=" << colnames.size() << "\n";
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

    int read_effect(const std::string &filename, std::vector<std::string> &trait, std::vector<std::string> &locus, std::vector<std::string> &allele, std::vector<double> &effect)
    {
        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "ERROR: can't open file for reading: " << filename << "\n";
            return 1;
        }

        size_t ln = 0;
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

            if (vs.size() < 3) {
                std::cerr << "ERROR: expected at least 3 columns at line: " << ln << "\n";
                return 1;
            }

            if (!curr.empty()) {
                std::vector<std::string> as;
                split(vs[1], ":/|", as);
                // !!! additive effects only !!!
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

    int make_qtl_allele_matrix(Genotype &gt, std::vector<std::string> &traits, std::vector< std::vector<double> > &effects)
    {
        std::cerr << "INFO: reading allele effect file...\n";

        std::vector<std::string> ae_trait, ae_locus, ae_allele;
        std::vector<double> ae_effect;

        if (read_effect(par.effect, ae_trait, ae_locus, ae_allele, ae_effect) != 0)
            return 1;

        std::cerr << "INFO: " << ae_locus.size() << " records were observed\n";

        auto qtl = unique(ae_locus);
        traits = stable_unique(ae_trait);

        int p = traits.size();
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
            effects.emplace_back(p*q, kNaN);
        }

        int n = ae_locus.size();
        for (int i = 0; i < n; ++i) {
            auto j = index(gt.loc, ae_locus[i]);
            if (j == gt.loc.size()) {
                std::cerr << "ERROR: can't find genotype for locus: " << ae_locus[i] << "\n";
                return 1;
            }
            const auto &allele = gt.allele[j];

            auto k = index(allele, ae_allele[i]);
            if (k == allele.size())
                continue;

            auto t = index(traits, ae_trait[i]);
            auto q = allele.size();
            effects[j][t*q + k] = ae_effect[i];
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

        size_t ln = 0;

        for (std::string line; std::getline(ifs, line); ) {
            ++ln;
            line.erase(line.find_last_not_of("\r\n") + 1);

            std::vector<std::string> vs;
            split(line, " \t", vs);
            if (vs.empty())
                continue;

            if (vs.size() < 3) {
                std::cerr << "ERROR: expected at least 3 columns at line: " << ln << "\n";
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
        std::cerr << "INFO: reading linkage map file...\n";

        std::vector<std::string> mloc, mchr;
        std::vector<double> mpos;

        if (read_map(par.map, mloc, mchr, mpos) != 0)
            return 1;

        std::cerr << "INFO: " << mloc.size() << " loci\n";

        std::map<std::string, double> chrlen;

        for (auto &e : unique(mchr)) {
            double len = 0.0;
            int n = mchr.size();
            for (int i = 0; i < n; ++i) {
                if (mchr[i] == e && mpos[i] > len)
                    len = mpos[i];
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

                auto j = index(mloc, gt.loc[i]);
                if (j == mloc.size() || mchr[j] != chr) {
                    std::cerr << "ERROR: can't find linkage map position for locus: " << gt.loc[i] << "\n";
                    return 1;
                }

                if (!std::isfinite(mpos[j]) || mpos[j] < 0.0) {
                    std::cerr << "ERROR: invalid linkage map position: " << mpos[j] << "\n";
                    return 1;
                }

                idx.push_back(i);
                pos.push_back(mpos[j]);
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

    void calc_add(const Genotype &gt, const std::vector< std::vector<double> > &effects, const std::vector< std::vector<allele_t> > &geno, std::vector< std::vector<double> > &pheno)
    {
        int m = gt.loc.size();
        int n = geno.empty() ? 0 : geno[0].size();

        int p = effects.empty() || gt.allele.empty() ? 0 : effects[0].size() / gt.allele[0].size();

        pheno.clear();
        for (int i = 0; i < p; ++i)
            pheno.emplace_back(n, kNaN);

        for (int j = 0; j < m; ++j) {
            int q = gt.allele[j].size();
            for (int i = 0; i < n; ++i) {
                int a = geno[j][i] - 1;
                for (int t = 0; t < p; ++t) {
                    auto& y = pheno[t][i];
                    auto x = effects[j][t*q+a];
                    if (std::isfinite(x)) {
                        if (std::isfinite(y))
                            y += x;
                        else
                            y = x;
                    }
                }
            }
        }
    }

    double quantile(int p, const std::vector<double> &v)
    {
        if (p <= 1)
            return v.front();

        if (p >= 100)
            return v.back();

        double x = v.size() * p / 100.0;
        auto j = static_cast<size_t>(x);

        if (x - j == 0.0)
            return (v[j - 1] + v[j]) / 2;

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

        void docross(double len, const std::vector<double> &pos, const std::vector<allele_t> &dad, const std::vector<allele_t> &mom, std::vector<allele_t> &kid)
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
        auto m = gt.loc.size();

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

        for (size_t j = 0; j < m; ++j) {
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
        auto m = gt.loc.size();

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

        for (size_t j = 0; j < m; ++j) {
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
        auto m = gt.loc.size();

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

        for (size_t j = 0; j < m; ++j) {
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

    void adjust_pred(int p1, int p2, const Phenotype &pt, const Phenotype &ptpred, std::vector< std::vector<double> > &ys)
    {
        int m = ys.size();
        for (int j = 0; j < m; ++j) {
            auto y1 = pt.dat[j][p1];
            auto y2 = pt.dat[j][p2];
            auto g1 = ptpred.dat[j][p1];
            auto g2 = ptpred.dat[j][p2];
            auto a = (y1 - g1 + y2 - g2) / 2;
            if (std::isfinite(a)) {
                for (auto &e : ys[j])
                    e += a;
            }
        }
    }

    void adjust_pred(int p1, int p2, int p3, const Phenotype &pt, const Phenotype &ptpred, std::vector< std::vector<double> > &ys)
    {
        int m = ys.size();
        for (int j = 0; j < m; ++j) {
            auto y1 = pt.dat[j][p1];
            auto y2 = pt.dat[j][p2];
            auto y3 = pt.dat[j][p3];
            auto g1 = ptpred.dat[j][p1];
            auto g2 = ptpred.dat[j][p2];
            auto g3 = ptpred.dat[j][p3];
            auto a = (y1 - g1 + y2 - g2 + y3 - g3) / 3;
            if (std::isfinite(a)) {
                for (auto &e : ys[j])
                    e += a;
            }
        }
    }

    void adjust_pred(int p1, int p2, int p3, int p4, const Phenotype &pt, const Phenotype &ptpred, std::vector< std::vector<double> > &ys)
    {
        int m = ys.size();
        for (int j = 0; j < m; ++j) {
            auto y1 = pt.dat[j][p1];
            auto y2 = pt.dat[j][p2];
            auto y3 = pt.dat[j][p3];
            auto y4 = pt.dat[j][p4];
            auto g1 = ptpred.dat[j][p1];
            auto g2 = ptpred.dat[j][p2];
            auto g3 = ptpred.dat[j][p3];
            auto g4 = ptpred.dat[j][p4];
            auto a = (y1 - g1 + y2 - g2 + y3 - g3 + y4 - g4) / 4;
            if (std::isfinite(a)) {
                for (auto &e : ys[j])
                    e += a;
            }
        }
    }

    void cross_pred_2(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm, const Phenotype &pt, const Phenotype &ptpred,
        const std::vector<std::string> &traits, const std::vector< std::vector<double> > &effects)
    {
        int n = gt.ind.size();
        int t = traits.size();

        std::vector< std::vector<allele_t> > geno;
        std::vector< std::vector<double> > ys, ys2;

        std::ofstream ofs(par.out + ".pred");
        if (!ofs) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
            return;
        }

        ofs << "P1\tP2";
        for (int i = 0; i < t; ++i) {
            auto s = traits[i];
            ofs << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs << "\t" << s << ".P" << e;
        }
        ofs << "\n";

        std::ofstream ofs2;

        if (!par.pheno.empty()) {
            ofs2.open(par.out + ".pred.adj");

            if (!ofs2) {
                std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
                return;
            }

            ofs2 << "P1\tP2";
            for (int i = 0; i < t; ++i) {
                auto s = traits[i];
                ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
                for (auto e : pct)
                    ofs2 << "\t" << s << ".P" << e;
            }
            ofs2 << "\n";
        }

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (lm.chr.empty())
                    simulate_indep(par.size, i, j, gt, geno);
                else
                    simulate_linkage(par.size, i, j, gt, lm, geno);

                calc_add(gt, effects, geno, ys);

                if (ofs2.is_open()) {
                    ys2 = ys;
                    adjust_pred(i, j, pt, ptpred, ys2);
                    ofs2 << gt.ind[i] << "\t" << gt.ind[j];
                    for (auto &v : ys2) {
                        ofs2 << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : pct)
                            ofs2 << "\t" << quantile(e, v);
                    }
                    ofs2 << "\n";
                }

                ofs << gt.ind[i] << "\t" << gt.ind[j];
                for (auto &v : ys) {
                    ofs << "\t" << mean(v) << "\t" << sd(v);
                    std::sort(v.begin(), v.end());
                    for (auto e : pct)
                        ofs << "\t" << quantile(e, v);
                }
                ofs << "\n";
            }
        }
    }

    void cross_pred_3(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm, const Phenotype &pt, const Phenotype &ptpred,
        const std::vector<std::string> &traits, const std::vector< std::vector<double> > &effects)
    {
        int n = gt.ind.size();
        int t = traits.size();

        std::vector< std::vector<allele_t> > geno;
        std::vector< std::vector<double> > ys, ys2;

        std::ofstream ofs(par.out + ".pred");
        if (!ofs) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
            return;
        }

        ofs << "P1\tP2\tP3";
        for (int i = 0; i < t; ++i) {
            auto s = traits[i];
            ofs << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs << "\t" << s << ".P" << e;
        }
        ofs << "\n";

        std::ofstream ofs2;

        if (!par.pheno.empty()) {
            ofs2.open(par.out + ".pred.adj");

            if (!ofs2) {
                std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
                return;
            }

            ofs2 << "P1\tP2\tP3";
            for (int i = 0; i < t; ++i) {
                auto s = traits[i];
                ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
                for (auto e : pct)
                    ofs2 << "\t" << s << ".P" << e;
            }
            ofs2 << "\n";
        }

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    simulate_indep(par.size, i, j, k, gt, geno);

                    calc_add(gt, effects, geno, ys);

                    if (ofs2.is_open()) {
                        ys2 = ys;
                        adjust_pred(i, j, k, pt, ptpred, ys2);
                        ofs2 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k];
                        for (auto &v : ys2) {
                            ofs2 << "\t" << mean(v) << "\t" << sd(v);
                            std::sort(v.begin(), v.end());
                            for (auto e : pct)
                                ofs2 << "\t" << quantile(e, v);
                        }
                        ofs2 << "\n";
                    }

                    ofs << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k];
                    for (auto &v : ys) {
                        ofs << "\t" << mean(v) << "\t" << sd(v);
                        std::sort(v.begin(), v.end());
                        for (auto e : pct)
                            ofs << "\t" << quantile(e, v);
                    }
                    ofs << "\n";
                }
            }
        }
    }

    void cross_pred_4(const std::vector<int> &pct, const Genotype &gt, const LinkageMap &lm, const Phenotype &pt, const Phenotype &ptpred,
        const std::vector<std::string> &traits, const std::vector< std::vector<double> > &effects)
    {
        int n = gt.ind.size();
        int t = traits.size();

        std::vector< std::vector<allele_t> > geno;
        std::vector< std::vector<double> > ys, ys2;

        std::ofstream ofs(par.out + ".pred");
        if (!ofs) {
            std::cerr << "ERROR: can't open file: " << par.out << ".pred\n";
            return;
        }

        ofs << "P1\tP2\tP3\tP4";
        for (int i = 0; i < t; ++i) {
            auto s = traits[i];
            ofs << "\t" << s << ".MEAN\t" << s << ".SD";
            for (auto e : pct)
                ofs << "\t" << s << ".P" << e;
        }
        ofs << "\n";

        std::ofstream ofs2;

        if (!par.pheno.empty()) {
            ofs2.open(par.out + ".pred.adj");

            if (!ofs2) {
                std::cerr << "ERROR: can't open file: " << par.out << ".pred.adj\n";
                return;
            }

            ofs2 << "P1\tP2\tP3\tP4";
            for (int i = 0; i < t; ++i) {
                auto s = traits[i];
                ofs2 << "\t" << s << ".MEAN\t" << s << ".SD";
                for (auto e : pct)
                    ofs2 << "\t" << s << ".P" << e;
            }
            ofs2 << "\n";
        }

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = j + 1; k < n; ++k) {
                    for (int l = k + 1; l < n; ++l) {
                        simulate_indep(par.size, i, j, k, l, gt, geno);

                        calc_add(gt, effects, geno, ys);

                        if (ofs2.is_open()) {
                            ys2 = ys;
                            adjust_pred(i, j, k, l, pt, ptpred, ys2);
                            ofs2 << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k] << "\t" << gt.ind[l];
                            for (auto &v : ys2) {
                                ofs2 << "\t" << mean(v) << "\t" << sd(v);
                                std::sort(v.begin(), v.end());
                                for (auto e : pct)
                                    ofs2 << "\t" << quantile(e, v);
                            }
                            ofs2 << "\n";
                        }

                        ofs << gt.ind[i] << "\t" << gt.ind[j] << "\t" << gt.ind[k] << "\t" << gt.ind[l];
                        for (auto &v : ys) {
                            ofs << "\t" << mean(v) << "\t" << sd(v);
                            std::sort(v.begin(), v.end());
                            for (auto e : pct)
                                ofs << "\t" << quantile(e, v);
                        }
                        ofs << "\n";
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
    CmdLine cmd("cross [options]");

    cmd.add("--vcf", "VCF genotype file", "");
    cmd.add("--map", "linkage map file", "");
    cmd.add("--effect", "allele effect file", "");
    cmd.add("--pheno", "phenotype file", "");
    cmd.add("--out", "output file", "cross.out");
    cmd.add("--type", "cross type 2/3/4", "2");
    cmd.add("--size", "sample size", "2000");
    cmd.add("--pct", "sample percentiles", "1,25,50,75,100");

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

    Genotype gt;
    Phenotype pt;
    LinkageMap lm;

    std::vector<int> pct = {1, 25, 50, 75, 100};
    parse_percentile(par.pct, pct);

    std::cerr << "INFO: reading genotype file...\n";
    if (read_vcf(par.vcf, gt) != 0)
        return 1;
    std::cerr << "INFO: " << gt.ind.size() << " individuals, " << gt.loc.size() << " loci\n";

    std::vector<std::string> traits;
    std::vector< std::vector<double> > effects;
    if (make_qtl_allele_matrix(gt, traits, effects) != 0)
        return 1;

    if (!par.pheno.empty()) {
        std::cerr << "INFO: reading phenotype file...\n";
        if (read_pheno(par.pheno, pt) != 0)
            return 1;
        std::cerr << "INFO: " << pt.ind.size() << " individuals, " << pt.phe.size() << " traits\n";

        Phenotype pt2;
        pt2.phe = traits;
        pt2.ind = gt.ind;

        for (auto &phe : pt2.phe) {
            auto j = index(pt.phe, phe);
            if (j == pt.phe.size()) {
                std::cerr << "WARNING: can't find phenotype data for trait: " << phe << "\n";
                pt2.dat.push_back(std::vector<double>(pt2.ind.size(), kNaN));
                continue;
            }

            std::vector<double> v;
            for (auto &ind : pt2.ind) {
                auto k = index(pt.ind, ind);
                if (k == pt.ind.size()) {
                    std::cerr << "WARNING: can't find phenotype data for individual: " << ind << "\n";
                    v.push_back(kNaN);
                }
                else
                    v.push_back(pt.dat[j][k]);
            }
            pt2.dat.push_back(v);
        }

        pt = pt2;
    }

    if (!par.map.empty()) {
        if (make_linkage_map(gt, lm) != 0)
            return 1;
    }

    Phenotype ptpred;
    ptpred.phe = traits;
    ptpred.ind = gt.ind;
    if (gt.ploidy != 1) {
        auto dat = gt.dat;
        int n = gt.ind.size();
        for (auto &v : dat) {
            for (int i = 0; i < n; ++i)
                v[i] = v[i * gt.ploidy];
            v.resize(n);
        }
        calc_add(gt, effects, dat, ptpred.dat);
    }
    else
        calc_add(gt, effects, gt.dat, ptpred.dat);

    if (par.type == 2)
        cross_pred_2(pct, gt, lm, pt, ptpred, traits, effects);
    else if (par.type == 3)
        cross_pred_3(pct, gt, lm, pt, ptpred, traits, effects);
    else if (par.type == 4)
        cross_pred_4(pct, gt, lm, pt, ptpred, traits, effects);

    write_qtl_allele_matrix(gt, traits, effects);

    write_pheno(ptpred, par.out + ".parents");

    return 0;
}
