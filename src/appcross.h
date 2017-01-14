#ifndef APPCROSS_H
#define APPCROSS_H

#include <map>
#include "main.h"

class AppCross
{
public:
    int run(int argc, char *argv[]);

private:
    int perform();

    void parse_percentile();

    void load_phenotype();

    void load_genotype();

    void load_linakge_map();

    void load_allele_effect();

    void save_allele_matrix() const;

    void update_chrlen(const vector<string> &chr, const vector<double> &pos);

    bool check_linkage_model() const;

    void make_genetic_map();

    void predict_parent();

    void cross_single_indep();

    void cross_single_linkage();

    void cross_threeway_linkage();

    void cross_threeway_indep();

    void cross_double_indep();

    void cross_double_linkage();

    void adjust_prediction(int p1, int p2, vector< vector<double> > &ys) const;

    void adjust_prediction(int p1, int p2, int p3, vector< vector<double> > &ys) const;

    void adjust_prediction(int p1, int p2, int p3, int p4, vector< vector<double> > &ys) const;

    void simulate_indep(int n, int p1, int p2, vector< vector<allele_t> > &geno) const;

    void simulate_linkage(int n, int p1, int p2, vector< vector<allele_t> > &geno) const;

    void simulate_indep(int n, int p1, int p2, int p3, vector< vector<allele_t> > &geno) const;

    void simulate_linkage(int n, int p1, int p2, int p3, vector< vector<allele_t> > &geno) const;

    void simulate_indep(int n, int p1, int p2, int p3, int p4, vector< vector<allele_t> > &geno) const;

    void simulate_linkage(int n, int p1, int p2, int p3, int p4, vector< vector<allele_t> > &geno) const;

    void calc_genotypic_value(const vector< vector<allele_t> > &geno, vector< vector<double> > &pheno) const;

    void crossover(double len, vector<double> &xo) const;

    void meiosis(const vector<double> &xo, const vector<double> &pos, vector<bool> &gam) const;

    void docross(double len, const vector<double> &pos, const vector<allele_t> &dad,
                 const vector<allele_t> &mom, vector<allele_t> &kid) const;

private:
    struct Params
    {
        string vcf;
        string ped;
        string hmp;
        string geno;
        string map;
        string pheno;
        string effect;
        string out;
        string type;
        string perc;
        int size = 2000;
        bool indep = false;
    };

    struct GenetMap
    {
        vector<string> chr;
        vector<double> len;
        vector< vector<size_t> > idx;
        vector< vector<double> > pos;
    };

    Params m_par;
    Genotype m_gt;
    Phenotype m_pt;
    Phenotype m_pred;
    GenetMap m_gmap;
    std::map<string, double> m_chrlen;
    vector<int> m_perc;
    vector<string> m_trait;
    vector< vector<double> > m_effect;
};

#endif // APPCROSS_H
