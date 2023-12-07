// Version 0.3 (2022-11-09)

#ifndef _export_lammps_hpp_
#define _export_lammps_hpp_

#include <string>
#include <vector>

// Structure for ion description
struct Ion
{
    Ion (const std::string &ion_name, double ion_charge) : label (ion_name), charge (ion_charge) {};
    std::string label;
    double charge;
};

// Parameters of Morse potentials
struct PotentialMorse
{
    PotentialMorse (int morse_ion1, int morse_ion2, double morse_D0, double morse_a, double morse_R0, double morse_cutoff) :
            ion1 (morse_ion1), ion2 (morse_ion2), D0 (morse_D0), a (morse_a), R0 (morse_R0), cutoff (morse_cutoff) {};
    int ion1;
    int ion2;
    double D0;
    double a;
    double R0;
    double cutoff;
};

// Parameters of Q*erfc potential
struct PotentialQerfc
{
    PotentialQerfc (int qerfc_ion1, int qerfc_ion2, double qerfc_rho, double qerfc_cutoff) :
            ion1 (qerfc_ion1), ion2 (qerfc_ion2), rho (qerfc_rho), cutoff (qerfc_cutoff) {};
    int ion1;
    int ion2;
    double rho;
    double cutoff;
};

struct Potential3Body
{
    Potential3Body (int p3b_cation, int p3b_anion1, int p3b_anion2, double p3b_k, double p3b_theta0, double p3b_cutoff_a1, double p3b_cutoff_a2) :
            cation_central (p3b_cation), anion1 (p3b_anion1), anion2 (p3b_anion2), k (p3b_k), theta0 (p3b_theta0), cutoff_a1 (p3b_cutoff_a1), cutoff_a2 (p3b_cutoff_a2) {}
    int cation_central;
    int anion1;
    int anion2;
    double k;
    double theta0;
    double cutoff_a1;
    double cutoff_a2;
};

/**
* Generates code for LAMMPS force field input script
* @param filename_core filename to be used with all files produced by this export script
* @param ions array of ions for which the potential has been calculated
* @param morse array of Morse pairwise interactions; numbers of ions are agreement with the order in the list of ions
* @param qerfc array of Qerfc pairwise interactions; numbers of ions are agreement with the order in the list of ions
* @return Array of strings - each to be streamed to seperate files (i.e. [0] LAMMPS input script, [1] 2-body potential table, [2] 3-body potential table)
*/
std::vector<std::string> ExportToLAMMPS (
        const std::string &filename_core,
        const std::vector<Ion> &ions,
        const std::vector<PotentialMorse> &morse,
        const std::vector<PotentialQerfc> &qerfc,
        const std::vector<Potential3Body> &cosine);

#endif
