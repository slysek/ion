// (C) Tomasz K. Pietrzak & Viktoriia Vlasenko
// Faculty of Physics, Warsaw University of Technology
// Version 0.3, 2022-11-09

#include "export_lammps.h"
#include <string>
#include <cmath>
#include <fstream>

//DODANY KOD

std::vector<double> Gradient(const std::vector<double>& y, const std::vector<double>& x) {
    size_t n = y.size();
    std::vector<double> gradient(n, 0.0);

    if (x.size() != n || n == 0) {
        return gradient;
    }

    for (size_t i = 1; i < n - 1; ++i) {
        double dx = x[i + 1] - x[i - 1];
        gradient[i] = (y[i + 1] - y[i - 1]) / dx;
    }

    gradient[0] = (y[1] - y[0]) / (x[1] - x[0]);
    gradient[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

    return gradient;
}

std::string GenerateQercfTable(const Ion& i1, const Ion& i2, const double rho, const double cutoff){

    int steps = 990;
    double max = 10;
    double r = 0.1;
    double r0 = 0.1;
    double dr = (10 - r)/steps;

    std::vector<double> lista_qrfc;
    std::vector<double> lista_x;

    double i = r;

    while(i <= max){
        lista_x.push_back(i);
        i += dr;
    }

    for(int j = 0; j < lista_x.size(); j++){
        lista_qrfc.push_back(14.39 * i1.charge * i2.charge * erfc(r/rho)/r);
        r += dr;
    }

    std::vector<double> grad = Gradient(lista_qrfc, lista_x);

    std::string record("QERFC");
    record += "_" + i1.label + "_" + i2.label;
    record += "\n";
    record += "N";
    record += "\t" + std::to_string(steps) + " R";
    record += "\t";
    record += std::to_string(r0) + " " + std::to_string(max);
    record += "\n";
    record += "\n";
    for(int k = 0; k < grad.size(); k++){
        record += std::to_string(k + 1);
        record += "\t";
        record += std::to_string(r0);
        r0 += dr;
        record += "\t";
        record += std::to_string(lista_qrfc[k]);
        record += "\t";
        record += std::to_string(-grad[k]);
        record += "\n";
    }

    return record;
}

//KONIEC DODANEGO KODU

std::vector<std::string> ExportToLAMMPS (
        const std::string &filename_core,
        const std::vector<Ion> &ions,
        const std::vector<PotentialMorse> &morse,
        const std::vector<PotentialQerfc> &qerfc,
        const std::vector<Potential3Body> &cosine)
{

    std::string lammps_script ("# --- BEGIN softBV generated LAMMPS code ---\n\n# IONS DEFINITIONS\n");

    // Atom definition section of LAMMPS input script
    int inum = 0;
    for (const Ion &i : ions)
    {
        ++inum;
        std::string record ("group\t");
        record += i.label;
        record += " type ";
        record += std::to_string (inum);
        record += '\n';
        lammps_script += record;
    }
    lammps_script += '\n';

    for (const Ion &i : ions)
    {
        std::string record ("set group\t");
        record += i.label;
        record += " charge ";
        record += std::to_string (i.charge);
        record += '\n';
        lammps_script += record;
    }
    lammps_script += '\n';

    // Force fields section of LAMMPS input script
    lammps_script += "# FORCE FIELDS\n";
    lammps_script += "pair_style hybrid/overlay table linear 991 morse 10.0\n\n";

    // Morse potential in LAMMPS script
    lammps_script += "# Morse potentials\n";
    for (const PotentialMorse &m : morse)
    {
            std::string record ("pair_coeff ");
            record += std::to_string (m.ion1 + 1) + ' ' + std::to_string (m.ion2 + 1); // dodalem +1 ponieważ pokazywało zle typy jonow
            record += " morse ";
            record += std::to_string (m.D0) + ' ' + std::to_string (m.a) + ' ' + std::to_string (m.R0) + ' ' + std::to_string (m.cutoff);
            record += "\t# ";
            record += ions[m.ion1].label + " - " + ions[m.ion2].label + '\n'; //TUTAJ ZMIENILEM bylo [m.ion1 + 1] i [m.ion2 + 1] usunalem te +1 poniewaz byl jakis wyciek pamieci
            lammps_script += record;
    }

    // Qerfc potentials in LAMMPS script
    lammps_script += "# Q*erfc potentials\n";
    for (const PotentialQerfc &q : qerfc)
    {
        std::string record ("pair_coeff ");
        std::string table_label ("QERFC_");
        table_label += ions [q.ion1].label + '_' + ions [q.ion2].label;
        record += std::to_string (q.ion1 + 1) + ' ' + std::to_string (q.ion2 + 1); // dodalem +1 ponieważ pokazywało zle typy jonow
        record += " table ";

        record += filename_core + ".table " + table_label;
        record += "\t# ";
        record += ions [q.ion1].label + " - " + ions [q.ion2].label + '\n'; //TUTAJ ZMIENILEM bylo [m.ion1 + 1] i [m.ion2 + 1] usunalem te +1 poniewaz byl jakis wyciek pamieci
        lammps_script += record;

        // TODO: (1) generate Q*erfc potential tables

        std::string r = GenerateQercfTable(ions[q.ion1], ions[q.ion2], q.rho, q.cutoff);
        std::string nazwa_pliku = "C:/Users/szymon/CLionProjects/ion/" + filename_core + ".table ";

        std::ofstream plik(nazwa_pliku, std::ios::app);
        plik << r;

    }

    // TODO: (2) Generate three-body potential tables and LAMMPS script


    lammps_script += "\n# --- END softBV generated LAMMPS code ---\n\n";

    // Prepare array with strings containing output files content
    std::vector<std::string> output_files;
    output_files.push_back (lammps_script); // LAMMPS input script - the first element of output array

    return output_files;
}



