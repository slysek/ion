#include "export_lammps.h"
#include <vector>
#include <iostream>
#include <fstream>

int main() {
    std::string filename_core = "lammps_simulation";

    std::vector<Ion> ions = {Ion("Li", +1.0), Ion("Cl", -1.0), Ion("Cu", +2.0)};

    std::vector<PotentialMorse> morse = {PotentialMorse(0, 1, 0.5, 1.5, 2.0, 3.0), PotentialMorse(1, 2, 0.5, 1.5, 2.0, 3.0)};

    std::vector<PotentialQerfc> qerfc = {PotentialQerfc(0, 0, 1.9650, 0.3), PotentialQerfc(1, 1, 2.250, 0.2)};
    std::vector<Potential3Body> cosine = {};

    auto output = ExportToLAMMPS(filename_core, ions, morse, qerfc, cosine);
    std::string name = filename_core + ".table"; // TUTAJ TRZEBA WPISAC SCIEZKE DO PLIKU .table
    std::fstream table_file(name);
    table_file << output[1];

    for (const auto& script : output) {
        std::cout << script << "\n";
    }

    return 0;
}