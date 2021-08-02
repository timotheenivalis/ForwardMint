//    This is the source code for ForwardMint version 0.3
//    ForwardMint simulates secondary contacts between two taxa (or species,
//    depending on definitions).
//    The software outputs introgression statistics and profiles for autosomal,
//    gonosomal and mitochondrial markers. Processes considered include
//    spatial structure with isolation by distance, multi-locus local adaptation,
//    reduced hybrid survival, recombination, mitochondrial selection,
//    sex-specific dispersal, sex-specific hybrid survival, asymmetric mating,
//    spatial invasion, swamping.
//
//	  Copyright (C) 2012-2021  Timothée Bonnet - timotheebonnetc@gmail.com
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#ifndef CNODES_H_INCLUDED
#define CNODES_H_INCLUDED
#include <vector>

class Cnodes
{
    public:
    Cnodes();
    int Id;//Identifiant, sert de clef quand on met dans un map
    long double BirthDate;//date d'apparition de la lignee
    int Parent;//identifiant de la lignee parentale
    short int habitat;//habitat 0 ou 1 des individus porteurs a la generation 0. Pour coalescence separee des deux especes avant contact secondaire.
    std::vector<int> Offspring;//identifiants des lignees engendrees
    std::vector<bool> MuState;//sequence de 0 et 1, ISM
    short int Allele;//id de l'etat allellique
};

#endif // CNODES_H_INCLUDED
