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

#ifndef CINDIVIDUS_H_INCLUDED
#define CINDIVIDUS_H_INCLUDED

#include <vector>
#include <map>

class Cindividus
{
    public:
    Cindividus();
    unsigned long Parents;//porte l'identifiant (Key) du couple parental
    std::vector<std::vector<bool> > AdaptationLocus;//doit etre de dim nb autosomes X 2, les alleles sont 0 ou 1
    std::vector<std::vector<int> > Genes;//de dim nb autosomes X 2 homologues + 2 gonosomes + Mt. Contient une clef de map contenant le code ISM de l'allele, ainsi que son habitat d'origine
    int AdaptationMt;//0 or 1, information about the habitat of origin of the mitochondrial haplotype
};

#endif // CINDIVIDUS_H_INCLUDED
