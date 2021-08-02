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

#ifndef CCOUPLES_H_INCLUDED
#define CCOUPLES_H_INCLUDED

#include <vector>
#include "Cindividus.h"

class Ccouples
{
    public:
    //methods
    Ccouples();
    Ccouples(int g, short int a, short int b, unsigned long K);
    //attributs
    unsigned short int generation;//numero de la generation d'existance, compris entre 0 et GenerationNumber-1
    unsigned short int x;//position sur x du deme sur lequel est le couple"
    unsigned short int y;//position sur y du deme sur lequel est le couple
    std::vector<Cindividus> Spouses;//contient les membres du couples, le [0] est la femelle, le [1] est le male
    unsigned long Key;//clef qui servira entre autre de clef dans la map graveyard
};

#endif // CCOUPLES_H_INCLUDED
