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

#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED
#include <vector>
#include <map>
#include "Ccouples.h"
#include "Cdemes.h"
#include "Cnodes.h"
#include "CAlleles.h"

std::vector<std::vector<Cdemes> > Finitialisation(unsigned long& Key);
std::map<long,Cnodes> FfullNodes(int const& GeneType);
int Fncoalescent(std::map<long,Cnodes>& Nodes, double& Ploidy);
int FMutations(std::map<long,Cnodes>& Nodes);
int FAddMutations(std::map<long,Cnodes>& Nodes, int& Ancestor, int& MuNumber);
long double FFactorial(int const& length);
int FAllelesStates(std::map<long,Cnodes>& Nodes,std::map<int,CAlleles>& AllelesGeneG);
int Ftransfer(std::map<long,Cnodes>& Nodes, std::vector<std::vector<Cdemes> >& Demes, int& Genetype, int const& gene);
int FMigrations(std::vector<std::vector<std::vector<double> > >& MigRates);
int FAcceptanceRate(double AcceptanceRate[3][3]);
int FInvasion(unsigned int const& years,std::vector<std::vector<Cdemes> >&NextGeneration, int& MovingLimit,double& SlideCompteur,int& FixedHabitatSlideDepth);
int FFiliation(std::vector<std::vector<Cdemes> >& Demes, unsigned int const& x, unsigned int const&y, unsigned int const& years,unsigned long& Key, std::vector<std::vector<Cdemes> >& NextGeneration, std::vector<std::vector<std::vector<double> > > const& MigRates, double const AcceptanceRate[3][3],int& MovingHybridNb,std::vector<std::map<int,CAlleles> >& Alleles);
double FChoosy(double& choosy,Ccouples& YoungCouples,double const AcceptanceRate[3][3]);
int FTranslateFitness(std::vector<long double>& Fecundity);
int FTranslateMigrationParameters();
long double FFitness(std::vector<std::vector<Cdemes> >const& Demes,unsigned int const& c,unsigned int const& OrigineX, unsigned int const& OrigineY,bool sex);
int FHangover(Ccouples& Parents,Cindividus& Spouse,int& sex,std::vector<std::map<int,CAlleles> >& Alleles);
Ccouples FRecombination(Ccouples& Parents);
int FForwardMutation(Cindividus& Spouse, int& sex,std::vector<std::map<int,CAlleles> >& Alleles);
std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > > FSampling(std::vector<std::vector<Cdemes> > const& Demes,std::vector<std::map<int,CAlleles> > & Alleles);
int FCorrectBounds(int& MovingLimit);
int FProbaID(std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > >& NodesGrid, unsigned int const& RUN, std::vector<std::vector<double> >& RunQIBD );
int FGenepopFile(std::vector<std::map<int,CAlleles> >& Alleles, std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > >& NodesGrid, unsigned int const& RUN, bool const& PreContact);
int FIntrogressionStats(std::vector<std::map<int,CAlleles> > const& Alleles, std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > > const& NodesGrid, unsigned int const& RUN, unsigned long const& years);

#endif // MAIN_H_INCLUDED
