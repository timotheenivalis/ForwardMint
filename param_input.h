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

#ifndef PARAM_INPUT_H_INCLUDED
#define PARAM_INPUT_H_INCLUDED
#include <vector>


//Cdemography demo;
extern unsigned int DemeSize;//demo.DemeSize;
extern double SexeRatio;//demo.SexeRatio; N'intervient pas pour le moment
extern double DemeSamplingRatio;//proportion de demes echantillonnes pour la coalescence
extern double IndMeanSample;//nombre moyen d'individus echantillonnes par demes echantillonnes
extern unsigned long GenerationNumber;//demo.GenerationNumber;
extern unsigned int DimX;//demo.DimX;
extern unsigned int DimY;//demo.DimY;
extern unsigned int Xlimit;
extern unsigned int HabitatSlideBegin; //
extern unsigned int HabitatSlideEnd; // si ==0, pas de changement d'habitat
extern int HabitatSlideDepth;
extern bool Swamping;
extern std::vector<long double> FitnessNormal;
extern std::vector<long double> FitnessMaladaptation;
extern std::vector<long double> FitnessHybridFemale;
extern std::vector<long double> FitnessHybridMale;
extern long double FitnessMt;
extern int HybridNb;
extern int DispMax;
extern std::vector<double> mFemale;//taux de migration sur une dimension pour les femelles
extern std::vector<double> geomFemale;//raison geometrique pour les femelles
extern std::vector<double> mMale;//taux de migration sur une dimension pour les males
extern std::vector<double> geomMale;//raison geometrique pour les males
extern bool HomogamyAllLoci;
extern double ChoosyFemale;
extern double MuRate;
extern double InterRecombiRate;
extern double IntraRecombiRate;
extern int AutLociNumber;
extern std::vector<double> AcceptRates;
extern long AllopatryLast;//temps entre la divergence et le contact secondaire, en nombre de generations
extern unsigned int RunNumber;
extern double LowHybridBound;
extern double HighHybridBound;
extern bool MigRatesCorrection;
extern bool WriteIdMatrix;
extern bool WriteIdentitiesProba;
extern bool WriteFstHe;
extern bool WriteGenepopFile;
extern bool WriteGenepopIntrog;
extern bool WriteGenepopOrigin;
extern bool WriteGenepopAlsoPreContact;
extern bool WriteIntrogProfile;
extern bool WriteIntrogStats;
extern unsigned long WritePeriod;
extern bool EdgeEffects;
extern unsigned long int _ptSamplingSeed;
extern bool pauseGP;
extern bool cinGetOnError;



int cmp_nocase(const std::string& s, const std::string& s2);
void rtrim(std::string *s);
int evaluateBool(bool &boolean, std::string buf);
int seeks_settings_file_name(const std::string cmdlinefilename,std::string& settingsfilename);
int read_settings_file(const std::string filename);





#endif // PARAM_INPUT_H_INCLUDED
