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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream> //stringstream


#include "param_input.h"


//variable globale du fichier
bool inputCheckBool=true;// si true, affiche tout ce qui est lu, ligne par ligne, dans le fichier de parametres.


using namespace std;


// utilitaire de comparaison de chaines insensible [`a] la casse
//attention sortie non intuitive! voir des usages prec['e]dents
int cmp_nocase(const string& s, const string& s2)
{
	string::const_iterator p = s.begin();
	string::const_iterator p2 = s2.begin();

	while(p != s.end() && p2 != s2.end())
        {
            if (toupper(*p) != toupper(*p2)) return((toupper(*p)<toupper(*p2)) ? -1 : 1);
            ++p; ++p2;
        }
	return((s2.size()==s.size()) ? 0 : (s2.size()<s.size()) ? -1 : 1);
}
// vire les blancs [`a] droite
void rtrim(string *s)
{
	while ((s->length()>0)  && (s->substr(s->length()-1,1)) == " ")
        {
			s->erase(s->length()-1,s->length());
        }
}

void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

int cmp_nocase_no__(const string& cs, const string& cs2) {

    string s(cs); // makes local copy that one could modify. Aletrnatively, (string& s, string& s2) [no const] would modify its arguments
    string s2(cs2);
    replaceAll(s,"_","");
    replaceAll(s2,"_","");

	return(cmp_nocase(s,s2));
}

/***********************************************************************************************/


int evaluateBool(bool &boolean, string buf) { // safe assignment of value `buf' to 'boolean'
	stringstream strstr(buf);
    string locstring;
    strstr>>locstring;
    if(cmp_nocase(locstring,"")==0 || cmp_nocase(locstring,"T")==0 ||
	   cmp_nocase(locstring,"True")==0 || cmp_nocase(locstring,"Yes")==0 ||
	   cmp_nocase(locstring,"Y")==0)
		boolean=true;
    else if(cmp_nocase(locstring,"F")==0 || cmp_nocase(locstring,"False")==0 ||
			cmp_nocase(locstring,"No")==0 || cmp_nocase(locstring,"N")==0)
		boolean=false;
    else {
        cout<<"(!) Suspicious specification for a boolean: "<< buf <<endl;
		cout<<"(!) Only \"\", \"T\", \"True\", \"Yes\", \"Y\", \"F\", \"False\", \"No\", and \"N\" are allowed"<<endl;
		cout<<"I exit..."<< endl;
		cin.get();
		exit(-1);
    }
	return(0);
}

/*********************************************************/

int seeks_settings_file_name(const string cmdlinefilename,string& settingsfilename) {
// conceived to be executed only if the file has been created.
string buf,var;
string::size_type pos;
ifstream settings(cmdlinefilename.c_str(),ios::in);
    if(!settings.is_open()) {
    	cout << "Unable to open file "<<cmdlinefilename<< endl;
    	cerr << "Unable to open file "<<cmdlinefilename<< endl;
    	cerr << "Possible cause: several processes try to read the same file\n on a client/file server architecture.\n";
    	cerr << "Further execution would ignore the command line option.\n";
    	cerr << "Use the CmdLineFileName setting to avoid this.\n";
    	cerr << "I exit for safety";
    	cin.get();
        exit(-1);
    }
    else
    do {
    		getline(settings,buf);
    		if(buf.length()==0) break;
    		while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
    		pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
    		var=buf.substr(0,pos).c_str();
    		if(cmp_nocase(var,"SettingsFile")==0) {
    			stringstream strstr(buf.substr(pos+1));
    			strstr>>settingsfilename;
                cout<<"Will take settings from file "<<settingsfilename<<endl;
    			settings.close();
    			return 0;
    		}
    } while(!settings.eof());
settings.close();
return 0;
}

/************************************************************/


int read_settings_file(const string filename) {
	string buf,var;
	//stringstring good mostly for string/numbers conversions
	string::size_type pos;
	//int tempo;
	char bidon;
	//double val;
	///RcodeStString.clear();
	int lineCounter=1;// to check whether cmdline is the first argument



	ifstream settings(filename.c_str(),ios::in);
	if(!settings.is_open()) {
		cout << "Unable to open file "<<filename<< endl;
		cerr << "Unable to open file "<<filename<< endl;
		cin.get();
		exit(-1);
	}
	else {
		cout<<"Reading settings file "<<filename<<endl;
		do {
			getline(settings,buf);
			if (inputCheckBool) {cout<<"Read line:\n"<<buf<<endl;}
			if(buf.length()==0) goto nextline;
			while((buf[0]==' ')||(buf[0]=='\t')) {buf.erase(0,1);}//vire les blancs initiaux
			pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
			//		pos=std::min(buf.find('='),buf.length());
			if ((buf[pos])=='=') while (buf[pos-1]==' ') {buf.erase(pos-1,1); pos--;}// vire les blancs avant le =
			var=buf.substr(0,pos).c_str();
			if (inputCheckBool) {cout<<"Parsed: "<<var<<endl;}
			if(var.length()==0) goto nextline;

			if(cmp_nocase(var,"InputCheck")==0) {//tordu
				inputCheckBool=true;
				goto nextline;
			}


			if(cmp_nocase(var,"cmdlinefilename")==0) {
				// meaningful on command line only (and processed at the beginning of main()), but written in <cmdline> text file where it should be ignored
				if (lineCounter>1) {
					cerr<<"(!)(!) Cmdlinefilename setting not first argument. Anything could happen! I exit.";
					cin.get();
					exit(-1);
				}
				goto nextline;
			}

            if(cmp_nocase_no__(var,"Pause")==0) {
    /* Pause determines only two correct contexts for cin.get():
           cerr<< error message + if(cinGetOnError) cin.get() + exit
    and
           cout<< some info + if(pauseGP) cin.get() + execution continues
    cinGetOnError is true at declaration in genepop.cpp but may then be set to false in latin.cpp
    */
                string locstring;
                stringstream strstr(buf.substr(pos+1));
                strstr>>locstring;
                if(cmp_nocase_no__(buf.substr(pos+1),"Final")==0) pauseGP=true;
                if(cmp_nocase_no__(locstring,"OnError")==0) cinGetOnError=true;
                if((cmp_nocase_no__(locstring,"Default")==0)) {cinGetOnError=false; pauseGP=false;}
                goto nextline;}

            if(cmp_nocase(var,"SamplingSeed")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>_ptSamplingSeed;
				goto nextline;}
			if(cmp_nocase(var,"DemeSize")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DemeSize;
				goto nextline;}
			if(cmp_nocase(var,"DemeSamplingRatio")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DemeSamplingRatio;
				goto nextline;}
            if(cmp_nocase(var,"IndMeanSample")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>IndMeanSample;
				goto nextline;}
			if(cmp_nocase(var,"GenerationNumber")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>GenerationNumber;
				goto nextline;}
			if(cmp_nocase(var,"DimX")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DimX;
				goto nextline;}
			if(cmp_nocase(var,"DimY")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DimY;
                goto nextline;}
			if(cmp_nocase(var,"Xlimit")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>Xlimit;
				goto nextline;}
            if(cmp_nocase(var,"HabitatSlideBegin")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>HabitatSlideBegin;
				goto nextline;}
            if(cmp_nocase(var,"HabitatSlideEnd")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>HabitatSlideEnd;
				goto nextline;}
            if(cmp_nocase(var,"HabitatSlideDepth")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>HabitatSlideDepth;
				goto nextline;}
            if(cmp_nocase(var,"Swamping")==0) {
                evaluateBool(Swamping,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"FitnessNormal")==0 || cmp_nocase(var,"FitnessNormal")==0) {
                FitnessNormal.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        FitnessNormal.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"FitnessMaladaptation")==0 || cmp_nocase(var,"FitnessMaladaptation")==0) {
                FitnessMaladaptation.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        FitnessMaladaptation.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"FitnessHybridFemale")==0 || cmp_nocase(var,"FitnessHybridFemale")==0) {
                FitnessHybridFemale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        FitnessHybridFemale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"FitnessHybridMale")==0 || cmp_nocase(var,"FitnessHybridMale")==0) {
                FitnessHybridMale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        FitnessHybridMale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"FitnessMaladaptMt")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>FitnessMt;
				goto nextline;}
            if(cmp_nocase(var,"HybridNb")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>HybridNb;
				goto nextline;}
            if(cmp_nocase(var,"DispMax")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>DispMax;
				goto nextline;}
            if(cmp_nocase(var,"mFemale")==0 || cmp_nocase(var,"mFemale")==0) {
                mFemale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        mFemale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"geomFemale")==0 || cmp_nocase(var,"geomFemale")==0) {
                geomFemale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        geomFemale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"mMale")==0 || cmp_nocase(var,"mMale")==0) {
                mMale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        mMale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"geomMale")==0 || cmp_nocase(var,"geomMale")==0) {
                geomMale.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        geomMale.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"HomogamyAllLoci")==0) {
                evaluateBool(HomogamyAllLoci,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"ChoosyFemale")==0) {
                stringstream strstr(buf.substr(pos+1));
                strstr>>ChoosyFemale;
                goto nextline;}
            if(cmp_nocase(var,"MuRate")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>MuRate;
                goto nextline;}
            if(cmp_nocase(var,"InterRecombiRate")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>InterRecombiRate;
                goto nextline;}
            if(cmp_nocase(var,"IntraRecombiRate")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>IntraRecombiRate;
                goto nextline;}
             if(cmp_nocase(var,"AutLociNumber")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>AutLociNumber;
                goto nextline;}
            if(cmp_nocase(var,"MigRatesCorrection")==0) {
                evaluateBool(MigRatesCorrection,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"LowHybridBound")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>LowHybridBound;
                goto nextline;}
            if(cmp_nocase(var,"HighHybridBound")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>HighHybridBound;
                goto nextline;}
            if(cmp_nocase(var,"WriteIdMatrix")==0) {
                evaluateBool(WriteIdMatrix,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteIdentitiesProba")==0) {
                evaluateBool(WriteIdentitiesProba,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteFstHe")==0) {
                evaluateBool(WriteFstHe,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteGenepopFile")==0) {
                evaluateBool(WriteGenepopFile,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteGenepopIntrog")==0) {
                evaluateBool(WriteGenepopIntrog,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteGenepopOrigin")==0) {
                evaluateBool(WriteGenepopOrigin,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteGenepopAlsoPreContact")==0) {
                evaluateBool(WriteGenepopAlsoPreContact,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteIntrogProfile")==0) {
                evaluateBool(WriteIntrogProfile,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WriteIntrogStats")==0) {
                evaluateBool(WriteIntrogStats,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"WritePeriod")==0) {
				stringstream strstr(buf.substr(pos+1));
				strstr>>WritePeriod;
				goto nextline;}
            if(cmp_nocase(var,"EdgeEffects")==0) {
                evaluateBool(EdgeEffects,buf.substr(pos+1));
                goto nextline;}
            if(cmp_nocase(var,"AcceptRates")==0 || cmp_nocase(var,"AcceptRates")==0) {
                AcceptRates.resize(0); // discards default values
                float value;
                string reste=buf.substr(pos+1);
                stringstream strstr(reste);
                    while (!strstr.eof()) {
                        strstr>>value;
                         while (!strstr.eof() && !(isdigit(bidon=strstr.peek())) && bidon!='.' && bidon!='-') strstr.get();
                        AcceptRates.push_back(value);
                    }
                strstr.clear(); // c'est le truc essentiel pour le r['e]util... snif
                goto nextline;}
            if(cmp_nocase(var,"RunNumber")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>RunNumber;
                goto nextline;}
            if(cmp_nocase(var,"AllopatryLast")==0)  {
                stringstream strstr(buf.substr(pos+1));
                strstr>>AllopatryLast;
                goto nextline;}

			if(cmp_nocase(var,"SettingsFile")==0) { //ignored when read from file ! cf seeks_settings_file_name()
				// meaningful on command line only but written in cmdline.txt which is read _after_ settingsfile
				goto nextline;
			}


			if ( ! (var[0]=='%' || var[0]=='#' || var[0]=='/' )) // exclude comments
				cout<<"(!) Unknown keyword \""<<var<<"\" in file "<<filename<<"."<<endl;
			// ici test on cinGetOnError cannot work since this setting may not yet have been set
		nextline: ;	// the pleasure of sin
			lineCounter++;
		} while(!settings.eof());
	}
	settings.close();
	return 0;
}
