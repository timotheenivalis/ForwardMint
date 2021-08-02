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

//pre-processor instructions
#include <iostream>
#include <sstream> //stringstream
#include <fstream> // ofstream
#include <vector>
#include <map>
#include <cmath>//fabs
#include <ctime>//clock
#include <iomanip>//setprecision
#include <cstdlib>//exit

using namespace std;

#include "MersenneTwister.h"
#include "Main.h"
#include "param_input.h"
#include "Ccouples.h"
#include "Cdemes.h"
#include "CAlleles.h"
#include "Cnodes.h"


unsigned int DemeSize=10;//individuals
double DemeSamplingRatio=0.1;//proportion de demes echantillonnes pour la coalescence
double IndMeanSample=1;//nombre moyen d'individus echantillonnes par demes echantillonnes
unsigned long GenerationNumber=5;
unsigned int DimX=10;
unsigned int DimY=10;
unsigned int Xlimit=5;//la transition entre les deux habitats. Habitat 0 pour X<Xlimit, habitat 1 pour X>=Xlimit
unsigned int HabitatSlideBegin=0; //
unsigned int HabitatSlideEnd=0; // si ==0, pas de changement d'habitat
int HabitatSlideDepth=0; //si ==0, pas de chandement d'habitat
bool Swamping=false; // if true, taxon 1 cannot enter habitat 0, meaning that there is a sustained  flow of pure genes 0 into background 1
vector<long double> FitnessNormal;//default=1 thank to FTranslateFitness
vector<long double> FitnessMaladaptation;//default=1 thank to FTranslateFitness
vector<long double> FitnessHybridFemale;//default=1 thank to FTranslateFitness
vector<long double> FitnessHybridMale;//default=1 thank to FTranslateFitness
long double FitnessMt=1.0;//fitness rate associated to the mitochondria
int HybridNb=-1;//il ne se passe rien quand vaut -1. Sinon, indique le nombre d'hybridations autorisees par simul
int DispMax=3;
vector<double> mFemale;
vector<double> geomFemale;
vector<double> mMale;
vector<double> geomMale;
bool HomogamyAllLoci=false;
double ChoosyFemale=0.5;//taux de femelles qui commencent la formation du couple, et choisissent donc leur partenaire selon le genotype.
double MuRate=5e-004;
double InterRecombiRate=0.;
double IntraRecombiRate=0.;
int AutLociNumber=1;
vector<double> AcceptRates;
long AllopatryLast=1000;//temps entre la divergence et le contact secondaire, en nombre de generations
unsigned int RunNumber=1;
double LowHybridBound=-1.;
double HighHybridBound=-1.;
bool MigRatesCorrection=false;
bool WriteIdMatrix=false;
bool WriteIdentitiesProba=false;
bool WriteFstHe=false;
bool WriteGenepopFile=false;
bool WriteGenepopIntrog=false;
bool WriteGenepopOrigin=false;
bool WriteGenepopAlsoPreContact=false;
bool WriteIntrogProfile=false;
bool WriteIntrogStats=false;
unsigned long WritePeriod=0;
bool EdgeEffects=true;
bool pauseGP=false;
bool cinGetOnError=false;
MTRand alea;
unsigned long int _ptSamplingSeed=67144630;

string cmdlinefilename="cmdlineArguments.txt";
string settingsfilename="AllForward.txt";//fichiers d'entree avec valeurs des parametres

//ofstream tmrca("C:/Users/Timothée/Documents/B2E/PodarcisCBGP/CodeAllForward/for_R/tmrca.txt");
//ofstream mutest("C:/Users/Timothée/Documents/B2E/PodarcisCBGP/CodeAllForward/for_R/mutest.txt");

int main(int argc, char *argv[])
{
    /**********************************************************///section lecture du fichier des paramètres
    clock_t start,end,startrun,endrun;
    double temps_ecoule(0.0);
    start=clock();
    if (argc>1)
     {
        // to give inline the name of the file in which command line is written

        string buf(argv[1]);
        string::size_type pos=std::min(buf.find('='),std::min(buf.find('\t'),buf.length()));
        string var=buf.substr(0,pos).c_str();
        if(cmp_nocase(var,"CmdlineFileName")==0) cmdlinefilename=buf.substr(pos+1);
        ofstream cmdline(cmdlinefilename.c_str(),ios::out);
        for (int it=1;it<=argc;it++) cmdline<<argv[it]<<endl;
        cmdline<<endl;
        cmdline.close();
        // seeks optional SettingsFile in cmdline
        seeks_settings_file_name(cmdlinefilename,settingsfilename);
     }
    read_settings_file(settingsfilename); //cf migraine.cpp... READS and SETS...
    if (argc>1) read_settings_file(cmdlinefilename);
    //remove(cmdlinefilename.c_str());// on le supprime si on prevoit de faire plein de trucs dans le meme dossier

    alea.seed(_ptSamplingSeed);

    vector<vector<double> > RunQIBD;//Contiendra les proba d'identite en fonction de la distance pour tous les Run, et permettra de calculer l'evolution de leurs moyennes. Dim1=#run; Dim2=distance
    if (WriteIdentitiesProba==true)
        {
            ofstream Stats("Identities.txt");
            Stats.close();//Ouverture et fermeture pour effacer le fichier avant l'utilisation en fin de programme
        }
    if (WriteIdMatrix==true)
        {
            ofstream Matrix("IdMatrix.txt");
            Matrix.close();//Ouverture et fermeture pour effacer le fichier avant l'utilisation en fin de programme
        }
    if (WriteIntrogProfile==true)
        {
            ofstream IntrogProfile("IntrogProfile.txt");
            IntrogProfile.close();//Ouverture et fermeture pour effacer le fichier avant l'utilisation en fin de programme
        }
    if (WriteIntrogStats==true)
        {
            ofstream IntrogStats("IntrogStats.txt");
            IntrogStats.close();//Ouverture et fermeture pour effacer le fichier avant l'utilisation en fin de programme
        }
    if (WriteFstHe==true)
        {
            ofstream FstHeFile("FstHeFile.txt");
            FstHeFile.close();//Ouverture et fermeture pour effacer le fichier avant l'utilisation en fin de programme
        }
    unsigned long WriteCounter(WritePeriod);

    //Parameter Modifications pre-run
    FTranslateFitness(FitnessNormal);
    FTranslateFitness(FitnessMaladaptation);
    FTranslateFitness(FitnessHybridFemale);
    FTranslateFitness(FitnessHybridMale);
    DemeSize=floor(DemeSize/2);
    FTranslateMigrationParameters();

    vector<vector<vector<double> > > MigRates;//pour contenir les taux de migration et les corriger. La premiere dimension correspond au sexe : 0=femelle, 1=male; the second dimension corresponds to taxon (0 or 1); the third dimension correspond to axial distance
    FMigrations(MigRates);
    for (unsigned int RUN(1);RUN<=RunNumber;RUN++)
        {
            startrun=clock();
            cout<<"Run number "<<RUN<<" begins"<<endl;
            WriteCounter = WritePeriod;
/***********************************************************************************///Initialisation et n-coalescent

            unsigned long Key(1);//sert d'identifiant couple
            vector<map<int,CAlleles> > Alleles;//contient pour tous les loci (sauf adaptation locale) les Identifiants IAM des alleles, leur habitat d'origine et la sequence ISM
            Alleles.resize(AutLociNumber+3);
            vector<vector<Cdemes> > Demes=Finitialisation(Key);//creation des demes, pleins d'individus, repartis dans l'espace
            double Ploidy(2.);
            int GeneType(0);//0 vaut Autosome, 1 Z, 2 W, 3 Mt/ Pour la fonction transfert
            for (int g(0);g<AutLociNumber+3;g++)//partie n-coalescent avant le contact secondaire
                {
                    if(g==AutLociNumber) {Ploidy=1.5;GeneType=1;}
                    if(g==AutLociNumber+1) {Ploidy=1.;GeneType=2;}
                    if(g==AutLociNumber+2) {Ploidy=1.;GeneType=3;}
                    map<long,Cnodes> Nodes=FfullNodes(GeneType);
                    Fncoalescent(Nodes,Ploidy);//determine les etats alleliques de tous les genes de tous les individus a la premiere generation
                    FMutations(Nodes);
                    FAllelesStates(Nodes,Alleles[g]);
                    Ftransfer(Nodes,Demes,GeneType,g);
                }
/***********************************************************************************///Pedigree

            vector<vector<Cdemes> > NextGeneration=Demes;//conteneur de la generation suivante, pour pouvoir transvaser
            //Taking a sample for Genepop before secondary contact
            if(WriteGenepopAlsoPreContact==true)
            {
                vector<vector<vector<vector<vector<int> > > > > NodesGridPre=FSampling(Demes, Alleles);//choisi l'echantillon d'individus genotypes
                FGenepopFile(Alleles,NodesGridPre,RUN,true);
                NodesGridPre.clear();
            }

            vector<vector<vector<vector<vector<int> > > > > NodesGrid; // Initialise the sample object for later time points writing.

            double AcceptanceRate[3][3];
            FAcceptanceRate(AcceptanceRate);//fixe les taux d'acceptation du partenaire en fonction du genotype lors de la formation des couples
            int MovingLimit(Xlimit);
            double SlideCompteur(0);
            int FixedHabitatSlideDepth(HabitatSlideDepth);
            int MovingHybridNb(HybridNb);//nombre d'hybridations autorisees dans ce run.

            if(WritePeriod>0)
            {
                NodesGrid=FSampling(Demes, Alleles);//choisi l'echantillon d'individus genotypes
                FIntrogressionStats(Alleles, NodesGrid, RUN, 0);
            }

            for (unsigned long years(1);years<GenerationNumber;years++)
                {
                    FInvasion(years,NextGeneration,MovingLimit,SlideCompteur,FixedHabitatSlideDepth);
                    for (unsigned int x(0);x<DimX;x++)
                        {
                            for (unsigned int y(0);y<DimY;y++)
                                {

                                    FFiliation(Demes,x,y, years, Key, NextGeneration, MigRates, AcceptanceRate,MovingHybridNb,Alleles);//initialiser une generation de juveniles, qui tireront chacun le couple dont ils sont issus, mettre les adultes au cimetiere

                                }
                        }

                    Demes=NextGeneration;

                   if((WritePeriod>0) && (years<(GenerationNumber-1)))
                    {
                        WriteCounter--;
                        if(WriteCounter==0)
                        {
                            NodesGrid=FSampling(NextGeneration, Alleles);//choisi l'echantillon d'individus genotypes
                            FIntrogressionStats(Alleles, NodesGrid, RUN, years);
                            WriteCounter = WritePeriod;
                        }
                    }
  //  cout<<"\r Generation "<<years<<" completed"<<flush;
                }
/***********************************************************************************///Sampling, calculations and output

            HabitatSlideDepth=FixedHabitatSlideDepth;//pour ecriture des fichiers
            NextGeneration.clear();
            NodesGrid=FSampling(Demes, Alleles);//choisi l'echantillon d'individus genotypes
            FCorrectBounds(MovingLimit);
            FProbaID(NodesGrid, RUN, RunQIBD);
            FGenepopFile(Alleles,NodesGrid,RUN,false);
            FIntrogressionStats(Alleles, NodesGrid, RUN, GenerationNumber);
            endrun=clock();
            temps_ecoule=(double)(endrun-startrun)/CLOCKS_PER_SEC;
            cout<<endl<<"run "<<RUN<<" took "<<temps_ecoule<<endl<<endl;
        }//end for (RUN<RunNumber)
    end=clock();
    temps_ecoule=(double)(end-start)/CLOCKS_PER_SEC;
    cout<<endl<<"temps total="<<temps_ecoule<<endl<<"\a\a\a";
    return 0;
}
/******************************************************************************************************************************************************************************/
//**************************************************************/FONCTIONS SECONDAIRES/***************************************************************************************//
/******************************************************************************************************************************************************************************/
vector<vector<Cdemes> > Finitialisation(unsigned long& Key)//initialisation d'individus
{
    vector<vector<Cdemes> > Demes;
    Demes.resize(DimX);
    vector<Ccouples> *BlockedCouples(0);//passage en pointeur pour gagner du temps
    for (unsigned int x(0);x<DimX;x++)
        {
            Demes[x].resize(DimY);
            for (unsigned int y(0);y<DimY;y++)
                {
                    vector<Ccouples> Couples(DemeSize);//cree un groupe d'individus vivant au meme endroit
                    Demes[x][y].Couples=Couples;
                    BlockedCouples=&Demes[x][y].Couples;
                    for (unsigned long i(0);i<(*BlockedCouples).size();i++)
                        {
                            (*BlockedCouples)[i].generation=0;
                            (*BlockedCouples)[i].Key=Key;
                            Key++;
                            (*BlockedCouples)[i].x=x;
                            (*BlockedCouples)[i].y=y;
							(*BlockedCouples)[i].Spouses.resize(2);
                            for (int j(0);j<2;j++)
                                {
                                    (*BlockedCouples)[i].Spouses[j].AdaptationLocus.resize(AutLociNumber);//selection seulement autosomale
                                    (*BlockedCouples)[i].Spouses[j].Genes.resize(AutLociNumber+3);//pour les autosomes, le Z, le W, la Mt
                                    for(int g(0);g<AutLociNumber;g++)//Autosomal
                                        {
                                            (*BlockedCouples)[i].Spouses[j].Genes[g].resize(2);
                                            (*BlockedCouples)[i].Spouses[j].AdaptationLocus[g].resize(2);
                                        }
                                    (*BlockedCouples)[i].Spouses[j].Genes[AutLociNumber].resize(2);    //Z
                                    for(int g(AutLociNumber+1);g<AutLociNumber+3;g++)// W and Mt
                                        {
                                            (*BlockedCouples)[i].Spouses[j].Genes[g].resize(1);
                                        }
                                }
                        }
                    if (x<Xlimit)//donne des genotypes differents aux individus des demes au dela d'une certaine limite geographique //habitat 0
                        {
                            Demes[x][y].habitat=0;
                            for (unsigned long i(0);i<(*BlockedCouples).size();i++)
                                {
                                    for (int j(0);j<2;j++)//sexe, 0=femela, 1=male
                                        {
                                            for (int g(0);g<AutLociNumber;g++)
                                            for (int k(0);k<2;k++)
                                                {
                                                    (*BlockedCouples)[i].Spouses[j].AdaptationLocus[g][k]=0;
                                                }
                                            (*BlockedCouples)[i].Spouses[j].AdaptationMt=0;
                                        }
                                }
                        }
                    else//habitat 1
                        {
                            Demes[x][y].habitat=1;
                            for (unsigned long i(0); i<(*BlockedCouples).size();i++)
                                {
                                    for (int j(0);j<2;j++)//sexe, 0=femela, 1=male
                                        {
                                            for (int g(0);g<AutLociNumber;g++)
                                            for (int k(0);k<2;k++)
                                                {
                                                    (*BlockedCouples)[i].Spouses[j].AdaptationLocus[g][k]=1;
                                                }
                                            (*BlockedCouples)[i].Spouses[j].AdaptationMt=1;
                                        }
                                }
                        }
                }
        }//end for ( long d(0);d<DemeNumber;d++)
    return Demes;
}//end Finitialisation()

/********************************************************/
map<long,Cnodes> FfullNodes(int const& GeneType)//rempli les listes de lignees de la generation 0.
{
    int A(4);
    if (GeneType==1) A=3;//Z
    if (GeneType==2) A=1;//W
    if (GeneType==3) A=2;//Mt (on en met aussi au mâle)
    long count(1);
    map<long,Cnodes> Nodes;
    for(unsigned int x(0);x<DimX;x++)
        {
            for(unsigned int y(0);y<DimY;y++)
                {
                    for(unsigned int DS(0);DS<DemeSize;DS++)
                        {
                            for(int gene(0);gene<A;gene++)
                                {
                                    Nodes[count].Id=count;
                                    if (x<Xlimit)
                                        {
                                            Nodes[count].habitat=0;
                                            Nodes[count].BirthDate=0.;
                                        }
                                    else
                                        {
                                            Nodes[count].habitat=1;
                                            Nodes[count].BirthDate=0.;
                                        }

                                    count++;
                                }
                        }
                }
        }
    return Nodes;
}// end map<long,Cnodes> FfullNodes(int const& GeneType)

/********************************************************/

int Fncoalescent(map<long,Cnodes>& Nodes, double& Ploidy)//n-coalescent avant le forward
{
    if(Ploidy==1.) Ploidy=0.5;//en fait c'est la taille efficace qui nous interesse plutot que la ploidie //////*****************//////////////???????????????????????????????????????????????????????
    map<long,Cnodes>::iterator it;
    vector<vector<vector<int> > > Present;//liste des lignees presentes a chaque generation. Dim1=habitat; Dim2=generation; Dim3=lignees;
    Present.resize(2);//pour habitat 1 et 2
    vector<vector<int> > Founders;//liste des lignees presentes avant la divergence des taxa. il peut y en avoir deux ou plus.
    Founders.resize(1);
    long c(-1);//identifiant des lignees a creer. On demarre a -1 car le numero de lignee le plus faible vaut 1 pour le moment, et on veut eviter parent==0
    long double NodeLength(0);//duree de la lignee (compte en negatif)
    long double Scale[2];
    for(int j(0);j<2;j++)//boucle sur les habitats
        {
            Scale[j]=Ploidy*2.*DemeSize*DimX*DimY;//Ploidy=1(mt) ou 2(n);La mise a l'echelle vaut N (genes)
        }
    if ((DimX>1) && (Xlimit<DimX))
        {
            Scale[0]*=(double(Xlimit)/double(DimX));//separer les deux habitats s'il y a
            Scale[1]*=((double(DimX-Xlimit))/double(DimX));
        }
    for (unsigned int i(0);i<Present.size();i++)//resize les listes de lignees par dates
        {
            Present[i].resize(1);
        }

    for (it=Nodes.begin();it!=Nodes.end();it++)//boucle de copie et d'initialisation de la generation 0
        {
            Present[(*it).second.habitat][0].push_back((*it).first);
        }
    double long time(0.);

    for (int j(0);j<2;j++)//Boucle sur les habitats
        {
            time=0.;
            long i(1);//compteur d'evenements

            while (time>=(-AllopatryLast) && (Present[j][i-1].size()>1))
                {
                    Present[j].push_back(Present[j][i-1]);//on copie la colonne precedante
                    NodeLength=((Scale[j]*2.0/((Present[j][i].size())*(Present[j][i].size()-1.0)))*log(alea()));//longueur de branche mise a l'echelle de la taille pop, en negatif
                    time+=NodeLength;

                    if (time>=(-AllopatryLast))
                        {
                            int a=floor(alea()*(Present[j][i].size()));//tirage de la premiere lignee descendante
                            if(a==Present[j][i].size()){a=Present[j][i].size()-1;}
                            Nodes.find(Present[j][i][a])->second.Parent=c;//met la nouvelle lignee en parent du premier descendant
                            Nodes[c].Id=c;//creation de la lignee parente
                            Nodes[c].habitat=Nodes.find(Present[j][i][a])->second.habitat;
                            Nodes[c].Offspring.push_back(Present[j][i][a]);//rajoute la premiere lignee descendante au parent
                            Present[j][i].erase(Present[j][i].begin()+a);//efface le premier element tire pour le second tirage
                            a=floor(alea()*(Present[j][i].size()));//tirage de la seconde lignee descendante
                            if(a==Present[j][i].size()){a=Present[j][i].size()-1;}
                            Nodes.find(Present[j][i][a])->second.Parent=c;//met la nouvelle lignee en parent du second descendant
                            Nodes[c].Offspring.push_back(Present[j][i][a]);//rajoute la seconde lignee descendante au parent
                            Present[j][i].erase(Present[j][i].begin()+a);//efface le second element pour la generation suivante
                            Present[j][i].push_back(c);//rajoute la lignee issue de la fusion des deux supprimees
                            Nodes[c].BirthDate=time;
    //if (Ploidy==2){if (Present[j][i].size()==1)            {                tmrca<<Nodes[c].BirthDate<<endl;            }}
                        }
                    c--;//le c n'est pas reinitialise a 0 entre les 2 pop => le premier noeud de la pop 2 a pour identifiant (pop1.size())
                    i++;
                }//end while (Present[j][i-1].size()>1)

        for (unsigned int k(0);k<Present[j][i-1].size();k++)
            {
                Founders[0].push_back(Present[j][i-1][k]);//rajoute les lignees presentes juste apres la divergence a la liste de celle presentes avant la divergence
            }
        }//end for (int j(0);j<2;j++)//Boucle sur les habitats
/******/ /*Maintenant on s'occupe de la coalescence avant la divergence, donc avec les deux taxa reunis*/ /***********/

    Scale[0]=Ploidy*2.0*DemeSize*DimX*DimY;//Additionne les tailles des deux taxa (NB, ne change pas s'il n'y avait qu'une pop)
    time=(-AllopatryLast);
    for (unsigned int i(1);i<Founders[0].size();i++)
        {
            Founders.push_back(Founders[i-1]);
            NodeLength=((Scale[0]*2.0/((Founders[i].size())*(Founders[i].size()-1.0)))*log(alea()));//longueur de branche mise a l'echelle de la taille pop
            time+=NodeLength;
            int a=floor(alea()*(Founders[i].size()));
            if(a==Founders[i].size()){a=Founders[i].size()-1;}
            Nodes.find(Founders[i][a])->second.Parent=c;
            Nodes[c].Id=c;
            Nodes[c].Offspring.push_back(Founders[i][a]);//rajoute la premiere lignee descendante au parent
            Founders[i].erase(Founders[i].begin()+a);//efface le premier element tire pour le second tirage
            a=floor(alea()*(Founders[i].size()));//tirage de la seconde lignee descendante
            if(a==Founders[i].size()){a=Founders[i].size()-1;}
            Nodes.find(Founders[i][a])->second.Parent=c;//met la nouvelle lignee en parent du second descendant
            Nodes[c].Offspring.push_back(Founders[i][a]);//rajoute la seconde lignee descendante au parent
            Founders[i].erase(Founders[i].begin()+a);//efface le second element pour la generation suivante
            Founders[i].push_back(Nodes[c].Id);//rajoute la lignee issue de la fusion des deux supprimees
            Nodes[c].BirthDate=time;
            c--;
        }

    return 0;
}// end int Fncoalescent(map<long,Cnodes>& Nodes, double& Ploidy)

/***************************************************************/
int FMutations(map<long,Cnodes>& Nodes)
{
    int MuNumber(0);//Compteur du nombre de sites mutes depuis le debut
    map<long,Cnodes>::iterator it;
    int OldestName;//contiendra l'identifiants de la plus vieille lignee...
    double OldestDate;//... et sa date
    OldestDate=2;//la BirthDate maximale+1, comme ça gere le cas ou une seule lignee (donc l'arbre se fini a la generation 1)
    for (it=Nodes.begin();it!=Nodes.end();it++)//recherche la plus vieille lignee
        {
            if ((*it).second.BirthDate<OldestDate)//le plus vieux a la plus petite birthdate
                {
                    OldestDate=(*it).second.BirthDate;
                    OldestName=(*it).first;
                }
        }
    FAddMutations(Nodes, OldestName, MuNumber);//Parcours l'arbre depuis le premier MRCA
    for (it=Nodes.begin();it!=Nodes.end();it++)//reparcours les lignees actuelles pour completer l'etat mutationnel
        {
            if((*it).second.BirthDate==0.)
                {
                    for (int i=(*it).second.MuState.size();i<MuNumber;i++)//met des 0 de size a MuNumber
                        {
                            (*it).second.MuState.push_back(0);
                        }
                }
        }
//mutest<<MuNumber<<endl;
return 0;
}// end int FMutations(map<long,Cnode>& Nodes)


/******************************************************************/
int FAddMutations(map<long,Cnodes>& Nodes, int& Oldest, int& MuNumber)
{
    long double NodeLength;
    double random;
    int k(0);
    long double lambda(0.);
    long double cursor(0.);
    int scion;
    for (unsigned int i(0); i<Nodes[Oldest].Offspring.size(); i++)//boucle sur les lignees descendantes, il y en a 2 ou exceptionnellement 0 (une seule lignee)
        {
            scion=Nodes.find(Oldest)->second.Offspring[i];//ID de la lignee descendante
            Nodes.find(scion)->second.MuState=Nodes.find(Oldest)->second.MuState;//Copie des etats mutationnels parents

            for (int n=Nodes.find(scion)->second.MuState.size();n<MuNumber;n++)//complete avec des 0 de size a MuNumber, ie les mutations qui sont intervenus chez les freres et cousins
                {
                    Nodes.find(scion)->second.MuState.push_back(0);
                }

            NodeLength=fabs(Nodes.find(Oldest)->second.BirthDate-Nodes.find(scion)->second.BirthDate);
            random=alea();
            cursor=0.;
            k=-1;//compteur du nombre du mutation sur la lignee focale. Part a -1 pour etre a 0 apres un tour de while
            lambda=MuRate*NodeLength;

            while (cursor<random)
                {
                    k++;
                    //binomial : cursor+=FFactorial(NodeLength)*pow(MuRate,k)*pow(ComplementMuRate,NodeLength-k)/(FFactorial(k)*FFactorial(NodeLength-k));
                    cursor+=pow(lambda,k)*exp(-lambda)/FFactorial(k);//Poisson
                }

            MuNumber+=k;//compte le nombre de sites modifies

            for (int n(0);n<k;n++)//ne tourne que si k>0, ie s'il y a eu au moins une mutation
                {
                    Nodes.find(scion)->second.MuState.push_back(1);//chaque nouveau site mute porte l'allele 1
                }

            if (Nodes.find(scion)->second.BirthDate!=0.)//Rappelle la fonction si la lignee descendante n'est pas de la derniere generation
                {
                    FAddMutations(Nodes,scion,MuNumber);
                }
        }
    return 0;
}//end FAddMutations()

/******************************************************************/
long double FFactorial(int const& length)//calcule (length!)
{
    long double sum;
    if (length==0)
        {
            sum=1;
        }
    else
        {
            int n(length);
            sum=length;
            while (n>2) //on s'arrete a n=3, car n-1=3-1=2, en desssous sum*=1 n'a pas d'effet.
                {
                    sum*=(n-1);
                    n--;
                }
        }
    return sum;
}//end long FFactorial(int const& length)//calcule length!

/******************************************************************/

int FAllelesStates(map<long,Cnodes>& Nodes,map<int,CAlleles>& AllelesGeneG)//converti la sequence ISM en ID allele contenu dans l'attribut Id de la classe Cnodes
{
    map<int,vector<bool> > Alleles;//contiendra en clef l'id d'allele, en attribut la sequence ISM correspondante
    map<long,Cnodes>::iterator it;
    int AlleleId(1);
    bool known(0);//on met a 0 pour le premier tour de boucle, car la seconde for ne tourne pas

    for (it=Nodes.begin();it!=Nodes.end();it++)
        {
            if((*it).second.BirthDate==0.)
                {
                    known=0;
                    for (unsigned int i(1);i<=AllelesGeneG.size();i++)//on parcours tous les alleles deja connus et si on en trouve un correspondant on l'attribue au noeud et on saute l'etape suivante
                        {

                            if ((*it).second.MuState==AllelesGeneG.find(i)->second.MuState)
                                {
                                    (*it).second.Allele=signed(i);
                                    known=1;
                                }
                        }
                    if (known==0)//si aucun allele connu ne correspond, on en cree un nouveau
                        {
                            (*it).second.Allele=AlleleId;
                            AllelesGeneG[(*it).second.Allele].MuState=(*it).second.MuState;
                            AllelesGeneG[(*it).second.Allele].habitat=(*it).second.habitat;
                            AlleleId++;
                        }
                }
        }//end for (it=Nodes.begin();it!=Nodes.end();it++)
    return 0;
}//end FAllelesStates(map<long,Cnode>& Nodes)

/******************************************************************/

int Ftransfer(map<long,Cnodes>& Nodes, vector<vector<Cdemes> >& Demes, int& Genetype, int const& gene)//met les alleles issues du n-coalescent dans la grille de demes
{
    long count(1);
    for(unsigned int x(0);x<DimX;x++)
        {
            for(unsigned int y(0);y<DimY;y++)
                {
                    for(unsigned int DS(0);DS<DemeSize;DS++)
                        {
                            if (Genetype==0)//ie Autosomes
                                {
                                    for(int ind(0);ind<2;ind++)
                                        {
                                            for(int p(0);p<2;p++)
                                                {
                                                    Demes[x][y].Couples[DS].Spouses[ind].Genes[gene][p]=Nodes[count].Allele;
                                                    count++;
                                                }
                                        }
                                }
                            if (Genetype==1) //Z
                                {
                                    Demes[x][y].Couples[DS].Spouses[0].Genes[gene][0]=Nodes[count].Allele;//le Z de la femelle
                                    Demes[x][y].Couples[DS].Spouses[0].Genes[gene][1]=0;//le Z fantome de la femelle
                                    count++;
                                    for(int p(0);p<2;p++)//les Z du male
                                        {
                                            Demes[x][y].Couples[DS].Spouses[1].Genes[gene][p]=Nodes[count].Allele;
                                            count++;
                                        }
                                }
                            if (Genetype==2)//W
                                {
                                    Demes[x][y].Couples[DS].Spouses[0].Genes[gene][0]=Nodes[count].Allele;//la femelle
                                    Demes[x][y].Couples[DS].Spouses[1].Genes[gene][0]=0;//le gene fantome du male
                                    count++;
                                }
                            if (Genetype==3)//Mt
                                {
                                    for (int ind(0);ind<2;ind++)
                                        {
                                            Demes[x][y].Couples[DS].Spouses[ind].Genes[gene][0]=Nodes[count].Allele;
                                            count++;
                                        }
                                }
                        }
                }
        }
    Nodes.clear();
    return 0;
}//end int Ftransfer(map<long,Cnodes>& Nodes, vector<vector<Cdemes> >& Demes, int const& Ploidy, int const& gene)

/******************************************************************/

int FMigrations(vector<vector<vector<double> > >& MigRates)//corrige les taux de migration pour la 2D et normalise selon l'immigration maximale, ie au centre de la grille
 {
    if(DimX==1){DispMax=0;}
    for (unsigned int taxon(0); taxon<2; taxon++)
    {
        if ((mFemale[taxon]>1)||(mMale[taxon]>1))
            {
                cerr<<"Max migration argument ="<<max(mFemale[taxon],mMale[taxon])<<"; value not feasible"<<endl;
                if (cinGetOnError==true) cin.get();
                exit(-1);
            }
        if ((geomFemale[taxon]>=1.0)||(geomMale[taxon]>=1.0))
            {
                cerr << "\nSorry, geometric shape parameters must be below 1." << endl;
                cerr << "If you want an Island model approximation, try g=0.999 for instance" << endl;
                cerr << "I exit" << endl;
                if (cinGetOnError==true) cin.get();
                exit(-1);
            }
    }
    //vector<long double> TEST;
    //TEST.resize(DimX+1);
    ofstream ImigRates("ImigrationRates.txt");
    ImigRates.close();

    double m(0.);
    double geom(0.);

    MigRates.resize(2);//pour femelle et male
    vector<double> *SexTaxonMigRates(0);
    for (unsigned int sex(0);sex<2;sex++)
        {
        MigRates[sex].resize(2);//for the two taxa

            for(unsigned int taxon(0); taxon<2; taxon++)
                {
                if (sex==0) // female
                    {
                        m=mFemale[taxon];
                        geom=geomFemale[taxon];
                    }else{
                        m=mMale[taxon];
                        geom=geomMale[taxon];
                    }
                SexTaxonMigRates=&MigRates[sex][taxon];
                (*SexTaxonMigRates).resize(DispMax+1);
                long double Normalize=(m/2)*(1-geom)/(geom-pow(geom,DispMax+1));//evite une repetition du calcul
                double norm(0);
                (*SexTaxonMigRates)[0]=1.-m;
                norm+=(*SexTaxonMigRates)[0]/2;
                for (int a(1);a<=DispMax;a++)
                    {
                        (*SexTaxonMigRates)[a]=pow(geom,a)*Normalize;
                        norm+=(*SexTaxonMigRates)[a];
                    }
                norm*=2.;//must egal 1

                if (MigRatesCorrection==true) // si on le souhaite, on corrige ces taux de migration pour la 2D, et le taux de migration maximal prend la valeur d'entree m
                    {
                        long double maxcumul(0.);
                        long double tmp;
                        for (int i=0;i<=DispMax;i++)
                            {
                                (*SexTaxonMigRates)[i]/=norm;//INUTILE!!!????
                            }
                        double cumul;
                        if(DimY>1)
                            {
                                for(unsigned int xpos=0;xpos<DimX;xpos++)
                                    {
                                        for(unsigned int ypos=0;ypos<DimY;ypos++)
                                            {
                                                cumul=0;
                                                for(int xxpos(-DispMax);xxpos<=DispMax;xxpos++)
                                                    {
                                                        for(int yypos(-DispMax);yypos<=DispMax;yypos++)
                                                            {
                                                                int X1=int(xpos)+xxpos;//il faut definir un int autre que celui de la boucle car on prend parfois l'oppose, la boucle est alors fausse et potentiellement infinie
                                                                int Y1=int(ypos)+yypos;
                                                                if(EdgeEffects==true)//Reflecting edges
                                                                    {
                                                                        if (X1<0)//bords reflectifs
                                                                            {X1=-(int(xpos)+xxpos);}
                                                                        if (X1>=int(DimX))
                                                                            {X1=2*(int(DimX)-1)-X1;}
                                                                        if (Y1<0)
                                                                            {Y1=-(int(ypos)+yypos);}
                                                                        if (Y1>=int(DimY))
                                                                            {Y1=2*(int(DimY)-1)-Y1;}
                                                                    }
                                                                else//Tore
                                                                    {
                                                                        if (X1<0)
                                                                            {X1=int(DimX)+(int(xpos)+xxpos);}
                                                                        if (X1>=int(DimX))
                                                                            {X1=int(xpos)+xxpos-int(DimX);}
                                                                        if (Y1<0)
                                                                            {Y1=int(DimY)+(int(ypos)+yypos);}
                                                                        if (Y1>=int(DimY))
                                                                            {Y1=int(ypos)+yypos-int(DimY);}
                                                                    }
                                                                if (!( X1==int(xpos) && Y1==int(ypos) ))
                                                                    {
                                                                         cumul+=(*SexTaxonMigRates)[fabs(xxpos)]*(*SexTaxonMigRates)[fabs(yypos)];
                                                                    }
                                                            }//end for(yypos
                                                    }//end for(xxpos
                                                if (cumul>maxcumul)
                                                    {
                                                        maxcumul=cumul;
                                                    }
                                            }//end for(ypos
                                    }//end for(xpos
                            }//end if(DimY>1)
                        else //1D...
                            {
                                for(unsigned int xpos=0;xpos<DimX;xpos++)
                                    {
                                        cumul=0;
                                        for(int xxpos(-DispMax);xxpos<=DispMax;xxpos++)
                                            {
                                                int X1=int(xpos)+xxpos;
                                                if(EdgeEffects==true)//Reflecting edges
                                                    {
                                                        if(X1<0)
                                                            {X1=-(int(xpos)+xxpos);}
                                                        if(X1>=int(DimX))
                                                            {X1=2*(int(DimX)-1)-X1;}
                                                    }
                                                else//Circle...
                                                    {
                                                        if (X1<0)
                                                            {X1=int(DimX)+(int(xpos)+xxpos);}
                                                        if (X1>=int(DimX))
                                                            {X1=int(xpos)+xxpos-int(DimX);}
                                                    }
                                                if(!(X1==int(xpos)))
                                                    {
                                                        cumul+=(*SexTaxonMigRates)[fabs(xxpos)];
                                                    }
                                            }//end for(xpos)
                                        if(cumul>maxcumul)
                                            {
                                                maxcumul=cumul;
                                            }
                                    }
                            }
                        tmp=m/maxcumul;
                        if(DimY>1)//2D
                            {
                                tmp=sqrt(tmp);
                            }
                        vector<vector<long double> > MigCheck;
                        MigCheck.resize(DimX);
                        for (unsigned int i(0);i<MigCheck.size();i++)
                            {
                                MigCheck[i].resize(DimY);
                            }
                        for(int i(0);i<=DispMax;i++)
                            {
                                (*SexTaxonMigRates)[i]*=tmp;
                            }
                            //check
                        double migra0=(*SexTaxonMigRates)[0];//constant, should be the final non-immigration rate
                        if(DimY>1)//2D
                            {
                                    for(unsigned int xpos(0);xpos<DimX;xpos++)
                                    {
                                        for (unsigned int ypos(0);ypos<DimY;ypos++)
                                            {
                                                cumul=0.;
                                                for (int xxpos(-DispMax);xxpos<=DispMax;xxpos++)
                                                    {
                                                        for(int yypos(-DispMax);yypos<=DispMax;yypos++)
                                                            {
                                                                int X1=int(xpos)+xxpos;//il faut definir un int autre que celui de la boucle car on prend parfois l'oppose, la boucle est alors fausse et potentiellement infinie
                                                                int Y1=int(ypos)+yypos;
                                                                if(EdgeEffects==true)//Reflecting edges
                                                                    {
                                                                        if (X1<0)//bords reflectifs
                                                                            {X1=-(int(xpos)+xxpos);}
                                                                        if (X1>=int(DimX))
                                                                            {X1=2*(int(DimX)-1)-X1;}
                                                                        if (Y1<0)
                                                                            {Y1=-(int(ypos)+yypos);}
                                                                        if (Y1>=int(DimY))
                                                                            {Y1=2*(int(DimY)-1)-Y1;}
                                                                    }
                                                                else//Tore...
                                                                    {
                                                                        if (X1<0)
                                                                            {X1=int(DimX)+(int(xpos)+xxpos);}
                                                                        if (X1>=int(DimX))
                                                                            {X1=int(xpos)+xxpos-int(DimX);}
                                                                        if (Y1<0)
                                                                            {Y1=int(DimY)+(int(ypos)+yypos);}
                                                                        if (Y1>=int(DimY))
                                                                            {Y1=int(ypos)+yypos-int(DimY);}
                                                                    }
                                                                if ((X1!=int(xpos))||(Y1!=int(ypos)))
                                                                    {
                                                                        cumul+=(*SexTaxonMigRates)[fabs(xxpos)]*(*SexTaxonMigRates)[fabs(yypos)];
                                                                    }

                                                            }//end for(yypos
                                                    }//end for(xxpos
                                                MigCheck[xpos][ypos]=migra0/(cumul+migra0);
                                            }
                                    }
                            }//end 2D
                        else//1D
                            {
                                int ypos=0;
                                for(unsigned int xpos(0);xpos<DimX;xpos++)
                                    {
                                        cumul=0;
                                        for(int xxpos(-DispMax);xxpos<=DispMax;xxpos++)
                                            {
                                                int X1=int(xpos)+xxpos;
                                                if(EdgeEffects==true)//Reflecting edges
                                                    {
                                                        if(X1<0)
                                                            {X1=-(int(xpos)+xxpos);}
                                                        if(X1>=int(DimX))
                                                            {X1=2*(DimX-1)-X1;}
                                                    }
                                                else//Circle...
                                                    {
                                                        if (X1<0)
                                                            {X1=int(DimX)+(int(xpos)+xxpos);}
                                                        if (X1>=int(DimX))
                                                            {X1=int(xpos)+xxpos-int(DimX);}
                                                    }
                                                if(!(X1==int(xpos)))
                                                    {
                                                        cumul+=(*SexTaxonMigRates)[fabs(xxpos)];
                                                    }
                                            }
                                        MigCheck[xpos][ypos]=migra0/(cumul+migra0);
                                    }//end for(xpos)
                            }//end 1D
                        ofstream ImigRates("ImigrationRates.txt", ios::app);
                        ImigRates<<"Sex="<<sex<<endl<<"x\ty\tImigrationRate"<<endl;
                        for(unsigned int x(0);x<DimX;x++)
                        for(unsigned int y(0);y<DimY;y++)
                            {
                                ImigRates<<x<<"\t"<<y<<"\t"<<MigCheck[x][y]<<endl;

                            }
                    }//end if (MigRatesCorrection==true)
                }//end for (unsigned int taxon(0);taxon<2;taxon++)
        }//end for (unsigned int sex(0);sex<2;sex++)
     return 0;
 }//end FMigration()

 /*******************************************************/

int FAcceptanceRate(double AcceptanceRate[3][3])
{
    if (AcceptRates.size()!=0)
        {
            AcceptanceRate[0][0]=AcceptRates[0];
            AcceptanceRate[0][1]=AcceptRates[1];
            AcceptanceRate[0][2]=AcceptRates[2];
            AcceptanceRate[1][0]=AcceptRates[3];
            AcceptanceRate[1][1]=AcceptRates[4];
            AcceptanceRate[1][2]=AcceptRates[5];
            AcceptanceRate[2][0]=AcceptRates[6];
            AcceptanceRate[2][1]=AcceptRates[7];
            AcceptanceRate[2][2]=AcceptRates[8];
        }
    else
    {
        for (int i(0);i<3;i++)
            {
                for (int j(0);j<3;j++)
                    {
                        AcceptanceRate[i][j]=1;//(1-fabs(i-j)/4);//vaut 1 si les genotypes sont identiques,0.75 si une difference d'allele, 0.5 si deux differences d'alleles
                    }
            }
    }
    return 0;
}//end FAcceptanceRate()
 /*******************************************************/

int FInvasion(unsigned int const& years,vector<vector<Cdemes> >& NextGeneration, int& MovingLimit, double& SlideCompteur, int& FixedHabitatSlideDepth)//realise le glissement d'habitat et permet l'invasion d'une espece sur l'autre
 {
    if((HabitatSlideDepth!=0)&&(MovingLimit>=0)&&(MovingLimit<signed(DimX)))//controle qu'il reste du mouvement a faire et que l'on n'est pas sorti de la grille
        {
              if(years>=HabitatSlideBegin)//inutile de borner par la fin du changement, c'est plus sur de laisser faire le nombre de pas voulu, avec l'intervale correspondant a la duree voulue
                {
                    if(SlideCompteur<1.)//permet d'espacer les evenements de changements
                        {
                            SlideCompteur=fabs((double(HabitatSlideEnd-HabitatSlideBegin)/double(FixedHabitatSlideDepth)));
                            if(HabitatSlideDepth>0)
                                {
                                    for(unsigned int y(0);y<DimY;y++)
                                        {
                                            NextGeneration[MovingLimit][y].habitat=0;
                                        }
                                    HabitatSlideDepth--;//positif, tend vers 0
                                    MovingLimit++;//faut le faire apres si on va vers la droite, car Xlimit est deja dans l'habitat 1

                                }
                            if(HabitatSlideDepth<0)
                                {
                                    MovingLimit--;//faut le faire avant si on va faire la gauche, car Xlimit est encore dans l'habitat 1
                                    for (unsigned int y(0);y<DimY;y++)
                                        {
                                            NextGeneration[MovingLimit][y].habitat=1;
                                        }
                                    HabitatSlideDepth++;//negatif, tend vers 0
                                    MovingLimit--;//faut le faire apres si on va vers la droite, car Xlimit est deja dans l'habitat 1
                                }
                        }
                    else
                        {
                            SlideCompteur--;
                        }
                }
        }
    return 0;
 }//end FInvasion()
 /*******************************************************/

int FFiliation(vector<vector<Cdemes> >& Demes, unsigned int const& x, unsigned int const& y, unsigned int const& years, unsigned long& Key, vector<vector<Cdemes> >& NextGeneration, vector<vector<vector<double> > > const& MigRates, double const AcceptanceRate[3][3],int& MovingHybridNb, vector<map<int,CAlleles> >& Alleles)
{
    //Conteneurs pour les femelles NB:On ne peut pas faire un vector size(2) que l'on push_back, il se remplit n'importe comment
    vector<long double> AddSumF;
    long double SumF(0);
    vector<long double> mxF(2);
    vector<long double> myF(2);
    //Conteneurs pour les males
    vector<long double> AddSumM;
    long double SumM(0);
    vector<long double> mxM(2);
    vector<long double> myM(2);
    //Conteneurs pour les deux
    vector<long double> AddSum;
    vector<int> Addc;
    vector<signed long> Addx;
    vector<signed long> Addy;

    const Cdemes *LockDeme(0); // CHECK IF IT IS A GOOD IDEA TO DEFINE HERE
    const Ccouples *LockCouple(0); // CHECK IF IT IS A GOOD IDEA TO DEFINE HERE


if (DispMax>(int(DimY)/2)&&DimY>1)        {cout<<"ERROR maximal dispersal beyond grid size"<<endl;}

    vector<Ccouples> YoungCouples;

    for (unsigned int i(0);i<DemeSize;i++)
        {
            Ccouples Juv=Ccouples(years, x, y, Key);
            Key++;
			Juv.Spouses.resize(2);
			for (int j(0);j<2;j++)
                {
                    Juv.Spouses[j].AdaptationLocus.resize(AutLociNumber);//seulement autosomal
                    Juv.Spouses[j].Genes.resize(AutLociNumber+3);//pour les autosomes, le Z, le W, la Mt
                    for(int g(0);g<AutLociNumber;g++)//Autosomal
                        {
                            Juv.Spouses[j].Genes[g].resize(2);
                            Juv.Spouses[j].AdaptationLocus[g].resize(2);
                        }
                    Juv.Spouses[j].Genes[AutLociNumber].resize(2);//Z
                    for(int g(AutLociNumber+1);g<AutLociNumber+3;g++)//W and Mt
                        {
                            Juv.Spouses[j].Genes[g].resize(1);
                        }
                }
            YoungCouples.push_back(Juv);
        }

    if(DimY>1)//2D
        {
            for (int X=-DispMax;X<=DispMax;X++)
                {
                    for (int Y=-DispMax;Y<=DispMax;Y++)
                        {
                            int X1=int(x)+X;//il faut definir un int autre que celui de la boucle car on prend parfois l'oppose, la boucle est alors fausse et potentiellement infinie
                            int Y1=int(y)+Y;

                            if (EdgeEffects==true)//Bords reflectifs
                                {
                                    if (X1<0)//bords reflectifs
                                        {X1=-(int(x)+X);}
                                    if (X1>=int(DimX))
                                        {X1=2*(int(DimX)-1)-X1;}
                                    if (Y1<0)
                                        {Y1=-(int(y)+Y);}
                                    if (Y1>=int(DimY))
                                        {Y1=2*(int(DimY)-1)-Y1;}
                                }
                            else//circle or tore
                                {
                                    if (X1<0)
                                        {X1=int(DimX)+(int(x)+X);}
                                    if (X1>=int(DimX))
                                        {X1=int(x)+X-int(DimX);}
                                    if (Y1<0)
                                        {Y1=int(DimY)+(int(y)+Y);}
                                    if (Y1>=int(DimY))
                                        {Y1=int(y)+Y-int(DimY);}
                                }


                            for (unsigned int taxon(0); taxon<2; taxon++)
                                {
                                    mxF[taxon]=MigRates[0][taxon][fabs(X)]; mxM[taxon]=MigRates[1][taxon][fabs(X)];//X et pas X1, c'est le mouvement total dont on veut le taux, pas le resultat du mouvement
                                    myF[taxon]=MigRates[0][taxon][fabs(Y)]; myM[taxon]=MigRates[1][taxon][fabs(Y)];
                                }



                            LockDeme=&Demes[X1][Y1];

                            for (unsigned int c(0);c<DemeSize;c++)
                                {

                                LockCouple=&(*LockDeme).Couples[c];

                                        double fitness=FFitness(Demes,c,X1,Y1,false);//false for female
                                        SumF+=fitness*(mxF[(*LockCouple).Spouses[0].AdaptationLocus[0][0]]+ mxF[(*LockCouple).Spouses[0].AdaptationLocus[0][1]])*(myF[(*LockCouple).Spouses[0].AdaptationLocus[0][0]]+ myF[(*LockCouple).Spouses[0].AdaptationLocus[0][1]])/4;
                                        fitness=FFitness(Demes,c,X1,Y1,true);//true for male
                                        SumM+=fitness*(mxM[(*LockCouple).Spouses[1].AdaptationLocus[0][0]]+ mxM[(*LockCouple).Spouses[1].AdaptationLocus[0][1]])*(myM[(*LockCouple).Spouses[1].AdaptationLocus[0][0]]+ myM[(*LockCouple).Spouses[1].AdaptationLocus[0][1]])/4;

                                        AddSumF.push_back(SumF);
                                        AddSumM.push_back(SumM);
                                        Addc.push_back(c);
                                        Addx.push_back(X1);
                                        Addy.push_back(Y1);
                                }//for (unsigned int c(0);c<DemeSize;c++)
                        }//end for (int Y=-DispMax;Y<=DispMax;Y++)
                }//end for (int X(-DispMax);X<DispMax;X++) Sum
        }//end 2D
    else //1D...
        {
            for (int X=-DispMax;X<=DispMax;X++)
                {
                    int X1=int(x)+X;//il faut definir un int autre que celui de la boucle car on prend parfois l'oppose, la boucle est alors fausse et potentiellement infinie
                    if (EdgeEffects==true)//Bords reflectifs
                        {
                            if (X1<0)//bords reflectifs
                                {X1=-(int(x)+X);}
                            if (X1>=int(DimX))
                                {X1=2*(int(DimX)-1)-X1;}
                        }
                    else//circle
                        {
                            if (X1<0)
                                {X1=int(DimX)+(int(x)+X);}
                            if (X1>=int(DimX))
                                {X1=int(x)+X-int(DimX);}
                        }

                    for (unsigned int taxon(0); taxon<2; taxon++)//there was a mistake with taxon<1 in V0.2
                        {
                            mxF[taxon]=MigRates[0][taxon][fabs(X)]; mxM[taxon]=MigRates[1][taxon][fabs(X)];//X et pas X1, c'est le mouvement total dont on veut le taux, pas le resultat du mouvement
                        }

                    LockDeme=&Demes[X1][0];

                    for (unsigned int c(0);c<DemeSize;c++)
                        {
                            LockCouple=&(*LockDeme).Couples[c];
                            double fitness=FFitness(Demes,c,X1,0,false);//false for female
                            SumF+=fitness*(mxF[(*LockCouple).Spouses[0].AdaptationLocus[0][0]]+ mxF[(*LockCouple).Spouses[0].AdaptationLocus[0][1]])/2;
                            fitness=FFitness(Demes,c,X1,0,true);//true for male
                            SumM+=fitness*(mxM[(*LockCouple).Spouses[1].AdaptationLocus[0][0]]+ mxM[(*LockCouple).Spouses[1].AdaptationLocus[0][1]])/2;
                            AddSumF.push_back(SumF); AddSumM.push_back(SumM);
                            Addc.push_back(c);
                            Addx.push_back(X1);
                            Addy.push_back(0);
                        }
                }
        }//end 1D

/*******///la, on tire un membre du couple dans la distri de proba creee juste avant, puis le second de la meme facon, et on a une certaine proba de le rejeter si les genotypes sont differents. On retire le second jusqu'a ce qu'il soit accepte
    long double randomC;//pour le choix du couple
    double random;//pour le test sur l'acceptation ou non du partenaire, et pour quel sexe tire en premier
    long h;
    long double w;
    Ccouples *BlockedCouple(0);//pointera vers les couples designes par le tirage, pour efficacite d'algorithme et d'ecriture
    for (unsigned int i(0); i<DemeSize;i++)
        {
            random=alea();
            int j;//choix du sexe tire en premier, ie celui qui "decide"
            if (ChoosyFemale<random)
                {
                    j=0;
                    AddSum=AddSumF;
                }
            else
                {
                    j=1;
                    AddSum=AddSumM;
                }
            randomC=alea()*AddSum[AddSum.size()-1];//nombre aleatoire compris entre 0 et Sum maximale
            h=0;
            w=0;
            while (w<randomC)
                {
                    w=AddSum[h];
                    h++;
                }
            h--;//on recule d'un pas pour prendre le dernier element avant le random

            BlockedCouple=&Demes[Addx[h]][Addy[h]].Couples[Addc[h]];//Couple focal en pointeur
            YoungCouples[i].Spouses[j].Parents=(*BlockedCouple).Key;
            FHangover((*BlockedCouple), YoungCouples[i].Spouses[j],j,Alleles);//transmet les genes d'adaptation locale des parents aux juveniles

            int k=fabs(1-j);//second membre du couple
            if (k==0)
                {
                    AddSum=AddSumF;
                }
            else
                {
                    AddSum=AddSumM;
                }
            double choosy;
            do {
                    choosy=0.;
                    randomC=alea()*AddSum[AddSum.size()-1];
                    h=0;
                    w=0.;
                    while (w<randomC)
                        {
                            w=AddSum[h];
                            h++;
                        }

                    h--;//on recule d'un pas pour prendre le dernier element avant le random
                    BlockedCouple=&Demes[Addx[h]][Addy[h]].Couples[Addc[h]];//Couple focal en pointeur
                    YoungCouples[i].Spouses[k].Parents=(*BlockedCouple).Key;

                    FHangover((*BlockedCouple), YoungCouples[i].Spouses[k],k,Alleles);//transmet les genes des parents aux juveniles
                    random=alea();
                    choosy=FChoosy(choosy,YoungCouples[i],AcceptanceRate);

                    if(MovingHybridNb==0)//in the case we start with a positive MovingHybridNb and we have reached the point were no more hybridization is posible.
                       {
                            if((YoungCouples[i].Spouses[0].AdaptationLocus[0]==YoungCouples[i].Spouses[0].AdaptationLocus[1])&&(YoungCouples[i].Spouses[1].AdaptationLocus[0]==YoungCouples[i].Spouses[1].AdaptationLocus[1])) //Test si on a affaire a des individus purs
                                {
                                    if(YoungCouples[i].Spouses[0].AdaptationLocus[0]!=YoungCouples[i].Spouses[1].AdaptationLocus[0])//teste si individus de types differents
                                        {
                                            choosy=0.;//choix impossible, cherche un autre partenaire!
                                        }
                                }
                       }
                }while (choosy<random);

            if(MovingHybridNb>0)
                {
                    if((YoungCouples[i].Spouses[0].AdaptationLocus[0]!=YoungCouples[i].Spouses[1].AdaptationLocus[0])||(YoungCouples[i].Spouses[0].AdaptationLocus[0]!=YoungCouples[i].Spouses[1].AdaptationLocus[1])||(YoungCouples[i].Spouses[0].AdaptationLocus[1]!=YoungCouples[i].Spouses[1].AdaptationLocus[1]))// Pour n'autoriser qu'un nombre limite d'hybridations
                         {
                             MovingHybridNb--;
                             cout<<"Only "<<MovingHybridNb<<" Hybridization left"<<endl;
                         }
                }
        }//end for ( int i(0);i<Landscape[x][y].juv.size();i++)

    NextGeneration[x][y].Couples=YoungCouples;//on transfere YoungCouples dans un plus gros conteneur pour pouvoir les passer a Demes ensuite.
    return 0;
}//end FFiliation()

 /*******************************************************/

double FChoosy(double& choosy,Ccouples& BlockedYoungCouple,double const AcceptanceRate[3][3])//calcul le nombre déterminant si on garde ce partenaire ou pas.
{
    double LocusChoosy;
    if (HomogamyAllLoci==true)// all loci used to compute choosy
        {
                for (int g(0);g<AutLociNumber;g++)
            {
                LocusChoosy=AcceptanceRate[BlockedYoungCouple.Spouses[0].AdaptationLocus[g][0]+BlockedYoungCouple.Spouses[0].AdaptationLocus[g][1]][BlockedYoungCouple.Spouses[1].AdaptationLocus[g][0]+BlockedYoungCouple.Spouses[1].AdaptationLocus[g][1]];// a va piocher une valeur dans le tableau MateChoice, en fonction des gènes des deux membres
                choosy+=LocusChoosy/AutLociNumber;//divided by chromosome number
            }
        }
    else //only the main locus used to compute choosy
        {
            choosy=AcceptanceRate[BlockedYoungCouple.Spouses[0].AdaptationLocus[0][0]+BlockedYoungCouple.Spouses[0].AdaptationLocus[0][1]][BlockedYoungCouple.Spouses[1].AdaptationLocus[0][0]+BlockedYoungCouple.Spouses[1].AdaptationLocus[0][1]];// a va piocher une valeur dans le tableau MateChoice, en fonction des gènes des deux membres
        }
    return choosy;
}//end FChoosy
 /*******************************************************/
 int FTranslateFitness(vector<long double>& Fitness)
 {
    if (Fitness.size()==0)//no input
        {
            for (int i(0);i<AutLociNumber;i++)//FitnessNormal always =1
                {
                    Fitness.push_back(1);
                }
        }
    if (Fitness.size()==1)//no differences between loci
        {
            for (int i(1);i<AutLociNumber;i++)//FitnessNormal always =1
                {
                    Fitness.push_back(Fitness[0]);
                }
        }
    if (Fitness.size()==2)//the first locus is different from the others
        {
            for (int i(2);i<AutLociNumber;i++)//FitnessNormal always =1
                {
                    Fitness.push_back(Fitness[1]);
                }
        }
    if ((Fitness.size()<unsigned(AutLociNumber) && Fitness.size()>2) || (Fitness.size()>unsigned(AutLociNumber)))
        {
            cerr<<"ERROR in FTranslateFitness(): Fitness.size="<<Fitness.size()<<". It should be equal to 0,1,2 or to AutLociNumber ("<<AutLociNumber<<")"<<endl;
            cerr<<"I exit"<<endl;
            if (cinGetOnError)
			cin.get();
            exit(-1);
        }
    //else Fitness.size()==AutLociNumber and need no transformation

    return 0;
 }// end int FTranslateFitness(vector<double>& Fitness)

  int FTranslateMigrationParameters()
  {

    if (mFemale.size()==0) // no input, default value to both taxa
        {
            for (int i(0);i<2;i++)
                {
                    mFemale.push_back(0.1);
                }
        }
    if (mFemale.size()==1) // the two taxa have the same parameter
        {
            mFemale.push_back(mFemale[0]);
        }

    if (geomFemale.size()==0) // no input, default value to both taxa
        {
            for (int i(0);i<2;i++)
                {
                    geomFemale.push_back(0.1);
                }
        }
    if (geomFemale.size()==1) // the two taxa have the same parameter
        {
            geomFemale.push_back(geomFemale[0]);
        }

    if (mMale.size()==0) // no input, default value to both taxa
        {
            for (int i(0);i<2;i++)
                {
                    mMale.push_back(0.1);
                }
        }
    if (mMale.size()==1) // the two taxa have the same parameter
        {
            mMale.push_back(mMale[0]);
        }

    if (geomMale.size()==0) // no input, default value to both taxa
        {
            for (int i(0);i<2;i++)
                {
                    geomMale.push_back(0.1);
                }
        }
    if (geomMale.size()==1) // the two taxa have the same parameter
        {
            geomMale.push_back(geomMale[0]);
        }
    return 0;
  }// end int FTranslateMigrationParameters()


  /*******************************************************/
long double FFitness(vector<vector<Cdemes> >const& Demes,unsigned int const& c,unsigned int const& OrigineX,unsigned int const& OrigineY,bool Sex)//calcule la fecondite moyenne d'un couple selon son genotype et l'habitat
{
    long double fitness(1.);
    long double locusfitness(0.);
    //long double SwampingMultiplier(1.);
    int genotype[2][AutLociNumber][2];
    const Cdemes *LockDeme(0);
    LockDeme=&Demes[OrigineX][OrigineY];
    const Ccouples *LockCouple(0);
    LockCouple=&(*LockDeme).Couples[c];
    for (int j(0);j<2;j++)//boucle sur les individus du couple
        {
            for (int g(0); g<AutLociNumber;g++)
            for (int k(0); k<2;k++)//boucle sur les alleles du locus d'adaptation locale
                {
                    genotype[j][g][k]=(*LockCouple).Spouses[j].AdaptationLocus[g][k];        //les x et y ne sont pas les memes que ceux de la fonction FFiliation, la difference est de a et de b respectivement.
                }
        }
    int Habitat=(*LockDeme).habitat;


    for (int g(0);g<AutLociNumber;g++)
        {
            for (int i(0);i<2;i++)//boucle sur les gametes de l'individu femelle (! ce n'est plus une boucle sur le sexe comme juste au dessus !)
                {
                    for (int j(0);j<2;j++)//boucle sur les gametes de l'individu male
                        {
                            if (genotype[0][g][i]==genotype[1][g][j])
                                if (genotype[0][g][i]==Habitat)
                                    {
                                        locusfitness+=FitnessNormal[g];
                                    }
                                else
                                    {
                                        locusfitness+=FitnessMaladaptation[g];
                                    }
                            else
                                {
                                    if(Sex==false)  locusfitness+=FitnessHybridFemale[g];
                                    else locusfitness+=FitnessHybridMale[g];
                                }
                        }
                }//end for (int i(0);i<2;i++) boucles sur gametes femelles
            locusfitness/=4;
            fitness*=locusfitness;
            locusfitness=0.;
        }
    if (Swamping==true && Habitat==0)// currently swamping only affects the first locus
        {
             for (int i(0);i<2;i++)//boucle sur le sexe
                {
                    for (int j(0);j<2;j++)//boucle sur les gametes
                        {
                            if(genotype[i][0][j]==1)
                                {
                                    fitness = 0.;
                                }
                        }
                }
        }

    if ((*LockCouple).Spouses[0].AdaptationMt==1)// selection on mitochondria is absolute, not habitat dependant
        {
            fitness*=FitnessMt;
        }

    if (fitness>1. || fitness<0.)
        {
            cerr<<"ERROR in fitness calculation: fitness="<<fitness<<". It should be between 0 and 1"<<endl;
            cerr<<"I exit"<<endl;
            if (cinGetOnError)
			cin.get();
            exit(-1);
        }
    return fitness;
}//end FFitness

/***********************************************************/

int FHangover(Ccouples& Parents, Cindividus& Spouse, int& sex,vector<map<int,CAlleles> >& Alleles)//transmet les genes des parents a l'individu focal
{
     Ccouples RecombinatedParents=FRecombination(Parents);//Simulate recombination for this mating without modifying the original chromosomes in case of another mating of this parental couple
     int randomFemale(0);
     if(alea()>=0.5)
     {
        randomFemale=1;
     }
     int randomMale(0);
     if(alea()>=0.5)
     {
        randomMale=1;
     }
    //local adaptation genes and Autosomal genes: All autosomal loci are now perfectly linked on an unique chromosome
    for (int g(0);g<AutLociNumber;g++)
        {
            Spouse.AdaptationLocus[g][0]=RecombinatedParents.Spouses[0].AdaptationLocus[g][randomFemale];
            Spouse.Genes[g][0]=RecombinatedParents.Spouses[0].Genes[g][randomFemale];

            Spouse.AdaptationLocus[g][1]=RecombinatedParents.Spouses[1].AdaptationLocus[g][randomMale];
            Spouse.Genes[g][1]=RecombinatedParents.Spouses[1].Genes[g][randomMale];
        }
     int random(0);
     if(alea()>=0.5)
     {
        random=1;
     }

    Spouse.Genes[AutLociNumber][0]=Parents.Spouses[1].Genes[AutLociNumber][random];//Chromosome Z vient de l'un des  2 Z du pere
    if (sex==0)//female
        {
            Spouse.Genes[AutLociNumber+1][0]=Parents.Spouses[0].Genes[AutLociNumber+1][0];//Chromosome W vient du W de la mere
            Spouse.Genes[AutLociNumber][1]=0;//second Z vide, pour verification
        }
    else//male
        {
            Spouse.Genes[AutLociNumber][1]=Parents.Spouses[0].Genes[AutLociNumber][0];//le second Chromosome Z vient du Z de la mere
            Spouse.Genes[AutLociNumber+1][0]=0;//W vide pour verification
        }
    //Mitochondria
    Spouse.Genes[AutLociNumber+2][0]=Parents.Spouses[0].Genes[AutLociNumber+2][0];
    Spouse.AdaptationMt=Parents.Spouses[0].AdaptationMt;

    FForwardMutation(Spouse,sex,Alleles);
    return 0;
}//end Hangover

/***********************************************************/

Ccouples FRecombination(Ccouples& Parents)
{
    Ccouples RecombinatedParents=Parents;//thus all information is already present even if there is no recombination
    Ccouples StorageRecombination=Parents;//thus all information is already present even if there is no recombination
    double lambda=AutLociNumber*IntraRecombiRate+(AutLociNumber-1)*InterRecombiRate;
    int ChromosomeLength=100*lambda;
    int IntraDist=IntraRecombiRate*100;
    int InterDist=InterRecombiRate*100;

    for (int i(0);i<2;i++)//sex
        {
            Cindividus *LockedStorage(&StorageRecombination.Spouses[i]);//fixed
            Cindividus *LockedRecombinatedParent(&RecombinatedParents.Spouses[i]);//modified
            //LockedParent=&Parents.Spouses[i];
            long double cursor=0.;
            int k=-1;
            long double random=alea();
            while (cursor<random)//number of recombination events
                {
                    k++;
                    cursor+=pow(lambda,k)*exp(-lambda)/FFactorial(k);//Poisson
                    //FFactorial(ChromosomeLength)*pow(RecombiRate,k)*pow(1.-RecombiRate,ChromosomeLength-k)/(FFactorial(k)*FFactorial(ChromosomeLength-k));//binomial
                }

            for (int r(0);r<k;r++)//places of recombinations
                {
                    double where=alea()*ChromosomeLength;
                    if (where>ChromosomeLength) where=ChromosomeLength;
                    int g(0);
                    double place(0);//where we are currently
                    bool selected(0); //oscillate between 0 (neutral marker) and 1 (selected locus)

                    while (place<where)
                        {
                            if (selected==0)//neutral marker
                                {
                                    (*LockedRecombinatedParent).Genes[g][0]=(*LockedStorage).Genes[g][1];
                                    (*LockedRecombinatedParent).Genes[g][1]=(*LockedStorage).Genes[g][0];
                                    place+=IntraDist;
                                }


                            if (selected==1)// selected locus
                                {
                                    (*LockedRecombinatedParent).AdaptationLocus[g][0]=(*LockedStorage).AdaptationLocus[g][1];
                                    (*LockedRecombinatedParent).AdaptationLocus[g][1]=(*LockedStorage).AdaptationLocus[g][0];
                                    place+=InterDist;
                                    g++;//we move to the next pair of marker+selected locus
                                }
                            selected=abs(1-selected);//inversion
                        }
                    (*LockedStorage)=(*LockedRecombinatedParent);// update the reference for the next recombination
                }
        }//end for(int i(0);i<2;i++)//sex

    return RecombinatedParents;
}// end Ccouples FRecombination(Ccouples& Parents)


/***********************************************************/

int FForwardMutation(Cindividus& Spouse, int& sex,vector<map<int,CAlleles> >& Alleles)
{
    long double random=alea();
    long double cursor=0.;
    int k=-1;//compteur du nombre du mutations. Part a -1 pour etre a 0 apres un tour de while
    int NodeLength=2*(AutLociNumber+1)+1;//somme longueur de branche de tous les genes d'un individu sur 1 generation
//long double lambda=MuRate*NodeLength;//taux de mutation multiplie par le nombre de gene ou une mutation est possible

    while (cursor<random)
        {
            k++;
            cursor+=FFactorial(NodeLength)*pow(MuRate,k)*pow(1.-MuRate,NodeLength-k)/(FFactorial(k)*FFactorial(NodeLength-k));//binomial
        }

    int newallele;
    map<int,CAlleles>::iterator it;
    map<int,bool> DejaVu; //will contain a list of the genes already mutated, as on gene can mutate only once by generation
    for (int d(0);d<(2*AutLociNumber+1)+1;d++)//fill in DejaVu
        {
            DejaVu[d]=0;
        }
    int g(0);

    for (int i(0);i<k;i++)//il faut choisir sur quels genes ont lieu ces mutations
        {
            do{
                g=int(alea()*(2*(AutLociNumber+1)+1));
            }while(DejaVu[g]==1);
            DejaVu[g]=1;//Now cannot mutate any more
            int G(0);
            int a (0);
            if (g<2*AutLociNumber)//si c'est un autosome
                {
                    a=g-2*floor(g/2);//homologue 0 ou 1
                    G=floor(g/2);
                }
            if ( ( ((g==2*AutLociNumber+1)||(g==2*AutLociNumber))&&(sex==1) ))//Z male
                {
                    a=g-2*floor(g/2);//homologue 0 ou 1
                    G=AutLociNumber;
                }
            if ((g==2*AutLociNumber)&&(sex==0))//Z female
                {
                    a=0;//pas le choix, il n'y en a qu'un
                    G=AutLociNumber;
                }
            if ((g==2*AutLociNumber+1)&&(sex==0))//W
                {
                    a=0;//pas le choix, il n'y en a qu'un
                    G=AutLociNumber+1;
                }
            if (g>(2*AutLociNumber+1))//mt
                {
                    a=0;//pas le choix, il n'y en a qu'un
                    G=AutLociNumber+2;
                }
            newallele=Alleles[G].size()+1;//le premier allele a pour code 1
            Alleles[G][newallele]=Alleles[G][Spouse.Genes[G][a]];//on copie l'etat d'origine, habitat et ISM
            Alleles[G][newallele].MuState.push_back(1);
            for (it=Alleles[G].begin();it!=Alleles[G].end()--;it++)//on prend pas le dernier element puisque l'on vient de le modifier manuellement. On rajoute un 0 a tous les autres
                {
                    (*it).second.MuState.push_back(0);
                }
            Spouse.Genes[G][a]=newallele;//on met le nouvel etat allelique de cet individu pour ce gene
        }
    return 0;
}//end FForwardMutation

/***********************************************************/
vector<vector<vector<vector<vector<int> > > > > FSampling(vector<vector<Cdemes> > const& Demes,vector<map<int,CAlleles> > & Alleles)
{
    vector<vector<int> > EmptyNode;//premiere dim gene, seconde dim alleles, Quand on les push_back on fait des ind
    EmptyNode.resize(2*AutLociNumber+3);//Autosomes + Mitochondria + W and Z Gonosomes + Adaptation locale(same numberas autosome)
    for (int p(0);p<2*AutLociNumber+3;p++)
        {
            if ((p<=AutLociNumber) || (p>=AutLociNumber+3))//Autosomes + Z Gonosome + Adaptation Locale
                {
                    EmptyNode[p].resize(2);//pour les aut
                        {
                            for(int a(0);a<2;a++)
                                {
                                    EmptyNode[p][a]=0;
                                }
                        }
                }
            else//pour la mt et W
                {
                    EmptyNode[p].resize(1);
                    EmptyNode[p][0]=0;
                }
        }
    vector<vector<vector<vector<vector<int> > > > > NodesGrid;//dimensions X,Y,Individu,Locus,Chromosome
    NodesGrid.resize(DimX);
    for (unsigned x(0);x<DimX;x++)
        {
            NodesGrid[x].resize(DimY);
        }
    int Sample(0);//compte le nb d'ind echantillonnes.
    double IMS=IndMeanSample*IndMeanSample/(IndMeanSample+1/(2*IndMeanSample));//correction pour la transformation du 0 en 1
    vector<vector<int> > IndList;//contient la liste des id de couples pour un x,y donne
    int Captured;
    int sex;
    int ind;
    for(unsigned int x(0);x<DimX;x++)
    for(unsigned int y(0);y<DimY;y++)
        {
            if (DemeSamplingRatio>alea())//Ce deme est-il echantillonne?
                {
                    IndList.resize(2);
                    for (int i(0);i<2;i++)//Dim 2 pour les femelles(0) et les males(1)
                        {
                            IndList[i].resize(DemeSize);
                            for (unsigned int j(0);j<DemeSize;j++)
                                {
                                    IndList[i][j]=j;
                                }
                        }
                    if (IMS>=DemeSize*2)//si on veut tout le monde, on a le droit d'echantilloner tout le monde!
                        {
                            Captured=DemeSize*2;
                        }
                    else
                        {
                            Captured=max(1,int(alea()*2*IMS));//si oui, combien d'individus preleve-t-on? (Parfois 0!) NB : on multiplie par deux pour avoir une distri uniforme entre 0 et 2*IMS, avec pour moyenne IMS
                        }
                    while((Captured>0)&&((IndList[0].size()!=0)||(IndList[1].size()!=0)))
                        {
                            Sample++;
                            if(alea()>=0.5)
                            {
                                sex=0;
                            }else{
                                sex=1;
                            }
                            if (IndList[sex].size()==0)
                                {
                                    sex=fabs(1-sex);
                                }
                            ind=int(alea()*(IndList[sex].size()-1));
                            //if (ind>(signed(IndList[sex].size())-1)) ind=IndList[sex].size()-1;
                            const Cindividus *Blockedind=&Demes[x][y].Couples[IndList[sex][ind]].Spouses[sex]; //Not 100 sure the const is correct
                            for (int p(0);p<AutLociNumber+3;p++)//Gene
                                {
                                    for(unsigned int a(0);a<EmptyNode[p].size();a++)//homologue (1 ou 2)
                                        {
                                            EmptyNode[p][a]=(*Blockedind).Genes[p][a];
                                        }
                                }
                            for (int p(0);p<AutLociNumber;p++)
                                {
                                    for (int a(0);a<2;a++)//locus d'adpatation locale
                                        {
                                            EmptyNode[AutLociNumber+3+p][a]=(*Blockedind).AdaptationLocus[p][a];
                                        }
                                }
                            NodesGrid[x][y].push_back(EmptyNode);
                            IndList[sex].erase(IndList[sex].begin()+ind);
                            Captured--;
                        }
                }
        }
    if(Sample==0)
        {
            cerr<<"Sampling too stringent, no individual sampled."<<endl;
        }

 //vector<map<int,CAlleles> > CurrentAlleles;
//CurrentAlleles.resize(AutLociNumber+3);
    for (int p(0);p<AutLociNumber+3;p++)//inutile de faire les loci d'adaptation. Recycle les numero d'alleles pour genepop
        {
            //int AlleleCount=1;
            //map<int,int> AlleleList;
            for(unsigned int x(0);x<DimX;x++)
            for(unsigned int y(0);y<DimY;y++)
            for(unsigned int i(0);i<NodesGrid[x][y].size();i++)
            for (unsigned int a(0);a<NodesGrid[x][y][i][p].size();a++)
                {
                    if (NodesGrid[x][y][i][p][a]!=0)//crucial to deal with haploid markers
                        {
            // if (AlleleList.count(NodesGrid[x][y][i][p][a])==0)
                    //{
            //CurrentAlleles[p][AlleleCount]=Alleles[p][NodesGrid[x][y][i][p][a]];
            //AlleleList[NodesGrid[x][y][i][p][a]]=Alleles[p].find(NodesGrid[x][y][i][p][a])->first;//AlleleCount;
//                                    NodesGrid[x][y][i][p][a]=AlleleCount;
//                                    AlleleCount++;
//                                }
//                            else
//                                {
                    NodesGrid[x][y][i][p][a]=Alleles[p].find(NodesGrid[x][y][i][p][a])->first;//AlleleList[NodesGrid[x][y][i][p][a]];
                    //}
                        }
                }
            //AlleleList.clear();
        }
    //Alleles=CurrentAlleles;
    return NodesGrid;
}//end FSampling()

/******************************************************************/

int FCorrectBounds(int& MovingLimit)//corrige les bornes de la zone consideree hybride si ces valeurs n'ont pas ete rentrees.
{
    if (LowHybridBound==-1.)
        {
            LowHybridBound=MovingLimit/2.;
        }
    if (HighHybridBound==-1.)
        {
            HighHybridBound=MovingLimit+(DimX-MovingLimit)/2.;
        }
    return 0;
}

/******************************************************************/

int FProbaID(vector<vector<vector<vector<vector<int> > > > >& NodesGrid, unsigned int const& RUN, vector<vector<double> >& RunQIBD )
{
    int step;//contera les pas separant des individus

    //on compare maintenant deux par deux, tous les alleles presents a tous les chromosomes de tous les genes de tous les individus de tous les y de tous les x
    if ((WriteIdentitiesProba==true)||(WriteIdMatrix==true))
        {
            vector<vector<vector<vector<vector<int> > > > > IdMatrix;//Utilise pour mettre en sortie la matrice des proba d'id (tout loci sauf mt) en fonction de x,y focal, et x',y' compare. La derniere dim contient les comptes globaux en 0 et les comptes d'identite en 1
            vector<vector<double> > IdentityIBD;
            vector<vector<double> > CompteurIBD;
            vector<vector<double> > QIBD;
            IdentityIBD.resize(AutLociNumber+3);
            CompteurIBD.resize(AutLociNumber+3);
            QIBD.resize(AutLociNumber+3);
            for (int p(0);p<AutLociNumber+3;p++)
                {
                    IdentityIBD[p].resize(DimX);
                    CompteurIBD[p].resize(DimX);
                    QIBD[p].resize(DimX);
                }
            IdMatrix.resize(DimX);
            for(int x(0);x<signed(DimX);x++)
                {
                   IdMatrix[x].resize(DimY);
                   for(int y(0);y<signed(DimY);y++)
                        {
                            IdMatrix[x][y].resize(DimX);
                            for(int dx(0);dx<signed(DimX);dx++)
                                {
                                    IdMatrix[x][y][dx].resize(DimY);
                                    for(int dy(0);dy<signed(DimY);dy++)
                                        {
                                            IdMatrix[x][y][dx][dy].resize(2);//Faut faire les resize avant d'entrer dans la partie de la boucle conditionnelle a l'echantillonage
                                        }
                                }
                            for(int i(0);i<signed(NodesGrid[x][y].size());i++)//boucle sur les individus echantillonnes en x,y
                                {
                                    for(int dx(0);dx<signed(DimX);dx++)
                                        {
                                            for(int dy(0);dy<signed(DimY);dy++)
                                                {
                                                    for(int di(0);di<signed(NodesGrid[dx][dy].size());di++)//boucle sur les individus echantillonnes en dx,dy
                                                        {
                                                            if (y==dy)//NB, on calcule en axial, ie sur une seule dimension a la fois
                                                                {
                                                                    step=fabs(x-dx);
                                                                    for(int p(0);p<AutLociNumber;p++)
                                                                        {
                                                                            for(int a(0);a<2;a++)//boucle sur les chromosomes du focal
                                                                                {
                                                                                    for(int da(0);da<2;da++)//boucle sur les chromosomes de la cible. Une seule comparaison interne et pas d'auto comparaison
                                                                                        {
                                                                                            //if ((x!=dx)||(y!=dy)||(i!=di)||(a!=da))//on ne fait pas la comparaison d'un gene avec lui meme
                                                                                            if(!((x==dx)&&(i==di)&&(a==da)))
                                                                                                {
                                                                                                    CompteurIBD[p][step]++;
                                                                                                    IdMatrix[x][y][dx][dy][0]++;//on le fait pour les Aut, pas pour la mt
                                                                                                    if(NodesGrid[x][y][i][p][a]==NodesGrid[dx][dy][di][p][da])
                                                                                                        {
                                                                                                            IdentityIBD[p][step]++;
                                                                                                            IdMatrix[x][y][dx][dy][1]++;//on le fait pour les Aut, pas pour la mt
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                    /***********************/ //Gonosomes, c'est un peu complique....
                                                                    if(NodesGrid[x][y][i][AutLociNumber+1][0]==0)//male, pas de lignees a ce locus, W. On va voir en Z
                                                                        {
                                                                            for(int a(0);a<2;a++)
                                                                                {
                                                                                    if(NodesGrid[dx][dy][di][AutLociNumber+1][0]==0)//male
                                                                                        {
                                                                                            for(int da(0);da<2;da++)
                                                                                                {
                                                                                                    if ((x!=dx)||(y!=dy)||(i!=di)||(a!=da))//on ne fait pas la comparaison d'un gene avec lui meme. On le fait pas quand sexes differents
                                                                                                        {
                                                                                                            CompteurIBD[AutLociNumber][step]++;
                                                                                                            if(NodesGrid[x][y][i][AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][da])
                                                                                                                {
                                                                                                                    IdentityIBD[AutLociNumber][step]++;
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                    else//femelle
                                                                                        {
                                                                                            CompteurIBD[AutLociNumber][step]++;
                                                                                            if(NodesGrid[x][y][i][AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][0])
                                                                                                {
                                                                                                    IdentityIBD[AutLociNumber][step]++;
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                    else//femelle
                                                                        {
                                                                            if(NodesGrid[dx][dy][di][AutLociNumber+1][0]==0)//male, comparaison du Z
                                                                                {
                                                                                    for(int da(0);da<2;da++)
                                                                                        {
                                                                                            CompteurIBD[AutLociNumber][step]++;
                                                                                            if(NodesGrid[x][y][i][AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][da])//Z
                                                                                                {
                                                                                                    IdentityIBD[AutLociNumber][step]++;
                                                                                                }
                                                                                        }
                                                                                }
                                                                            else//femelle, comparaison du Z et du W
                                                                                {
                                                                                    CompteurIBD[AutLociNumber][step]++;
                                                                                    if ((x!=dx)||(y!=dy)||(i!=di))//Un seul test necessaire pour savoir si on est sur des individus differents
                                                                                        {
                                                                                            CompteurIBD[AutLociNumber][step]++;
                                                                                            if(NodesGrid[x][y][i][AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][0])//Z
                                                                                                {
                                                                                                    IdentityIBD[AutLociNumber][step]++;
                                                                                                }

                                                                                            CompteurIBD[AutLociNumber+1][step]++;
                                                                                            if(NodesGrid[x][y][i][AutLociNumber+1][1]==NodesGrid[dx][dy][di][AutLociNumber+1][1])//W
                                                                                                {
                                                                                                    IdentityIBD[AutLociNumber+1][step]++;
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }

                                                                    if ((x!=dx)||(y!=dy)||(i!=di))//on ne fait pas la comparaison d'un gene avec lui meme//la on fait pareil sur la mitochondrie
                                                                        {
                                                                            CompteurIBD[AutLociNumber+2][step]++;
                                                                            if(NodesGrid[x][y][i][AutLociNumber+2][0]==NodesGrid[dx][dy][di][AutLociNumber+2][0])
                                                                                {
                                                                                    IdentityIBD[AutLociNumber+2][step]++;
                                                                                }
                                                                        }
                                                                }// end if(y==dy)
                                                            if (x==dx)////////PAS PRET POUR LES GONOSOMES HEIN??
                                                                {
                                                                    step=fabs(y-dy);
                                                                    if (step!=0)
                                                                        {
                                                                            for(int p(0);p<AutLociNumber;p++)
                                                                                {
                                                                                    for(int a(0);a<2;a++)//boucle sur les chromosomes du focal
                                                                                        {
                                                                                            for(int da(0);da<2;da++)//boucle sur les chromosomes de la cible. Pas d'auto comparaison et une seule comparaison interne
                                                                                                {
                                                                                                    if(!((y==dy)&&(i==di)&&(a==da)))
                                                                                                        {
                                                                                                            CompteurIBD[p][step]++;
                                                                                                            IdMatrix[x][y][dx][dy][0]++;//on le fait pour les Aut, pas pour la mt
                                                                                                            if(NodesGrid[x][y][i][p][a]==NodesGrid[dx][dy][di][p][da])
                                                                                                                {
                                                                                                                    IdentityIBD[p][step]++;
                                                                                                                    IdMatrix[x][y][dx][dy][1]++;//on le fait pour les Aut, pas pour la mt
                                                                                                                }
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                            CompteurIBD[AutLociNumber][step]++;
                                                                            if(!((y==dy)&&(i==di)))
                                                                                {
                                                                                    if(NodesGrid[x][y][i][AutLociNumber+2][0]==NodesGrid[dx][dy][di][AutLociNumber+2][0])
                                                                                        {
                                                                                            IdentityIBD[AutLociNumber+2][step]++;
                                                                                        }
                                                                                }
                                                                        }
                                                                }//end if (dx==x)
                                                        }
                                                }
                                        }
                                }
                        }
                }//end for(int x(0);x<signed(DimX);x++)

            //On passe a l'ecriture des Q(distance)
            if (WriteIdentitiesProba==true)
                {
                    vector<double> QSum;//Moyenne des proba d'identite en fonction de la distance, sur tous les loci autosomaux
                    QSum.resize(DimX);
                    for (int p(0);p<=AutLociNumber;p++)
                        {
                            for (int x(0);x<signed(DimX);x++)
                                {
                                    if (CompteurIBD[p][x]!=0)
                                        {
                                            QIBD[p][x]=double(IdentityIBD[p][x]/CompteurIBD[p][x]);
                                        }
                                    else
                                        {
                                            QIBD[p][x]=0;
                                        }
                                    if(p!=AutLociNumber)//on ne prend pas la mitochondrie
                                        {
                                            QSum[x]+=QIBD[p][x];
                                        }
                                    else
                                        {
                                            QSum[x]/=AutLociNumber;//rien a voir avec la mt, mais on profite de la boucle sur x pour faire la moyenne maintenant
                                        }
                                }
                        }//end for (int p(0);p<=AutLociNumber;p++)
                    ofstream Stats("Identities.txt",ios::app);//Fichier pour les proba d'identites !!!Atention, le fichier est efface en debut de programme!!!
                    Stats<<setprecision(5);
                    RunQIBD.push_back(QSum);
                    if(RUN==1)
                        {
                            Stats<<"run\t\t\t";
                            for (unsigned int x(0);x<DimX;x++)
                                {
                                    Stats<<"Q"<<x<<"\t\tmean Q"<<x<<"\t\t\t";//prepare l'en-tete du fichier de proba d'id
                                }
                            Stats<<endl;
                        }
                    Stats<<RUN<<"\t\t";
                    for (unsigned int x(0);x<DimX;x++)
                        {
                            for (int r(RUN-1);r>0;r--)
                                {
                                    QSum[x]+=RunQIBD[r-1][x];
                                }
                            QSum[x]/=RUN;
                            Stats<<fixed<<RunQIBD[RUN-1][x]<<"\t\t"<<QSum[x]<<"\t\t";//ecrit la proba d'id pour ce x et ce run, et la moyenne des proba d'id pour ce x, sur tous les runs deja faits
                        }
                    Stats<<endl;
                    Stats.close();
                }

            //On passe a l'ecriture de la matrice des identites point par point
            if (WriteIdMatrix==true)
                {
                    ofstream Matrix("IdMatrix.txt",ios::app);//fichier pour la matrice des identites. !!!Atention, le fichier est efface en debut de programme!!!
                    Matrix<<setprecision(5);
                    for(unsigned int x(0);x<DimX;x++)
                        {
                            for(unsigned int y(0);y<DimY;y++)
                                {
                                    for(unsigned int dx(0);dx<DimX;dx++)
                                        {
                                            for(unsigned int dy(0);dy<DimY;dy++)
                                                {
                                                    Matrix<<x<<" "<<y<<"\t"<<dx<<" "<<dy<<"\t";
                                                    if(IdMatrix[x][y][dx][dy][0]!=0)
                                                        {
                                                            Matrix<<fixed<<double(IdMatrix[x][y][dx][dy][1])/double(IdMatrix[x][y][dx][dy][0]);
                                                        }
                                                    else
                                                        {
                                                            Matrix<<"nan";
                                                        }
                                                    Matrix<<endl;
                                                }
                                        }
                                }
                        }
                }
        }//end if((WriteIdProba==true)||(Write IdMatrix==true))

        //ENSUITE C'est pour les Fst par habitats. Oblige de le faire a part, car on veut pas les comparaisons en axial

    if (WriteFstHe==true)
        {
            vector<double> Heterozygosities;
            vector<double> IdQ1;
            vector<double> IdQ2;
            vector<vector<double> > CompteurQ1;
            vector<vector<double> > CompteurQ2;
            vector<vector<double> > FstH;
            vector<int> Haplotypes;
            Heterozygosities.resize(AutLociNumber+3);
            IdQ1.resize(AutLociNumber+3);// pour Fsts et 1-Qintra
            CompteurQ1.resize(AutLociNumber+3);
            IdQ2.resize(AutLociNumber+3);
            CompteurQ2.resize(AutLociNumber+3);
            FstH.resize(AutLociNumber+3);
            Haplotypes.resize(AutLociNumber+3);
            for (int p(0);p<AutLociNumber+3;p++)
                {
                    CompteurQ1[p].resize(2);//entre les deux taxa
                    FstH[p].resize(2);//entre les deux taxa)
                    CompteurQ2[p].resize(2);
                    Haplotypes[p]=0;
                }
            vector<vector<int> > *FixedInd;
            for(int x(0);x<signed(DimX);x++)
                {
                    if((x<LowHybridBound)||(x>HighHybridBound))
                        for(int y(0);y<signed(DimY);y++)
                            {
                                for(int i(0);i<signed(NodesGrid[x][y].size());i++)//boucle sur les individus echantillonnes en x,y
                                    {
                                        FixedInd=&NodesGrid[x][y][i];//Gain de temps
                                        //compte des haplotypes
                                        for (int a(0);a<2;a++)
                                            {
                                                for(int p(0);p<AutLociNumber;p++)
                                                    {
                                                        if (Haplotypes[p]<(*FixedInd)[p][a])
                                                            {
                                                                Haplotypes[p]=(*FixedInd)[p][a];
                                                            }
                                                    }
                                                if((Haplotypes[AutLociNumber]<(*FixedInd)[AutLociNumber][a])&&((*FixedInd)[AutLociNumber][a]!=0))//case male
                                                    {
                                                        Haplotypes[AutLociNumber]=NodesGrid[x][y][i][AutLociNumber][a];
                                                    }
                                            }
                                        if((Haplotypes[AutLociNumber+1]<(*FixedInd)[AutLociNumber+1][0])&&((*FixedInd)[AutLociNumber+1][0]!=0))//case female
                                            {
                                                Haplotypes[AutLociNumber+1]=(*FixedInd)[AutLociNumber+1][0];
                                            }
                                        if(Haplotypes[AutLociNumber+2]<(*FixedInd)[AutLociNumber+2][0])//Mt
                                            {
                                                Haplotypes[AutLociNumber+2]=(*FixedInd)[AutLociNumber+2][0];
                                            }
                                        for(int dx(0);dx<signed(DimX);dx++)
                                            {
                                                if((dx<LowHybridBound)||(dx>HighHybridBound))
                                                    {
                                                        for(int dy(0);dy<signed(DimY);dy++)
                                                        for(int di(0);di<signed(NodesGrid[dx][dy].size());di++)//boucle sur les individus echantillonnes en dx,dy
                                                            {
                                                                for(int p(0);p<AutLociNumber;p++)
                                                                    {
                                                                        for(int a(0);a<2;a++)//boucle sur les chromosomes du focal
                                                                        for(int da(0);da<2;da++)//boucle sur les chromosomes de la cible
                                                                            {
                                                                                if ((x!=dx)||(y!=dy)||(i!=di)||(a!=da))//on ne fait pas la comparaison d'un gene avec lui meme. On le fait pas quand sexes differents
                                                                                    {
                                                                                        if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                                            {
                                                                                                CompteurQ1[p][0]++;
                                                                                                if((*FixedInd)[p][a]==NodesGrid[dx][dy][di][p][da])
                                                                                                    {
                                                                                                        CompteurQ1[p][1]++;
                                                                                                    }
                                                                                            }
                                                                                        else
                                                                                            {
                                                                                                CompteurQ2[p][0]++;
                                                                                                if((*FixedInd)[p][a]==NodesGrid[dx][dy][di][p][da])
                                                                                                    {
                                                                                                        CompteurQ2[p][1]++;
                                                                                                    }
                                                                                            }
                                                                                    }
                                                                            }
                                                                    }//end autosomes
                                                            if((*FixedInd)[AutLociNumber+1][0]==0)//male, pas de lignees a ce locus, W. On va voir en Z
                                                                {
                                                                    for(int a(0);a<2;a++)
                                                                        {
                                                                            if(NodesGrid[dx][dy][di][AutLociNumber+1][0]==0)//male
                                                                                {
                                                                                    for(int da(0);da<2;da++)
                                                                                        {
                                                                                             if ((x!=dx)||(y!=dy)||(i!=di)||(a!=da))//on ne fait pas la comparaison d'un gene avec lui meme. On le fait pas quand sexes differents
                                                                                                {
                                                                                                     if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                                                        {
                                                                                                            CompteurQ1[AutLociNumber][0]++;
                                                                                                            if((*FixedInd)[AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][da])
                                                                                                                {
                                                                                                                    CompteurQ1[AutLociNumber][1]++;
                                                                                                                }
                                                                                                        }
                                                                                                    else
                                                                                                        {
                                                                                                            CompteurQ2[AutLociNumber][0]++;
                                                                                                            if((*FixedInd)[AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][da])
                                                                                                                {
                                                                                                                    CompteurQ2[AutLociNumber][1]++;
                                                                                                                }

                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
                                                                            else//di = female
                                                                                {
                                                                                    if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                                        {
                                                                                            CompteurQ1[AutLociNumber][0]++;
                                                                                            if((*FixedInd)[AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][0])
                                                                                                {
                                                                                                    CompteurQ1[AutLociNumber][1]++;
                                                                                                }
                                                                                        }
                                                                                    else
                                                                                        {
                                                                                            CompteurQ2[AutLociNumber][0]++;
                                                                                            if((*FixedInd)[AutLociNumber][a]==NodesGrid[dx][dy][di][AutLociNumber][0])
                                                                                                {
                                                                                                    CompteurQ2[AutLociNumber][1]++;
                                                                                                }

                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                            else//i = femelle
                                                                {
                                                                    if(NodesGrid[dx][dy][di][AutLociNumber+1][0]==0)//male, comparaison du Z
                                                                        {
                                                                            for(int da(0);da<2;da++)
                                                                                {
                                                                                    if(NodesGrid[x][y][i][AutLociNumber][0]!=NodesGrid[dx][dy][di][AutLociNumber][da])//comparaison necessaire car on peut etre sur le meme individu. On le fait pas quand on a des ind de sexe differents
                                                                                        {
                                                                                             if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                                                {
                                                                                                    CompteurQ1[AutLociNumber][0]++;
                                                                                                    if((*FixedInd)[AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][da])
                                                                                                        {
                                                                                                            CompteurQ1[AutLociNumber][1]++;
                                                                                                        }
                                                                                                }
                                                                                            else
                                                                                                {
                                                                                                    CompteurQ2[AutLociNumber][0]++;
                                                                                                    if((*FixedInd)[AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][da])
                                                                                                        {
                                                                                                            CompteurQ2[AutLociNumber][1]++;
                                                                                                        }

                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                    else// di = female, comparaison du Z et du W
                                                                        {
                                                                             if ((x!=dx)||(y!=dy)||(i!=di))//on ne fait pas la comparaison d'un gene avec lui meme. On le fait pas quand sexes differents
                                                                                {
                                                                                   if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                                        {
                                                                                            CompteurQ1[AutLociNumber][0]++;
                                                                                            if((*FixedInd)[AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][0])
                                                                                                {
                                                                                                    CompteurQ1[AutLociNumber][1]++;
                                                                                                }
                                                                                            CompteurQ1[AutLociNumber+1][0]++;
                                                                                            if((*FixedInd)[AutLociNumber+1][0]==NodesGrid[dx][dy][di][AutLociNumber+1][0])
                                                                                                {
                                                                                                    CompteurQ1[AutLociNumber+1][1]++;
                                                                                                }

                                                                                        }
                                                                                    else
                                                                                        {
                                                                                            CompteurQ2[AutLociNumber][0]++;
                                                                                            if((*FixedInd)[AutLociNumber][0]==NodesGrid[dx][dy][di][AutLociNumber][0])
                                                                                                {
                                                                                                    CompteurQ2[AutLociNumber][1]++;
                                                                                                }
                                                                                            CompteurQ2[AutLociNumber+1][0]++;
                                                                                            if((*FixedInd)[AutLociNumber+1][0]==NodesGrid[dx][dy][di][AutLociNumber+1][0])
                                                                                                {
                                                                                                    CompteurQ2[AutLociNumber+1][1]++;
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                 }//end gonosomes
                                                             if ((x!=dx)||(y!=dy)||(i!=di))//on ne fait pas la comparaison d'un gene avec lui meme. On le fait pas quand sexes differents
                                                                {
                                                                    if(((x<LowHybridBound)&&(dx<LowHybridBound))||((x>HighHybridBound)&&(dx>HighHybridBound)))
                                                                        {
                                                                            CompteurQ1[AutLociNumber+2][0]++;
                                                                            if((*FixedInd)[AutLociNumber+2][0]==NodesGrid[dx][dy][di][AutLociNumber+2][0])
                                                                                {
                                                                                    CompteurQ1[AutLociNumber+2][1]++;
                                                                                }
                                                                        }
                                                                    else
                                                                        {
                                                                            CompteurQ2[AutLociNumber+2][0]++;
                                                                            if((*FixedInd)[AutLociNumber+2][0]==NodesGrid[dx][dy][di][AutLociNumber+2][0])
                                                                                {
                                                                                    CompteurQ2[AutLociNumber+2][1]++;
                                                                                }
                                                                        }
                                                                }//end mitochondria
                                                        }
                                                }//end if((dx<LowHybridBound)||dx>HighHybridBound)
                                        }//end for(int i(0);i<signed(NodesGrid[x][y].size());i++)//boucle sur les individus echantillonnes en x,y
                                }
                        }// end if(x<LowHybridBound)||(x>HighHybridBound)
                }//end for(int x(0);x<signed(DimX);x++)

            ofstream FstHeFile("FstHeFile.txt", ios::app);
            FstHeFile<<setprecision(3);
            for (int p(0);p<AutLociNumber+3;p++)
                {
                    IdQ1[p]=CompteurQ1[p][1]/CompteurQ1[p][0];
                    IdQ2[p]=CompteurQ2[p][1]/CompteurQ2[p][0];
                }
            if (RUN==1)//header une seule fois
                {
                    FstHeFile<<"#AllForward output file"<<endl;
                    FstHeFile<<"#Simulation Parameters :"<<endl;
                    FstHeFile<<"#DemeSize="<<DemeSize*2<<"\n#DimX="<<DimX<<"\n#DimY="<<DimY<<"\n#Xlimit="<<Xlimit<<"\n#Generation Number="<<GenerationNumber<<"\n#Allopatry last="<<AllopatryLast<<endl;// DemeSize must be multiplied by 2 because for the user it is the number of ind and for me the number of couple
                    FstHeFile<<"#HabitatSlideBegin"<<HabitatSlideBegin<<"\n#HabitatSlideEnd"<<HabitatSlideEnd<<"\n#HabitatSlideDepth"<<HabitatSlideDepth<<endl;
                    FstHeFile<<"#DemeSamplingRatio="<<DemeSamplingRatio<<"\n#IndMeanSampled="<<IndMeanSample<<endl;
                    FstHeFile<<"#dispmax="<<DispMax<<"\n#mFemale="<<mFemale[0]<< ", "<<mFemale[1] <<"\n#geomFemale="<<geomFemale[0]<<", "<<geomFemale[1]<<"\n#mMale="<<mMale[0]<<", "<<mMale[1]<<"\n#geomMale="<<geomMale[0]<<", "<<geomMale[1]<<"\n#EdgeEffects="<<EdgeEffects<<endl;
                    FstHeFile<<"#FitnessNormal\t"<<"FitnessMaladaptation\t"<<"FitnessHybridFemale\t"<<"FitnessHybridMale\t"<<endl;
                    for (int i(0);i<AutLociNumber;i++)
                        {
                            FstHeFile<<"#\t"<<FitnessNormal[i]<<"\t"<<FitnessMaladaptation[i]<<"\t"<<FitnessHybridFemale[i]<<"\t"<<FitnessHybridMale[i]<<endl;
                        }
                    FstHeFile<<"#Mitochondrial_Maladapation_Fitness="<<FitnessMt<<endl;
                    FstHeFile<<"#HybridNB="<<HybridNb<<"  (-1=infinity)"<<endl;
                    FstHeFile<<"#ChoosyFemale="<<ChoosyFemale<<endl;
                    FstHeFile<<"#AcceptRates=";
                    for (int i(0);i<9;i++)
                        {
                            FstHeFile<<AcceptRates[i]<<" ; ";
                        }
                    FstHeFile<<endl;
                    FstHeFile<<"#MuRate="<<MuRate<<endl;
                    FstHeFile<<"#InterRecombiRate="<<InterRecombiRate<<endl<<"#IntraRecombiRate="<<IntraRecombiRate<<endl;
                    FstHeFile<<"#Stats calculated within "<<0<<" ; "<<LowHybridBound<<" and "<<HighHybridBound<<" ; "<<DimX-1<<endl<<endl;
                    FstHeFile<<"Run\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            FstHeFile<<"FstLoc"<<p<<"\t";
                        }
                    FstHeFile<<"LocZ\t";
                    FstHeFile<<"LocW\t";
                    FstHeFile<<"LocMt\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            FstHeFile<<"HeLoc"<<p<<"\t";
                        }
                    FstHeFile<<"1-QiZ\t";
                    FstHeFile<<"1-QiW\t";
                    FstHeFile<<"1-QiMt\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            FstHeFile<<"Nh"<<p<<"\t";
                        }
                    FstHeFile<<"NhZ\tNhW\tNhMt\t";
                    FstHeFile<<endl;
                }
            FstHeFile<<RUN<<"\t";
            for (int p(0);p<AutLociNumber+3;p++)
                {
                    FstH[p][0]=(IdQ1[p]-IdQ2[p])/(1-IdQ2[p]);
                    FstHeFile<<fixed<<FstH[p][0]<<"\t";
                }
            for (int p(0);p<AutLociNumber+3;p++)// heterozygoties ou 1-Qintra
                {
                    Heterozygosities[p]=(1.0-IdQ1[p]);
                    FstHeFile<<fixed<<Heterozygosities[p]<<"\t";
                }
            for (int p(0);p<AutLociNumber+3;p++)//nombre d'haplotypes
                {
                    FstHeFile<<Haplotypes[p]<<"\t";
                }
            FstHeFile<<endl;
            FstHeFile.close();
        }//end if(writeFstHe==true)


    return 0;
}// end FProbaID(vector<map<int,CAlleles> >& Alleles, vector<vector<vector<vector<vector<int> > > > >& NodesGrid, unsigned int const& RUN, vector<vector<double> >& RunQIBD )


/******************************************************************/
int FGenepopFile(vector<map<int,CAlleles> >& Alleles, vector<vector<vector<vector<vector<int> > > > >& NodesGrid, unsigned int const& RUN, bool const& PreContact)//doit ecrire un fichier traitable par Genepop
{

    string Contact("");
    if ((PreContact==true) && (WriteGenepopAlsoPreContact==true))
        {
          Contact = "PreContact";
        }


    if (WriteGenepopFile==true)
        {
            stringstream stst;
            stst << RUN; //transforme entier en strstr
            string fileName("GenepopFile_");
            fileName+=Contact;
            fileName+=stst.str();
            fileName+=".txt";
            ofstream GPFile(fileName.c_str());
            GPFile<<"This file has been generated by the AllForward software"<<endl;
            for(int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AutLocus"<<p<<endl;
                }
            GPFile<<"ZLocus"<<endl;
            GPFile<<"WLocus"<<endl;
            GPFile<<"MtLocus"<<endl;
            for (int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AdaptLocus"<<p<<endl;
                }
            for(int x(0);x<signed(DimX);x++)
                {
                    for(int y(0);y<signed(DimY);y++)
                        {
                            if(NodesGrid[x][y].size()!=0)
                                {
                                    GPFile<<"pop"<<endl;
                                    for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                        {
                                            GPFile<<x<<" "<<y<<" , ";
                                            /*******************************///AUTOSOMES
                                            for(int p(0);p<AutLociNumber;p++)
                                                {
                                                    for(int a(0);a<2;a++)
                                                        {
                                                            if (NodesGrid[x][y][i][p][a] > 999)
                                                                {
                                                                    cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                    cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                    cerr << "Choose a lower mutation rate or ..." << endl;
                                                                    cerr << "\nand start the program again...." << endl;
                                                                    if (cinGetOnError==true) cin.get();
                                                                    exit(-1);
                                                                }
                                                            if (NodesGrid[x][y][i][p][a]<100)
                                                                {
                                                                    GPFile<<0;
                                                                    if(NodesGrid[x][y][i][p][a]<10)
                                                                        {
                                                                            GPFile<<0;
                                                                        }
                                                                }
                                                            GPFile<<NodesGrid[x][y][i][p][a];//
                                                        }
                                                    GPFile<<" ";
                                                }
                                            /*******************************///Z chromosome
                                            for(int a(0);a<2;a++)//si un seul chromosome 000 = donnee manquante (cf Genepop 4.1.1 tutoriel)
                                                {
                                                    if (NodesGrid[x][y][i][AutLociNumber][a]!=0)//en cas de femelle
                                                        {
                                                            if (NodesGrid[x][y][i][AutLociNumber][a]> 999)
                                                                {
                                                                    cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                    cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                    cerr << "Choose a lower mutation rate or ..." << endl;
                                                                    cerr << "\nand start the program again...." << endl;
                                                                    if (cinGetOnError==true) cin.get();
                                                                    exit(-1);
                                                                }
                                                            if (NodesGrid[x][y][i][AutLociNumber][a]<100)
                                                                {
                                                                    GPFile<<0;
                                                                    if(NodesGrid[x][y][i][AutLociNumber][a]<10)
                                                                        {
                                                                            GPFile<<0;
                                                                        }
                                                                }
                                                            GPFile<<NodesGrid[x][y][i][AutLociNumber][a];
                                                        }
                                                    else GPFile<<0<<0<<0;
                                                }
                                            GPFile<<" ";
                                            /*******************************///W chromosome
                                            if(NodesGrid[x][y][i][AutLociNumber+1][0]!=0)
                                                 {
                                                    if(NodesGrid[x][y][i][AutLociNumber+1][0]>999)
                                                        {
                                                                    cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                    cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                    cerr << "Choose a lower mutation rate or ..." << endl;
                                                                    cerr << "\nand start the program again...." << endl;
                                                                    if (cinGetOnError==true) cin.get();
                                                                    exit(-1);
                                                        }
                                                    if(NodesGrid[x][y][i][AutLociNumber+1][0]<100)
                                                        {
                                                            GPFile<<0;
                                                            if(NodesGrid[x][y][i][AutLociNumber+1][0]<10)
                                                                {
                                                                    GPFile<<0;
                                                                }
                                                        }
                                                    GPFile<<NodesGrid[x][y][i][AutLociNumber+1][0]<<" ";
                                                 }
                                                else GPFile<<0<<0<<0<<" ";
                                            /*******************************///MITOCHONDRIA
                                            if(NodesGrid[x][y][i][AutLociNumber+2][0]>999)
                                                {
                                                            cerr << "An allelic state has been found to be larger than 999," << endl;
                                                            cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                            cerr << "Choose a lower mutation rate or ..." << endl;
                                                            cerr << "\nand start the program again...." << endl;
                                                            if (cinGetOnError==true) cin.get();
                                                            exit(-1);
                                                }
                                            if (NodesGrid[x][y][i][AutLociNumber+2][0]<100)
                                                {
                                                    GPFile<<0;
                                                        if(NodesGrid[x][y][i][AutLociNumber+2][0]<10)
                                                            {
                                                                GPFile<<0;
                                                            }
                                                }
                                            GPFile<<NodesGrid[x][y][i][AutLociNumber+2][0]<<" ";
                                            /******************************///Adaptation Locale
                                            for (int p(0);p<AutLociNumber;p++)
                                                {
                                                    for (int k(0);k<2;k++)
                                                        {
                                                            GPFile<<0<<0<<NodesGrid[x][y][i][AutLociNumber+3+p][k]+1;
                                                        }
                                                    GPFile<<" ";
                                                }
                                            GPFile<<endl;
                                        }//end for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                }
                        }
                }// end for(int x(0);x<signed(DimX);x++)
            GPFile.close();
        }//end of genepop writing
 ///Fichier genepop special, sans la zone de contact hybride, et avec seulement 2 pop
    if (WriteGenepopIntrog==true)
        {
            stringstream stst;
            stst << RUN; //transforme entier en strstr
            string fileName("GenepopIntrog_");
            fileName+=Contact;
            fileName+=stst.str();
            fileName+=".txt";
            ofstream GPFile(fileName.c_str());
            GPFile<<"This file has been generated by the AllForward software"<<endl;
            for(int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AutLocus"<<p<<endl;
                }
            GPFile<<"ZLocus"<<endl;
            GPFile<<"WLocus"<<endl;
            GPFile<<"MtLocus"<<endl;
            for (int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AdaptLocus"<<p<<endl;
                }
            bool FirstPop(false), SecondPop(false);
            for(int x(0);x<signed(DimX);x++)
                {
                    if ((x<LowHybridBound)||(x>HighHybridBound))
                        {
                            if((x<LowHybridBound) && (FirstPop==false))
                                {
                                    GPFile<<"pop"<<endl;
                                    FirstPop=true;
                                }
                            if((x>HighHybridBound) && (SecondPop==false))
                                {
                                    GPFile<<"pop"<<endl;
                                    SecondPop=true;
                                }
                            for(int y(0);y<signed(DimY);y++)
                                {
                                    if(NodesGrid[x][y].size()!=0)
                                        {

                                            for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                                {
                                                    GPFile<<x<<" "<<y<<" , ";
                                                    /*******************************///AUTOSOMES
                                                    for(int p(0);p<AutLociNumber;p++)
                                                        {
                                                            for(int a(0);a<2;a++)
                                                                {
                                                                    if (NodesGrid[x][y][i][p][a] > 999)
                                                                        {
                                                                            cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                            cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                            cerr << "Choose a lower mutation rate or ..." << endl;
                                                                            cerr << "\nand start the program again...." << endl;
                                                                            if (cinGetOnError==true) cin.get();
                                                                            exit(-1);
                                                                        }
                                                                    if (NodesGrid[x][y][i][p][a]<100)
                                                                        {
                                                                            GPFile<<0;
                                                                            if(NodesGrid[x][y][i][p][a]<10)
                                                                                {
                                                                                    GPFile<<0;
                                                                                }
                                                                        }
                                                                    GPFile<<NodesGrid[x][y][i][p][a];//
                                                                }
                                                            GPFile<<" ";
                                                        }
                                                    /*******************************///Z chromosome
                                                    for(int a(0);a<2;a++)//si un seul chromosome 000 = donnee manquante (cf Genepop 4.1.1 tutoriel)
                                                        {
                                                            if (NodesGrid[x][y][i][AutLociNumber][a]!=0)//en cas de femelle
                                                                {
                                                                    if (NodesGrid[x][y][i][AutLociNumber][a] > 999)
                                                                        {
                                                                            cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                            cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                            cerr << "Choose a lower mutation rate or ..." << endl;
                                                                            cerr << "\nand start the program again...." << endl;
                                                                            if (cinGetOnError==true) cin.get();
                                                                            exit(-1);
                                                                        }
                                                                    if (NodesGrid[x][y][i][AutLociNumber][a]<100)
                                                                        {
                                                                            GPFile<<0;
                                                                            if(NodesGrid[x][y][i][AutLociNumber][a]<10)
                                                                                {
                                                                                    GPFile<<0;
                                                                                }
                                                                        }
                                                                    GPFile<<NodesGrid[x][y][i][AutLociNumber][a];
                                                                }
                                                            else GPFile<<0<<0<<0;
                                                        }
                                                    GPFile<<" ";
                                                    /*******************************///W chromosome
                                                    if(NodesGrid[x][y][i][AutLociNumber+1][0]!=0)
                                                         {
                                                            if(NodesGrid[x][y][i][AutLociNumber+1][0]>999)
                                                                {
                                                                            cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                            cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                            cerr << "Choose a lower mutation rate or ..." << endl;
                                                                            cerr << "\nand start the program again...." << endl;
                                                                            if (cinGetOnError==true) cin.get();
                                                                            exit(-1);
                                                                }
                                                            if(NodesGrid[x][y][i][AutLociNumber+1][0]<100)
                                                                {
                                                                    GPFile<<0;
                                                                    if(NodesGrid[x][y][i][AutLociNumber+1][0]<10)
                                                                        {
                                                                            GPFile<<0;
                                                                        }
                                                                }
                                                            GPFile<<NodesGrid[x][y][i][AutLociNumber+1][0]<<" ";
                                                         }
                                                        else GPFile<<0<<0<<0<<" ";
                                                    /*******************************///MITOCHONDRIA
                                                    if(NodesGrid[x][y][i][AutLociNumber+2][0]>999)
                                                        {
                                                                    cerr << "An allelic state has been found to be larger than 999," << endl;
                                                                    cerr << "Genepop files do not accept alleles with more than 3 digits" << endl;
                                                                    cerr << "Choose a lower mutation rate or ..." << endl;
                                                                    cerr << "\nand start the program again...." << endl;
                                                                    if (cinGetOnError==true) cin.get();
                                                                    exit(-1);
                                                        }
                                                    if (NodesGrid[x][y][i][AutLociNumber+2][0]<100)
                                                        {
                                                            GPFile<<0;
                                                                if(NodesGrid[x][y][i][AutLociNumber+2][0]<10)
                                                                    {
                                                                        GPFile<<0;
                                                                    }
                                                        }
                                                    GPFile<<NodesGrid[x][y][i][AutLociNumber+2][0]<<" ";
                                                    /******************************///Adaptation Locale
                                                    for (int p(0);p<AutLociNumber;p++)
                                                        {
                                                            for (int k(0);k<2;k++)
                                                                {
                                                                    GPFile<<0<<0<<NodesGrid[x][y][i][AutLociNumber+3+p][k]+1;
                                                                }
                                                            GPFile<<" ";
                                                        }
                                                    GPFile<<endl;
                                                }//end for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                        }
                                }
                        }
                }// end for(int x(0);x<signed(DimX);x++)
            GPFile.close();
        }//end of GenepopIntrog writing

    if (WriteGenepopOrigin==true)//On fait ca pour n'avoir que les habitats d'origine des noeuds, a la place des alleles
        {
            stringstream stst;
            stst << RUN; //transforme entier en strstr
            string fileName("GenepopOrigin_");
            fileName+=Contact;
            fileName+=stst.str();
            fileName+=".txt";
            ofstream GPFile(fileName.c_str());
            GPFile<<"This file has been generated by the AllForward software"<<endl;
            for(int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AutLocus"<<p<<endl;
                }
            GPFile<<"ZLocus"<<endl;
            GPFile<<"WLocus"<<endl;
            GPFile<<"MtLocus"<<endl;
            for (int p(0);p<AutLociNumber;p++)
                {
                    GPFile<<"AdaptLocus"<<p<<endl;
                }
            for(int x(0);x<signed(DimX);x++)
                {
                    for(int y(0);y<signed(DimY);y++)
                        {
                            if(NodesGrid[x][y].size()!=0)
                                {
                                    GPFile<<"pop"<<endl;
                                    for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                        {
                                            GPFile<<x<<" "<<y<<" , ";
                                            for(int p(0);p<AutLociNumber;p++)
                                                {
                                                    for(int a(0);a<2;a++)
                                                        {
                                                            GPFile<<0<<Alleles[p][NodesGrid[x][y][i][p][a]].habitat+1;//
                                                        }
                                                    GPFile<<" ";
                                                }
                                            GPFile<<0<<0<<Alleles[AutLociNumber][NodesGrid[x][y][i][AutLociNumber][0]].habitat+1;//premier Z
                                            if (NodesGrid[x][y][i][AutLociNumber][1]!=0)
                                                {
                                                        GPFile<<0<<0<<Alleles[AutLociNumber][NodesGrid[x][y][i][AutLociNumber][1]].habitat+1<<" ";//second Z
                                                        GPFile<<0<<0<<0<<" ";//W fantome
                                                }
                                            else
                                                {
                                                    GPFile<<0<<0<<0<<" ";//Z fantome
                                                    GPFile<<0<<0<<Alleles[AutLociNumber+1][NodesGrid[x][y][i][AutLociNumber+1][0]].habitat+1<<" ";//W
                                                }
                                            GPFile<<0<<0<<Alleles[AutLociNumber+2][NodesGrid[x][y][i][AutLociNumber+2][0]].habitat+1<<" ";//Mt

                                            for (int p(0);p<AutLociNumber;p++)
                                                {
                                                    for (int k(0);k<2;k++)
                                                        {
                                                            GPFile<<0<<NodesGrid[x][y][i][AutLociNumber+3][k]+1;// local adaptation gene
                                                        }
                                                    GPFile<<" ";
                                                }
                                            GPFile<<endl;
                                        }//end for(unsigned int i(0);i<NodesGrid[x][y].size();i++)//boucle sur les individus
                                }
                        }
                }// end for(int x(0);x<signed(DimX);x++)
        }//end if (WriteGenepopOrigin==true)
    return 0;
}//  end FGenepopFile(vector<map<int,CAlleles> >& Alleles, vector<vector<vector<vector<vector<int> > > > >& NodesGrid, unsigned int const& RUN)


/******************************************************************/
int FIntrogressionStats(vector<map<int,CAlleles> > const& Alleles, vector<vector<vector<vector<vector<int> > > > > const& NodesGrid, unsigned int const& RUN, unsigned long const& years)
{
    if (WriteIntrogProfile==true)
        {
            vector<vector<double> > Type0Prop;//contiendra proportion d'alleles type 0 par locus et par x
            Type0Prop.resize(DimX);
            vector<vector<int> > TotalGenes;//nombre de genes par locus et par x
            TotalGenes.resize(DimX);

            for (unsigned int x(0); x<DimX;x++)//Pour le profil d'introgression
                {
                    Type0Prop[x].resize(2*AutLociNumber+3);//autosomes plus Z, W, la mt et les adapt locale
                    TotalGenes[x].resize(2*AutLociNumber+3);//autosomes plus Z, W, la mt et les adapt locale
                    for (unsigned int y(0);y<DimY;y++)
                        {
                            for( unsigned int i(0);i<NodesGrid[x][y].size();i++)
                                {

                                    for(int p(0);p<AutLociNumber;p++)
                                        {
                                            for(unsigned int a(0);a<2;a++)
                                                {
                                                    TotalGenes[x][p]++;
                                                    if (Alleles[p].find(NodesGrid[x][y][i][p][a])->second.habitat==0)
                                                        {
                                                            Type0Prop[x][p]++;
                                                        }
                                                }
                                        }
                                    //Z Chromosome
                                    TotalGenes[x][AutLociNumber]++;
                                    if(Alleles[AutLociNumber].find(NodesGrid[x][y][i][AutLociNumber][0])->second.habitat==0)
                                        {
                                            Type0Prop[x][AutLociNumber]++;
                                        }
                                    if(NodesGrid[x][y][i][AutLociNumber][1]!=0)//male
                                        {
                                            TotalGenes[x][AutLociNumber]++;
                                            if(Alleles[AutLociNumber].find(NodesGrid[x][y][i][AutLociNumber][0])->second.habitat==0)
                                                {
                                                    Type0Prop[x][AutLociNumber]++;
                                                }
                                        }
                                    else//female
                                        {//W Chromosome
                                            TotalGenes[x][AutLociNumber+1]++;
                                            if(Alleles[AutLociNumber+1].find(NodesGrid[x][y][i][AutLociNumber+1][0])->second.habitat==0)
                                                {
                                                    Type0Prop[x][AutLociNumber+1]++;
                                                }
                                        }
                                    //Mitochondria
                                    TotalGenes[x][AutLociNumber+2]++;
                                    if (Alleles[AutLociNumber+2].find(NodesGrid[x][y][i][AutLociNumber+2][0])->second.habitat==0)
                                        {
                                            Type0Prop[x][AutLociNumber+2]++;
                                        }
                                    //Local Adaptation
                                    for(int p(0);p<AutLociNumber;p++)
                                        {
                                            for(int a(0);a<2;a++)
                                                {
                                                    TotalGenes[x][AutLociNumber+3+p]++;
                                                    if(NodesGrid[x][y][i][AutLociNumber+3+p][a]==0)
                                                        {
                                                            Type0Prop[x][AutLociNumber+3+p]++;
                                                        }
                                                }
                                        }
                                }
                        }
                }// end for (unsigned int x(0); x<DimX;x++)
            vector<double> IntroX;
            IntroX.resize(DimX);
            vector<double> IntroLocus;
            IntroLocus.resize(2*AutLociNumber+3);
            ofstream IntrogProfile("IntrogProfile.txt", ios::app);
            IntrogProfile<<setprecision(3);
            if (RUN==1 && (years==0 || WritePeriod==0))//header une seule fois
                {
                    IntrogProfile<<"#AllForward output file"<<endl;
                    IntrogProfile<<"#Simulation Parameters :"<<endl;
                    IntrogProfile<<"#DemeSize="<<DemeSize*2<<"\n#DimX="<<DimX<<"\n#DimY="<<DimY<<"\n#Xlimit="<<Xlimit<<"\n#Generation Number="<<GenerationNumber<<"\n#Allopatry last="<<AllopatryLast<<endl;// DemeSize must be multiplied by 2 because for the user it is the number of ind and for me the number of couple
                    IntrogProfile<<"#HabitatSlideBegin"<<HabitatSlideBegin<<"\n#HabitatSlideEnd"<<HabitatSlideEnd<<"\n#HabitatSlideDepth"<<HabitatSlideDepth<<endl;
                    IntrogProfile<<"#DemeSamplingRatio="<<DemeSamplingRatio<<"\n#IndMeanSampled="<<IndMeanSample<<endl;
                    IntrogProfile<<"#dispmax="<<DispMax<<"\n#mFemale="<<mFemale[0]<< ", "<<mFemale[1] <<"\n#geomFemale="<<geomFemale[0]<<", "<<geomFemale[1]<<"\n#mMale="<<mMale[0]<<", "<<mMale[1]<<"\n#geomMale="<<geomMale[0]<<", "<<geomMale[1]<<"\n#EdgeEffects="<<EdgeEffects<<endl;
                    IntrogProfile<<"#FitnessNormal\t"<<"FitnessMaladaptation\t"<<"FitnessHybridFemale\t"<<"FitnessHybridMale\t"<<endl;
                    for (int i(0);i<AutLociNumber;i++)
                        {
                            IntrogProfile<<"#\t"<<FitnessNormal[i]<<"\t"<<FitnessMaladaptation[i]<<"\t"<<FitnessHybridFemale[i]<<"\t"<<FitnessHybridMale[i]<<endl;
                        }
                    IntrogProfile<<"#Mitochondrial_Maladapation_Fitness="<<FitnessMt<<endl;
                    IntrogProfile<<"#HybridNB="<<HybridNb<<"  (-1=infinity)"<<endl;
                    IntrogProfile<<"#ChoosyFemale="<<ChoosyFemale<<endl;
                    IntrogProfile<<"#AcceptRates=";
                    for (int i(0);i<9;i++)
                        {
                            IntrogProfile<<AcceptRates[i]<<" ; ";
                        }
                    IntrogProfile<<endl;
                    IntrogProfile<<"#MuRate="<<MuRate<<endl;
                    IntrogProfile<<"#InterRecombiRate="<<InterRecombiRate<<endl<<"#IntraRecombiRate="<<IntraRecombiRate<<endl;
                    IntrogProfile<<endl<<endl;
                    IntrogProfile<<"Run\tYear\tx\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            IntrogProfile<<"Locus"<<p<<"\t";
                        }
                    IntrogProfile<<"LocusZ\t";
                    IntrogProfile<<"LocusW\t";
                    IntrogProfile<<"LocusMt\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            IntrogProfile<<"LocusAdapt"<<p<<"\t";
                        }
                    IntrogProfile<<endl;
                }
            for (unsigned int x(0);x<DimX;x++)
                {
                    if(TotalGenes[x][0]!=0)//sinon, c'est que le deme n'est pas echantillonne
                        {
                            IntrogProfile<<RUN<<"\t"<<years<<"\t"<<x<<"\t";
                            for (int p(0);p<2*AutLociNumber+3;p++)//dont Z, W et la Mt
                                {
                                    if(TotalGenes[x][p]!=0)
                                        {
                                            Type0Prop[x][p]/=TotalGenes[x][p];
                                            IntroX[x]+=Type0Prop[x][p]/(AutLociNumber+1);
                                            IntroLocus[p]+=Type0Prop[x][p]/DimX;
                                            IntrogProfile<<fixed<<Type0Prop[x][p]<<"\t";
                                        }
                                    else IntrogProfile<<fixed<<"NA\t";
                                }
                            IntrogProfile<<endl;
                        }
                }
            IntrogProfile.close();
        }// end if (WriteIntrogProfile==true)

    //Autres stats intorgression
    if(WriteIntrogStats==true)
        {
            vector<vector<double> > SpIntroLocus;//Pourcentage de copies introgresses par locus, chez les deux especes, dans les zones non hybrides
            SpIntroLocus.resize(2);
            vector<vector<double> > SpTotGenes;//total de copies par locus, chez les deux especes, dans les zones non hybrides
            SpTotGenes.resize(2);
            for (int s(0);s<2;s++)
                {
                    SpIntroLocus[s].resize(AutLociNumber+3);
                    SpTotGenes[s].resize(AutLociNumber+3);
                        for(int p(0);p<AutLociNumber+3;p++)
                            {
                                SpIntroLocus[s][p]=0.;
                                SpTotGenes[s][p]=0.;
                            }
                }
            for (unsigned int x(0);x<DimX;x++)
            for (unsigned int y(0);y<DimY;y++)
                {
                    for( unsigned int i(0);i<NodesGrid[x][y].size();i++)
                        {
                            if (NodesGrid[x][y][i][AutLociNumber+3][0]==NodesGrid[x][y][i][AutLociNumber+3][1])//la c'est pour le % de loci avec de l'introgression, par espece, et pour une zone a priori libre d'hybrides. On elimine quand meme les heterozygotes au locus d'adaptation locale
                                {
                                    for(int p(0);p<AutLociNumber;p++)
                                        {
                                            for(unsigned int a(0);a<NodesGrid[x][y][i][p].size();a++)
                                                {
                                                    if ((x<LowHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==0))
                                                        {
                                                            SpTotGenes[0][p]++;
                                                            if (Alleles[p].find(NodesGrid[x][y][i][p][a])->second.habitat==1)//ie pas dans son habitat
                                                                {
                                                                    SpIntroLocus[0][p]++;
                                                                }
                                                        }
                                                    if ((x>HighHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==1))
                                                        {
                                                            SpTotGenes[1][p]++;
                                                            if (Alleles[p].find(NodesGrid[x][y][i][p][a])->second.habitat==0)//ie pas dans son habitat
                                                                {
                                                                    SpIntroLocus[1][p]++;
                                                                }
                                                        }
                                                }
                                        }
                                    int A(2);//Sex chromosomes
                                    if (NodesGrid[x][y][i][AutLociNumber][1]==0)//female
                                        {
                                            if ((x<LowHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==0))//W chromosome
                                                {
                                                    SpTotGenes[0][AutLociNumber+1]++;
                                                    if (Alleles[AutLociNumber+1].find(NodesGrid[x][y][i][AutLociNumber+1][0])->second.habitat==1)//ie pas dans son habitat
                                                        {
                                                            SpIntroLocus[0][AutLociNumber+1]++;
                                                        }
                                                }
                                            if ((x>HighHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==1))
                                                {
                                                    SpTotGenes[1][AutLociNumber+1]++;
                                                    if (Alleles[AutLociNumber+1].find(NodesGrid[x][y][i][AutLociNumber+1][0])->second.habitat==0)//ie pas dans son habitat
                                                        {
                                                            SpIntroLocus[1][AutLociNumber+1]++;
                                                        }
                                                }
                                            A=1;
                                        }
                                    for (int a(0);a<A;a++)//Z chromosome
                                        {
                                            if ((x<LowHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==0))
                                                {
                                                    SpTotGenes[0][AutLociNumber]++;
                                                    if (Alleles[AutLociNumber].find(NodesGrid[x][y][i][AutLociNumber][a])->second.habitat==1)//ie pas dans son habitat
                                                        {
                                                            SpIntroLocus[0][AutLociNumber]++;
                                                        }
                                                }
                                            if ((x>HighHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==1))
                                                {
                                                    SpTotGenes[1][AutLociNumber]++;
                                                    if (Alleles[AutLociNumber].find(NodesGrid[x][y][i][AutLociNumber][a])->second.habitat==0)//ie pas dans son habitat
                                                        {
                                                            SpIntroLocus[1][AutLociNumber]++;
                                                        }
                                                }
                                        }
                                    //Mitochondria
                                    if ((x<LowHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==0))
                                        {
                                            SpTotGenes[0][AutLociNumber+2]++;
                                            if (Alleles[AutLociNumber+2].find(NodesGrid[x][y][i][AutLociNumber+2][0])->second.habitat==1)//ie pas dans son habitat
                                                {
                                                    SpIntroLocus[0][AutLociNumber+2]++;
                                                }
                                        }
                                    if ((x>HighHybridBound)&&(NodesGrid[x][y][i][AutLociNumber+3][0]==1))
                                        {
                                            SpTotGenes[1][AutLociNumber+2]++;
                                            if (Alleles[AutLociNumber+2].find(NodesGrid[x][y][i][AutLociNumber+2][0])->second.habitat==0)//ie pas dans son habitat
                                                {
                                                    SpIntroLocus[1][AutLociNumber+2]++;
                                                }
                                        }
                                }//end Si pas hybrides
                        }
                }//end for (unsigned int x(0); x<DimX;x++)
            vector<unsigned int> HybridZone;
            HybridZone.resize(2);
            for (int h(0);h<2;h++) HybridZone[h]=Xlimit;
            for (unsigned int x(0);x<DimX;x++)
            for (unsigned int y(0);y<DimY;y++)
                {
                    for (unsigned int i(0);i<NodesGrid[x][y].size();i++)
                        {
                            if (NodesGrid[x][y][i][AutLociNumber+3][0]!=NodesGrid[x][y][i][AutLociNumber+3][1])//adapt locale
                                {
                                    if(x<HybridZone[0]) HybridZone[0]=x;
                                    if(x>HybridZone[1]) HybridZone[1]=x;
                                }
                        }
                }
            vector<double> ProportionLociIntrog;
            ProportionLociIntrog.resize(2);
            for (int s(0);s<2;s++)
                {
                    for (int p(0);p<AutLociNumber+2;p++)//Tout les loci sauf la mitochondrie
                        {
                            if (SpIntroLocus[s][p]>0.)
                                {
                                    ProportionLociIntrog[s]+=(1./(AutLociNumber+2));
                                }
                            if (SpTotGenes[s][p]>0)
                                {
                                    SpIntroLocus[s][p]/=SpTotGenes[s][p];
                                }
                        }
                    if (SpTotGenes[s][AutLociNumber+2]>0)//mitochondria
                        {
                            SpIntroLocus[s][AutLociNumber+2]/=SpTotGenes[s][AutLociNumber+2];
                        }
                }

            ofstream IntrogStats("IntrogStats.txt", ios::app);
            IntrogStats<<setprecision(3);
            if (RUN==1 && (years==0 || WritePeriod==0))
                {
                    IntrogStats<<"#AllForward output file"<<endl;
                    IntrogStats<<"#Simulation Parameters :"<<endl;
                    IntrogStats<<"#DemeSize="<<DemeSize*2<<"\n#DimX="<<DimX<<"\n#DimY="<<DimY<<"\n#Xlimit="<<Xlimit<<"\n#Generation Number="<<GenerationNumber<<"\n#Allopatry last="<<AllopatryLast<<endl;// DemeSize must be multiplied by 2 because for the user it is the number of ind and for me the number of couple
                    IntrogStats<<"#HabitatSlideBegin"<<HabitatSlideBegin<<"\n#HabitatSlideEnd"<<HabitatSlideEnd<<"\n#HabitatSlideDepth"<<HabitatSlideDepth<<endl;
                    IntrogStats<<"#DemeSamplingRatio="<<DemeSamplingRatio<<"\n#IndMeanSampled="<<IndMeanSample<<endl;
                    IntrogStats<<"#dispmax="<<DispMax<<"\n#mFemale="<<mFemale[0]<< ", "<<mFemale[1] <<"\n#geomFemale="<<geomFemale[0]<<", "<<geomFemale[1]<<"\n#mMale="<<mMale[0]<<", "<<mMale[1]<<"\n#geomMale="<<geomMale[0]<<", "<<geomMale[1]<<"\n#EdgeEffects="<<EdgeEffects<<endl;
                    IntrogStats<<"#FitnessNormal\t"<<"FitnessMaladaptation\t"<<"FitnessHybridFemale\t"<<"FitnessHybridMale\t"<<endl;
                    for (int i(0);i<AutLociNumber;i++)
                        {
                            IntrogStats<<"#\t"<<FitnessNormal[i]<<"\t"<<FitnessMaladaptation[i]<<"\t"<<FitnessHybridFemale[i]<<"\t"<<FitnessHybridMale[i]<<endl;
                        }
                    IntrogStats<<"#Mitochondrial_Maladapation_Fitness="<<FitnessMt<<endl;
                    IntrogStats<<"#HybridNB="<<HybridNb<<"  (-1=infinity)"<<endl;
                    IntrogStats<<"#ChoosyFemale="<<ChoosyFemale<<endl;
                    IntrogStats<<"#AcceptRates=";
                    for (int i(0);i<9;i++)
                        {
                            IntrogStats<<AcceptRates[i]<<" ; ";
                        }
                    IntrogStats<<endl;
                    IntrogStats<<"#MuRate="<<MuRate<<endl;
                    IntrogStats<<"#InterRecombiRate="<<InterRecombiRate<<endl<<"#IntraRecombiRate="<<IntraRecombiRate<<endl;
                    IntrogStats<<"#Stats calculated within "<<0<<" ; "<<LowHybridBound<<" and "<<HighHybridBound<<" ; "<<DimX-1<<endl<<endl;
                    IntrogStats<<"#Mean number of introgressed genes by locus and by run"<<endl;
                    IntrogStats<<"Run\tYear\tSp\t";
                    for (int p(0);p<AutLociNumber;p++)
                        {
                            IntrogStats<<"Aut"<<p<<"\t";
                        }
                    IntrogStats<<"Zchr\tWchr\t";
                    IntrogStats<<"Nucl_Intro\t";
                    IntrogStats<<"Mt\t";
                    IntrogStats<<"HybridZone"<<endl;
                }

            for (int s(0);s<2;s++)
                {
                    IntrogStats<<RUN<<"\t"<<years<<"\t"<<s<<"\t";
                    for (int p(0);p<AutLociNumber+2;p++)//tout les loci nucleaires
                        {
                            if (SpTotGenes[s][p]>0)
                                {
                                    IntrogStats<<fixed<<SpIntroLocus[s][p]<<"\t";
                                }
                            else IntrogStats<<"NA\t";
                        }
                    IntrogStats<<fixed<<ProportionLociIntrog[s]<<"\t";
                    IntrogStats<<fixed<<SpIntroLocus[s][AutLociNumber+2]<<"\t";//mitochondria
                    IntrogStats<<HybridZone[s]<<endl;//limites de la zone hybride
                }
            IntrogStats.close();
        }// end if(WriteIntrogStats==true)
    return 0;
}// end int FIntrogressionStats(map<int,Cnode>& MtNodes, vector<map<int,Cnode> >& AutNodes, map<int,Cnode>& ZNodes, map<int,Cnode>& WNodes, vector<vector<vector<vector<vector<int> > > > >& NodesGrid, unsigned int const& RUN, vector<vector<vector<vector<bool> > > >& AdaptNodes)
