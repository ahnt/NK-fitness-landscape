#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <time.h>
#include <math.h>

using namespace std;

#define randDouble ((double)rand()/(double)RAND_MAX)
#define discreteMode
#define ED 5

vector<double> fitnessesForUIntGenotype;
map<vector<double>,double> fitnessForGenotypes;
map<vector<double>,int> genotypeCount;

double **NKTable;
unsigned int bitMask;
int N,K,A;
int update=0,updates;
double replacementRate;
double mutationRate;
unsigned int nrOfGenotypes;

class tMutation{
public:
    int from,to,where;
    void setup(int f, int w, int t){from=f; to=t; where=w;}
};

class tAgent{
public:
#ifdef discreteMode
    vector<int> genome;
#else
    vector<double> genome;
#endif
    double fitness;
    int born;
    int nrPointingAtMe;
    tAgent *ancestor;
    void inherit(tAgent *from);
    tAgent();
    ~tAgent();
    void setupDiscrete(unsigned int newGenotype);
    void setupProbabilistic(void);
    void saveLOD(FILE* lod);
    void show(void);
    void computeOwnFitness(void);
};

vector<tAgent*> lodBuffer;


double genotypeToFitness(unsigned int genotype);
void showNKTable(void);
void doEpistasis(FILE *f);

int main (int argc, char * const argv[]) {
    int i,j;
    double maxFitness;
    double lowest=1.0;
    unsigned int IDofLowest=0;
    vector<tAgent*> population;
	N=atoi(argv[1]);				//N
	K=atoi(argv[2]);				//K
	A=atoi(argv[3]);				//popsize;
    nrOfGenotypes=(unsigned int)1<<(unsigned int)N;
	updates=atoi(argv[4]);			//updates
	replacementRate=atof(argv[5]);	//replacement rate
	mutationRate=atof(argv[6]);		//mutation rate
	bitMask=(1<<K)-1;
	cout<<N<<" "<<K<<" "<<bitMask<<" "<<nrOfGenotypes<<endl;
	cout<<A<<" "<<RAND_MAX<<endl;
	srand((unsigned int)time(NULL));
    if(strcmp(argv[8],"RAND")==0){
        //make the NKTable
        cout<<"make the table"<<endl;
        NKTable=NULL;
        NKTable=(double**)malloc(sizeof(double*)*N);
        for(i=0;i<N;i++){
            NKTable[i]=(double*)malloc(sizeof(double)*(1<<K));
            for(j=0;j<(1<<K);j++)
                //don't produce zeros as entries
                do{
                    NKTable[i][j]=randDouble;
                }while(NKTable[i][j]==0.0);
        }
        //showNKTable();
        //make all fitnesses for genotype table
        cout<<"make all fitnesses"<<endl;
        fitnessesForUIntGenotype.resize(nrOfGenotypes);
        cout<<fitnessesForUIntGenotype.size()<<endl;
        for(unsigned int u=0;u<nrOfGenotypes;u++){
            fitnessesForUIntGenotype[u]=genotypeToFitness(u);
            // cout<<genotypeToFitness(u)<<endl;
            if(fitnessesForUIntGenotype[u]<lowest){
                lowest=fitnessesForUIntGenotype[u];
                IDofLowest=u;
            }
        }
        FILE *f=fopen("fls.txt","w+t");
        for(i=0;i<nrOfGenotypes;i++)
            fprintf(f, "%f\n",fitnessesForUIntGenotype[i]);
        fclose(f);
        f=fopen("NKT_fls.txt","w+t");
        for(i=0;i<N;i++){
            for(j=0;j<(1<<K);j++)
                fprintf(f,"%f   ",NKTable[i][j]);
            fprintf(f,"\n");
        }
        fclose(f);
    } else{
        fitnessesForUIntGenotype.clear();
        FILE *f=fopen(argv[8],"r+t");
        for(unsigned int u=0;u<nrOfGenotypes;u++){
            float d;
            fscanf(f,"%f\n",&d);
            fitnessesForUIntGenotype.push_back((double)d);
            if((double)d<lowest){
                lowest=(double)d;
                IDofLowest=u;
            }
        }
        fclose(f);
        
    }
    //initialize the population
    population.clear();
    for(i=0;i<A;i++){
        tAgent *dA=new tAgent;
#ifdef discreteMode
        dA->setupDiscrete(IDofLowest);
#else
        dA->setupProbabilistic();
#endif
        population.push_back(dA);
    }
    cout<<lowest<<" "<<population[0]->fitness<<" "<<IDofLowest<<endl;
    //main loop
	for(update=1;update<updates;update++){
        maxFitness=0.0;
        for(i=0;i<A;i++)
            if(population[i]->fitness>maxFitness)
                maxFitness=population[i]->fitness;
        for(i=0;i<A;i++)
            if(randDouble<replacementRate)
            {
                population[i]->nrPointingAtMe--;
                if(population[i]->nrPointingAtMe==0)
                    delete population[i];
                tAgent *dA=new tAgent;
                do{j=rand()%A;} while((i==j)||(population[j]->born==update)||(randDouble>population[j]->fitness/maxFitness));
                dA->inherit(population[j]);
                population[i]=dA;
            }
        cout<<update<<" "<<maxFitness<<" "<<genotypeCount.size()<< endl;
	}
    //save LOD
    FILE *lod=fopen(argv[7],"w+t");
    population[1]->saveLOD(NULL);
    //do Epistasis
    doEpistasis(lod);
    fclose(lod);
    return 0;
}


//regular functions
double genotypeToFitness(unsigned int genotype){
	int i,j,n;
	double fitness=0.0;
	n=0;
	for(j=0;j<K;j++)
		n=(n<<1)+((genotype>>j)&1);
	for(i=0;i<N;i++){
		fitness+=log(NKTable[i][n]);
		n=((n<<1)+((genotype>>((i+K)%N))&1))&bitMask;
	}
	fitness=exp(fitness/(double)N);
    /* //testfitness maker
    fitness=0.0;
    for(i=0;i<N;i++)
        fitness+=1.0/(double)N*(double)((genotype>>i)&1);
     */
	return fitness;
}

void showNKTable(void){
	int i,j;
	for(j=0;j<(1<<K);j++){
		cout<<j<<": ";
		for(i=0;i<N;i++)
			cout<<NKTable[i][j]<<" ";
		cout<<endl;
	}
    
}

void doEpistasis(FILE *f){
    int i,j,n,k;
    map<string,tMutation> mutationsMap;
    vector<string> mutations;
    set<string> setM;
    set<int> setI;
    i=0;
    while(i<(lodBuffer.size()-1)){
        k=0;
        for(j=0;j<N;j++){
            if(lodBuffer[i]->genome[j]!=lodBuffer[i+1]->genome[j])
                k++;
        }
        if(k==0){
            cout<<"--"<<endl;
            lodBuffer.erase(lodBuffer.begin()+i);
        }
        else
            i++;
    }
    for(i=0;i<lodBuffer.size()-ED;i++){
        mutations.clear();
        mutationsMap.clear();
        setM.clear();
        setI.clear();
        for(j=0;j<ED;j++){
            for(n=0;n<N;n++)
                if(lodBuffer[i+j]->genome[n]!=lodBuffer[i+j+1]->genome[n]){
                    string M;
                    tMutation tM;
                    char c[1000];
                    sprintf(c,"%i-%i-%i",lodBuffer[i+j]->genome[n],n,lodBuffer[i+j+1]->genome[n]);
                    tM.setup(lodBuffer[i+j]->genome[n],n,lodBuffer[i+j+1]->genome[n]);
                    M.assign(c);
                    mutations.push_back(M);
                    setM.insert(M);
                    setI.insert(n);
                    mutationsMap[M]=tM;
                }
        }
        fprintf(f,"%i   %f",lodBuffer[i]->born,lodBuffer[i]->fitness);
        for(j=0;j<N;j++)
            fprintf(f," %i",lodBuffer[i]->genome[j]);
        if((setI.size()!=ED)||(setM.size()!=mutations.size())||(mutationsMap.size()!=ED)){
            for(j=0;j<1<<ED;j++)
                fprintf(f,"    0.0");
            fprintf(f,"\n");
        }
        else{
                
            for(j=0;j<(1<<ED);j++){
                tAgent *A=new tAgent;
                A->genome=lodBuffer[i]->genome;
                for(k=0;k<ED;k++)
                    if((j&(1<<k))!=0){
                        //apply mutation
                        A->genome[mutationsMap[mutations[k]].where]=mutationsMap[mutations[k]].to;
                    }
                A->computeOwnFitness();
                fprintf(f,"    %f",A->fitness);
            }
            fprintf(f,"\n");
        }
    }
}


// tAgent definitions
void tAgent::inherit(tAgent *from){
    genome.resize(from->genome.size());
    from->nrPointingAtMe++;
    ancestor=from;
    bool doit=false;
    for(int i=0;i<from->genome.size();i++)
        if(randDouble<mutationRate){
#ifdef discreteMode
            genome[i]=rand()&1;
#else
            genome[i]=randDouble;
            doit=true;
#endif
        }
        else
            genome[i]=from->genome[i];
#ifdef discreteMode
    unsigned int G=0;
    for(int i=0;i<genome.size();i++)
        if(genome[i]==1)
            G|=(unsigned int)1<<(unsigned int)i;
    fitness=fitnessesForUIntGenotype[G];
#else
    if(doit){
        fitness=0.0;
        for(unsigned int u=0;u<nrOfGenotypes;u++){
            double p=1.0;
            for(int i=0;i<N;i++)
                if(((u>>i)&1)==1)
                    p*=genome[i];
                else
                    p*=(1.0-genome[i]);
            fitness+=fitnessesForUIntGenotype[u]*p;
        }
        fitnessForGenotypes[genome]=fitness;
    }
    else
        fitness=from->fitness;
#endif
}
tAgent::tAgent(){
    nrPointingAtMe=1;
    ancestor=0;
    born=update;
    fitness=0;
}

tAgent::~tAgent(){
    if(ancestor!=NULL){
        ancestor->nrPointingAtMe--;
        if(ancestor->nrPointingAtMe==0)
            delete ancestor;
    }
    /*
    if(--genotypeCount[genome]==0){
        genotypeCount.erase(genotypeCount.find(genome));
        fitnessForGenotypes.erase(fitnessForGenotypes.find(genome));
    }*/
}

void tAgent::setupDiscrete(unsigned int newGenotype){
    genome.resize(N);
    for(int i=0;i<N;i++){
        if(((newGenotype>>i)&1)==1)
            genome[i]=1;
        else
            genome[i]=0;
    }
    fitness=fitnessesForUIntGenotype[newGenotype];
}

void tAgent::setupProbabilistic(void){
    genome.resize(N);
    for(int i=0;i<N;i++)
        genome[i]=0.5;
}

void tAgent::saveLOD(FILE* lod){
    if(ancestor!=NULL){
        ancestor->saveLOD(lod);
        if(ancestor->nrPointingAtMe!=1){
            //fclose(lod);
            //exit(0);
        }
    }
    lodBuffer.push_back(this);
    if(lod!=NULL){
        fprintf(lod,"%i %f  ",born,fitness);
        for(int i=0;i<genome.size();i++){
#ifdef discreteMode
            fprintf(lod,"%i ",(int)genome[i]);
#else
            fprintf(lod,"%f ",genome[i]);
#endif
        }
        fprintf(lod,"\n");
    }
}

void tAgent::show(void){
    printf("%i   %f ",born,fitness);
    for(int i=0;i<N;i++){
#ifdef discreteMode
        printf("%i",(int)genome[i]);
#else
    printf("    %f",genome[i]);
#endif
    }
    printf("\n");
}
void tAgent::computeOwnFitness(void){
#ifdef discreteMode
    unsigned int G=0;
    for(int i=0;i<genome.size();i++)
        if(genome[i]==1)
            G|=(unsigned int)1<<(unsigned int)i;
    fitness=fitnessesForUIntGenotype[G];
#else
    if(true){
        fitness=0.0;
        for(unsigned int u=0;u<nrOfGenotypes;u++){
            double p=1.0;
            for(int i=0;i<N;i++)
                if(((u>>i)&1)==1)
                    p*=genome[i];
                else
                    p*=(1.0-genome[i]);
            fitness+=fitnessesForUIntGenotype[u]*p;
        }
        fitnessForGenotypes[genome]=fitness;
    }
    else
        fitness=from->fitness;
#endif
}

/*
 #ifdef discreteMode
 #else
 #endif
*/