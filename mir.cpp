#include "mir.h"
#include "stdio.h"
#include "string.h"
#include <algorithm>
#include <fstream>
#include <math.h>
#include <time.h>
#include <random>
#include <sstream>

using namespace lite;

// consts
int LOG_FREQ = 1000;
float Mir::maxSubstance = 1000;
float Org::maxEAT = 0.6;
float Org::maxEnergy = 1000;
int Org::maxAge = 3000;
const char Mir::alphabet[4] = {'A','T','G','C'};
int Mir::maxId = 0;

//////////////////////////////////////////////////////////////////////////
// Org

int Org::maxId = 0;

Org::Org(Mir* mir)
{
	this->mir = mir;
	id = to_string (maxId++) + "_" + to_string(mir->age);
	age = 0;
	energy = 0;
	parent = NULL;
	alive = tru/;
	name = "";
	SNPrate = 0;
}

Org::~Org()
{
}

Org* Org::divide()
{
	Org* newbie = new Org(mir);
	newbie->genome = genome;
	newbie->SNPrate = SNPrate;
	energy /= 2;
	newbie->energy = energy;
	// add mutagenes
	int sumLength = 0;
	for(int g = 0; g < newbie->genome.size(); g++) sumLength += newbie->genome[g].seq.size();

	float meanSNPs = sumLength*SNPrate;
	std::mt19937 generator;
	std::normal_distribution<double> normal(meanSNPs, meanSNPs);
	int SNPcount = (int)normal(generator);
	for(int i = 0; i < SNPcount; i++)
	{
		int g = rand()%newbie->genome.size();
		int pos = rand()%newbie->genome[g].seq.size();
		newbie->genome[g].seq[pos] = Mir::alphabet[rand()%Mir::alphabetLength];
	}
	// 
	if(mir->bPhyloLog){
		newbie->parent = this;
		children.push_back(newbie);
	}
	return newbie;
}

float Org::meanFit()
{
	float fit = 0;
	for(int g = 0; g < genome.size(); g++) fit += genome[g].fit;
	return fit/genome.size();
}

bool Org::anyLivingChild(){
	for(auto c : children){
		if(c->alive || c->anyLivingChild()){
			return(true);
		}
	}
	return false;
}

void Org::maybeDelete(){
	if(alive) return;
	cout << "void Org::maybeDelete\n"<< name;
	for(auto c : children){
		c->parent = parent;
	}
	if(parent != NULL){
		auto me = find(parent->children.begin(), parent->children.end(), this); 
		if(me != parent->children.end()){
			parent->children.erase(me);
		}
	}
	if(!anyLivingChild()){
		
		if(parent != NULL) parent->maybeDelete();
		if(alive){ children.clear(); return; }		
		delete this;
	}
}

// :)
//////////////////////////////////////////////////////////////////////////
// Mir

Mir::Mir(int argc, char **argv)
{
	bPhyloLog = true;
	id = ++maxId;
	// files
	paramFile = (char *)"params.txt";
	constFile = (char *)"consts";
	popLogFile = (char *)"populationLog.txt";

	cout << "argc = " << argc << endl;
	// seed
	if(argc >= 2) randomSeed = atoi(argv[1]); else randomSeed = 1; //time(NULL);
	srand (randomSeed);
	// params, logs
	if(argc >= 3) paramFile = argv[2];
	if(argc >= 4) constFile = argv[3];
	if(argc >= 5) popLogFile = argv[4];
	cout << "params file: " << paramFile << " const file: " << constFile << " pop out file: " << popLogFile << endl;

	printf("new mir: %d\n", id);
	// consts
	age = 0;
	w = 100; h =100;
	NSubstances = 2;
	NSubstanceSources = 10;
	minSubstance_dE = -5;
	maxSubstance_dE = 5;
	orgStartCount = 1000;
	initialOrgGenes = 1;
	initialOrgEnergy = 70;
	energyDivide = initialOrgEnergy*10;
	initialGeneLength = 20;
	diffusion = 0.5;
	substanceDegrade = 0.99;
	expressionCost = 1;
	initialSNPrate = 0.1;
	SourceLifetime = 1000;
	Org::maxAge = 3000; // -1 - no max age
	MirLifetime = 100000000;
    sourceMaxIntensity = 100;
    bSaveGenomes = false;
    minAgeToDivide = 100;

    // not in config yet!
	minSourceRadius = 1;
	maxSourceRadius = 40;
	bSaveGenomes = true; // EXPERIMENTAL
    // log
    divideLogOn = false;

	// LOAD params
	loadConfig();
	// resize
	dE.resize(NSubstances,NSubstances);
	goldSeqs.resize(NSubstances,NSubstances);
	substances.resize(w,h, NSubstances);
	orgs.resize(w,h);
	sources.reserve(NSubstanceSources);
	// populate
	adam = NULL;
	printf("Mir constructed\n");
	// files

}

void Mir::init()
{
	populate_dE();
	populateSources();
	populateGoldSeqs();
	nullPole();
	populateOrgs();
	openLogFiles();
}

void Mir::loadConfig()
{
	FILE *sf;
	sf = fopen(paramFile, "r");
	if(!sf) { cout << "no config file - working on defaults.\n"; return; }
	char buf[1024]; float f;
	while(fscanf(sf, "%s %f", buf, &f) == 2)
	{
		if(strcmp(buf, "width") == 0) { w = (int)f; h = w; }
		//if(strcmp(buf, "height") == 0) h = (int)f;
		if(strcmp(buf, "substances") == 0) NSubstances = (int)f;
		if(strcmp(buf, "sources") == 0) NSubstanceSources = (int)f;
		if(strcmp(buf, "genes") == 0) initialOrgGenes = (int)f;
		if(strcmp(buf, "geneLength") == 0) initialGeneLength = (int)f;
		if(strcmp(buf, "orgs") == 0) orgStartCount = (int)f;
		if(strcmp(buf, "min_dE") == 0) minSubstance_dE = f;
		if(strcmp(buf, "max_dE") == 0) maxSubstance_dE = f;
		if(strcmp(buf, "startEnergy") == 0) initialOrgEnergy = f;
		if(strcmp(buf, "energyToDivide") == 0) energyDivide = f;
		if(strcmp(buf, "minAgeToDivide") == 0) minAgeToDivide = (int)f;
		if(strcmp(buf, "diffusion") == 0) diffusion = f;
		if(strcmp(buf, "substanceDegrade") == 0) substanceDegrade = f;

		if(strcmp(buf, "expressionCost") == 0) expressionCost = f;
		if(strcmp(buf, "SNPrate") == 0) initialSNPrate = f;
		if(strcmp(buf, "sourceRadius") == 0) { minSourceRadius = (int)f; maxSourceRadius = minSourceRadius; }
        if(strcmp(buf, "sourceMaxIntensity") == 0) sourceMaxIntensity = f;


	 if(strcmp(buf, "maxAge") == 0) Org::maxAge = (int)f;
	 if(strcmp(buf, "SourceLifetime") == 0) SourceLifetime = (int)f;
	}
	fclose(sf);
	sf = fopen(constFile, "r");
	if(!sf) { cout << "no consts file - working on defaults.\n"; return; }
	while(fscanf(sf, "%s %f", buf, &f) == 2)
	{
		if(strcmp(buf, "MirLifetime") == 0) MirLifetime = (int)f;
		if(strcmp(buf, "substances") == 0) NSubstances = (int)f;
		if(strcmp(buf, "LogFreq") == 0) LOG_FREQ = (int)f;
		if(strcmp(buf, "genes") == 0) initialOrgGenes = (int)f;
		if(strcmp(buf, "saveGenomes") == 0) bSaveGenomes = (bool)f;
}
	fclose(sf);
}


Mir::~Mir()
{
}

void Mir::degrade()
{
	for(int i = 0; i < orgsVector.size(); i++)
		delete orgsVector[i];
	orgsVector.clear();
	sources.clear();
	nullPole();
	closeLogFiles();
	age = 0;
	delete adam;
}

void Mir::tic()
{
	sourcesEmmit();
	sourceReincarnate();
	//diffuse();
	orgEat();
	orgDie();
	orgDivide();
	if(orgsVector.size() == 0) { id++; populateOrgs(); }
	logPopulation();
	age++;
}

void Mir::main()
{
	int echoTimer = 0;
	while(age <= MirLifetime)
	{
	 tic();
	 if(++echoTimer  > 1000) {
	     printf("mean fit: %f\torgs: %d\tmirAge:%d\n", meanEnzymeFit(), (int)orgsVector.size(), age);
	     echoTimer = 0;
	 }
	}
	if(bPhyloLog){
		giveNames(adam);
		saveGenomes();
	}
}

Org* Mir::org(int x, int y)
{
	for(int o = 0; o < orgsVector.size(); o++)
	{
		if(orgsVector[o] == orgs(x,y))
		{
			return orgsVector[o];
		}
	}
	return NULL;
}

//////////////////////////////////////////////////////////////////////////////

void Mir::reportDE()
{
	printf("*** dE matrix ***\n");
	for(int i = 0; i < NSubstances; i++){
		for(int j = 0; j < NSubstances; j++){
			printf("%f\t", dE(i,j));
		}
		printf("\n");
	}
}

//////////////////////////////////////////////////////////////////////////////
// init

void Mir::nullPole()
{
	// null orgs and substances
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			orgs(i,j) = NULL;
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			for(int s = 0; s < NSubstances; s++)
	substances(i,j,s) = 0;
}

void Mir::populateOrgs()
{
	id = maxId++; // experimental
	//////////////////////////////////////////////////////////////////////////
	orgs.resize(w,h);
	orgsVector.clear();
	orgsVector.reserve(orgStartCount);
	delete(adam);
	adam = new Org(this);
	adam->name = "adam";
	adam->alive = false;	
	for(int i = 0; i < orgStartCount; i++)
	{
		int x = rand()%w;
		int y = rand()%h;
		if(orgs(x,y) != NULL) continue;
		Org* newbie = new Org(this);
		newbie->x = x;
		newbie->y = y;
		newbie->energy = rand()%(int)initialOrgEnergy;
		// make genome
		newbie->genome.resize(initialOrgGenes);
		for(int g = 0; g < initialOrgGenes; g++)
		{
			newbie->genome[g].seq = randomSeq(initialGeneLength);
			//cout << newbie->genome[g].seq << endl;
			//newbie->genome[g].seq = goldSeqs(0,1);
			determineEnzyme(newbie->genome[g]);
		}
		newbie->SNPrate = initialSNPrate;
		newbie->parent = adam;
		
		adam->children.push_back(newbie);
		orgs(x,y) = newbie;
		orgsVector.push_back(newbie);
	}
}

void Mir::createRandomSource(SubstanceSource& source)
{
	int goods = goodSubstances.size();
	if(goods <= 0) return;
	source.x = rand()%w;
	source.y = rand()%h;
	int deltaRadius = maxSourceRadius - minSourceRadius;
	if(deltaRadius > 0) source.radius = rand()%deltaRadius + minSourceRadius;
        else source.radius = maxSourceRadius;
	source.substanceId = goodSubstances[ rand()%goods ];
	source.intensity = rand()%sourceMaxIntensity;
	source.age = 0;
}

void Mir::populateSources()
{
	SubstanceSource source;
	if(goodSubstances.size() == 0) return;
	for(int s = 0; s < NSubstanceSources; s++)
	{
		createRandomSource(source);
		sources.push_back(source);
	}
}

void Mir::populate_dE()
{
	for(int i = 0; i < NSubstances; i++)
		for(int j = 0; j <= i; j++){
			if(i == j)
	 dE(i,j) = 0;
			else
	dE(i,j) = (rand()%(100*(int)(maxSubstance_dE - minSubstance_dE)))/100.0 + minSubstance_dE;
	dE(j,i) = -dE(i,j);
		}
	reportDE();
	goodSubstances.clear();
	for(int i =0 ; i < NSubstances; i++)
		if(!badSubstance(i)) goodSubstances.push_back(i);
	cout << "good substances: " << goodSubstances.size() << endl;
}

string Mir::randomSeq(int length)
{
	string tmp;
	tmp.resize(length);
	for(int i = 0; i < length; i++)
	{
		tmp[i] = alphabet[rand()%alphabetLength];
	}
	return tmp;
}

void Mir::populateGoldSeqs()
{
	for(int i = 0; i < NSubstances; i++)
		for(int j = 0; j < NSubstances; j++){
			if(i == j) continue;
            goldSeqs(i,j) = randomSeq(initialGeneLength);
		}
}

float stringDist(string a, string b)
{
	int sizeDiff = a.size() - b.size();
	int minSize = min(a.size(), b.size());
	int matched = 0;
	int count = 0;
	for(int i = 0; i < minSize; i++)
	{
		if(i%3 == 0) continue; // synonymus position
		if(a[i] == b[i]) matched++;
		count++;
	}
	return (1 - (matched/(float)count) + 0.2*sizeDiff);
}

void Mir::determineEnzyme(Gene& gene)
{
	float minDist = 100000000;
	for(int i = 0; i < NSubstances; i++)
		for(int j = 0; j < NSubstances; j++)
		{
			if(i == j) continue;
			float dist = stringDist(gene.seq, goldSeqs(i,j));
			if(dist == 0) { gene.in = i; gene.out = j; gene.fit = 1; return; }
			if(dist < minDist){ minDist = dist; gene.in = i; gene.out = j; }
		}
	gene.fit = max(1-minDist, 0.0f);
}

//////////////////////////////////////////////////////////////////////////

void Mir::sourcesEmmit()
{
	for(int s = 0; s < sources.size(); s++)
	{
		SubstanceSource source = sources[s];
        // 2bd: add check 4 maximun subst value
        //if(substances(source.x, source.y, source.substanceId) >= maxSubstance)
		//	substances(source.x, source.y, source.substanceId) = maxSubstance;
        int radius = source.radius;
        if(radius == 0) substances(source.x, source.y, source.substanceId) += source.intensity;
        else
        {
            for(int dx = -radius; dx < radius; dx++)
                for(int dy = -radius; dy < radius; dy++)
                {
                   int X, Y;
                   X = source.x + dx; Y = source.y + dy;
                   putToWorld(X, Y);
                   substances(X, Y, source.substanceId) += source.intensity;
                }
        }
		sources[s].age++;
	}
}

void Mir::sourceReincarnate()
{
	if(SourceLifetime < 0) return;
	for(int s = 0; s < sources.size(); s++)
	{
		if(rand()%SourceLifetime == 1)
		{
			createRandomSource(sources[s]);
		}
	}
}

void Mir::diffuse()
{ // substanceDegradeDegrade
	lite::array<float[1]> neighbourSum(NSubstances);
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
		{
			int count = 0, I, J;
			neighbourSum = 0;
			for(int d_x = -1; d_x <= 1; d_x++)
                for(int d_y = -1; d_y <= 1; d_y++)
                {
                    if(d_x == 0 && d_y == 0) continue;
                    I = i + d_x; J = j + d_y;
                    putToWorld(I, J);
                    neighbourSum = neighbourSum + substances[row(I)][row(J)];
                    count++;
                }
			neighbourSum /= count;
			substances[row(i)][row(j)] = (1 - diffusion)*substances[row(i)][row(j)] + neighbourSum*diffusion;
			substances[row(i)][row(j)] *= substanceDegrade;
		}
}

void Mir::orgEat()
{
	for(int o = 0; o < orgsVector.size(); o++)
	{
		Org* org = orgsVector[o];
		if(!org->alive) cout << "deadman!!\n";
		for(int g = 0; g < org->genome.size(); g++)
		{
			Gene* G = &(org->genome[g]);
			float eated = min(substances(org->x, org->y, G->in), Org::maxEAT);
			org->energy += G->fit*eated*dE(G->in, G->out);
			substances(org->x, org->y, G->in) -= eated;
			substances(org->x, org->y, G->out) += eated*G->fit;

			// gene expression cost
			org->energy -= expressionCost;
		}
		org->age++;
		if(org->energy >= Org::maxEnergy) org->energy = Org::maxEnergy;
	}
}

//bool energyLow(Org *org)
//{
//	if(org->energy > 0) return false;
//	return true;
//}
/*
void Mir::orgDie()
{
	//org->age++;

	vector<Org*>::iterator dead = remove_if(orgsVector.begin(),orgsVector.end(), energyLow);
	vector<Org*>::iterator d = dead;
	while(d != orgsVector.end())
	{
			orgs((*d)->x, (*d)->y) = NULL;
			delete *d;
			d++;
	}
	orgsVector.erase(dead, orgsVector.end());

}
*/

void Mir::orgDie()
{
	int o = orgsVector.size() - 1;
	while(o >= 0 && o != orgsVector.size())
	{
		if(orgsVector[o]->energy <= 0 || orgsVector[o]->age > Org::maxAge)
		{
			Org* org = orgsVector[o];
			//fprintf(necroLog, "%d\t%d\t%f\n", id, org->age, org->meanFit());
			orgs(org->x, org->y) = NULL;
			orgsVector.erase(orgsVector.begin() + o);
			org->alive = false;
			if(bPhyloLog) org->maybeDelete(); else delete org;
		}
		else o--;
	}
}


void Mir::orgDivide()
{
	int count = orgsVector.size();
	if(count < 1) return;
	vector< lite::array<int[2]> > places; // free nearby cells

	for(int o = 0; o < orgsVector.size(); o++)
	{
		Org* org = orgsVector[o];
		if(org->energy > energyDivide && org->age > minAgeToDivide)
 		{
			// place
			places.clear();
			places.reserve(8);
			int x = org->x, y = org->y, X, Y;
			for(int d_x = -1; d_x <= 1; d_x++)
			{
                for(int d_y = -1; d_y <= 1; d_y++)
                {
                    if(d_x == 0 && d_y == 0) continue;
                    X = x + d_x; Y = y + d_y;
                    putToWorld(X, Y);
                    lite::array<int[2]> tmp(X,Y);
                    if(orgs(X,Y) == NULL)
                        places.push_back(tmp);
                }
			}
			if(places.size() == 0) continue;
			lite::array<int[2]> place = places[rand()%places.size()];
			//if(orgs(place[0],place[1]) != NULL) continue;
			Org* newbie = org->divide();
			newbie->x = place[0];
			newbie->y = place[1];
			// enzymes
			for(int g = 0; g < newbie->genome.size(); g++)
			{
                determineEnzyme(newbie->genome[g]);
			}
			orgs(newbie->x,newbie->y) = newbie;
			orgsVector.push_back(newbie);
			//if(divideLogOn) fprintf(divideLog, "%d -> %d;\n%d -> %d;\n", org->id, newbie->id, org->id, orgNewId);
	    }
	}
}


void Mir::putToWorld(int &x, int &y)
{
	if(x >= w) x -= w;
	if(y >= h) y -= h;

	if(y < 0) y += h;
//	if(x >= w) x = w-1;
//	if(y >= h) y = h-1;
//	if(x < 0) x = 0;
//	if(y < 0) y = 0;
}

//////////////////////////////////////////////////////////////////////////
// stat helpers

float Mir::meanEnzymeFit()
{
    if(orgsVector.size() == 0) return 0;
	float meanFit = 0;
	int count = 0;
	for(int o = 0; o < orgsVector.size(); o++)
	{
		for(int g = 0; g < orgsVector[o]->genome.size(); g++)
		{
			meanFit +=  orgsVector[o]->genome[g].fit;
			count++;
		}
	}
	return meanFit/count;
}

void Mir::saveGenomes()
{
    char fname[256];
	sprintf(fname, "MirAge_%d.fasta", age);
	FILE* fw = fopen(fname, "w");

	for(int o = 0; o < orgsVector.size(); o++)
	{
		for(int g = 0; g < orgsVector[o]->genome.size(); g++)
		{
			fprintf(fw, ">%s|gene%d\n", orgsVector[o]->name.c_str() , g);
			fprintf(fw, "%s\n", orgsVector[o]->genome[g].seq.c_str());
		}
	}
	fclose(fw);
}



void Mir::giveNames(Org* org)
{
	for(int i = 0; i < org->children.size(); i++)
	{
		std::stringstream ss;
		ss << org->name  << "_" << i;
		org->children[i]->name = ss.str();	
		giveNames(org->children[i]);
	}
}

bool Mir::badSubstance(int id)
{
	// energetically non profit to eat
	for(int i = 0; i < NSubstances; i++)
		if(dE(id, i) > 0) return false;
	return true;
}

void Mir::openLogFiles()
{
	// 4 future maybe close-open logs from time to time 4 monitoring
	//necroLog = fopen("necroLog", "w");
	//divideLog = fopen("divideLog.dot", "w");
	populationLog = fopen(popLogFile, "w");
	if(!populationLog) { cout << "problem creating logs... :(\n"; return; }
	//fprintf(necroLog, "id\tage\tfit\n");
	//fprintf(divideLog, "id\tfit1\tfit2\tage\tenergy\n");
	//fprintf(divideLog, "graph History {\n");
	fprintf(populationLog, "id\torgs\tmeanFit\tmedianFit\tmaxFit\n");
}

void Mir::closeLogFiles()
{
	//fclose(necroLog);
	//fprintf(divideLog, "}\n");
	//fclose(divideLog);
	fclose(populationLog);
}

void Mir::logPopulation()
{
	if((age % LOG_FREQ) != 0 || orgsVector.size() == 0) return;
	float meanFit = 0;
	vector<float> fits; fits.resize(orgsVector.size());
	for(int o = 0; o < orgsVector.size(); o++)
	{
		meanFit += orgsVector[o]->meanFit();
		fits[o] = orgsVector[o]->meanFit();
	}
	sort(fits.begin(), fits.end());
	float medianFit = fits[orgsVector.size()/2];

	meanFit /= orgsVector.size();
	fprintf(populationLog, "%d\t%d\t%f\t%f\t%f\n", id, (int)orgsVector.size(), meanFit, medianFit, fits[orgsVector.size()-1]);
	fflush(populationLog);
	if (bSaveGenomes) saveGenomes();
}
