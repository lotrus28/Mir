#ifndef MIR_H
#define MIR_H

#include <vector>
#include <string>
#include <deque>
#include "litearray.hpp"
#include "genomeDistanceLoger.h"
using namespace std;


//==========================================================

class Mir;

struct Gene
{
	Gene() { };
	string seq;
	int in, out; // its enzyme
	float fit; // how good can do reaction
};

class Soul
{
public:
	Soul();
	~Soul() {};
	string name;
	bool alive;
	deque<Soul*> children;
	Soul* parent;
public:
	void die();
	bool anyLivingChild();
	void maybeDelete();
	void deleteAll();
};


class Org
{
public:
	Org(Mir* mir);
	~Org();
public:
	Org* divide();
	float meanFit();
public:
	bool anyLivingChild();
	void maybeDelete();
public:
	vector<Gene> genome;
	string id;
	static int maxId;
	float energy;
	int age;
	int x, y;
	static float maxEAT;
	static float maxEnergy;
	static int maxAge;
	float SNPrate;
	Soul* soul;
private:
    Mir* mir;
};

//==========================================================

struct SubstanceSource
{
	int x, y, age;
	int radius;
	int substanceId;
	float intensity;
	SubstanceSource()
	{
		x = 0; y = 0; intensity = 1.0; radius = 0;
	}
};

//==========================================================

struct Cell
{
	int x, y; // needed?
	vector<float> substances;
	vector<Org*> orgs; // needed?
	void init(int NSubstances)
	{
		substances.resize(NSubstances);
	}
	friend ostream& operator<<(ostream& os, const Cell& cell)
	{
	 //os << (cell.orgs.size() > 0? 1 : 0);
	 return os;
	}

};

// ****************************************************
// class MIR

class Mir
{
public:
	Mir(int argc, char **argv);
	void init();
	void deinit();
	void loadConfig();
	~Mir();
public:
	void main();
	void tic();
public: // stat
	Org* org(int x, int y);
	float meanEnzymeFit();
	void saveGenomes();
	void giveNames(Soul* soul);
	void calcTaylorMeanVariance(float& mean, float& var);
private:
	void sourcesEmmit();
	void diffuse();
	void orgEat();
	void orgDie();
	void orgDivide();
	void sourceReincarnate();
	void logPopulation();

	// helpers
	void putToWorld(int &x, int &y);
	void nullPole();
	void populateOrgs(); // make their genomes
	void populate_dE();
	void populateGoldSeqs();
	void populateSources();
	void createRandomSource(SubstanceSource& source);
	void openLogFiles();
	void closeLogFiles();
	string randomSeq(int length = 10);
	void determineEnzyme(Gene& gene);
	bool badSubstance(int id);
	void reportDEMatrix();
// data
public:
	static int maxId;
	int id, age;
	int w,h;
	static const char alphabet[4];
	static const int alphabetLength = 4;
	int randomSeed;
public: // Pole
	lite::array <float[1][1][1]> substances;
	lite::array <Org*[1][1]> orgs;
	lite::array <vector<Gene>[1][1]> deadGenomes;
	vector<Org*> orgsVector;
	vector<SubstanceSource> sources;
	vector<int> goodSubstances; // energy profit
	Org* adam;
public: // Energetics
	int NSubstances;
	float minSubstance_dE, maxSubstance_dE;
	int minSourceRadius, maxSourceRadius;
	int sourceMaxIntensity;
	float substanceDegrade;
	lite::array<float[1][1]> dE;
	lite::array<string[1][1]> goldSeqs;
	float diffusion;
	static float maxSubstance;
	float expressionCost;
	float energyDivide;
	int minAgeToDivide;
	int SourceLifetime;
	int MirLifetime;
private: // at the beginning
	int orgStartCount;
	int initialOrgGenes;
	float initialOrgEnergy;
	int NSubstanceSources;
	int initialGeneLength;
	float initialSNPrate;
private: // Logs and files
	FILE *necroLog, *divideLog, *populationLog;
	bool divideLogOn;
	char *paramFile, *popLogFile, *constFile;
	bool bSaveGenomes;
public:
	bool bPhyloLog; 
};

#endif // MIR_H
