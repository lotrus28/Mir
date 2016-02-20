#include <stdio.h>
#include "genomeDistanceLoger.h"
#include "mir.h"

float genomeDistanceLoger::calcDistance(Org* a, Org* b)
{
	float sumDistance = 0;
	int count = 0;
	for(int ai = 0; ai < a->genome.size(); ai++)
		for(int bi = 0; bi < b->genome.size(); bi++)
		{
			if(a->genome[ai].in == b->genome[bi].in &&
			  a->genome[ai].out == b->genome[bi].out)
			{
				sumDistance += calcDistance(a->genome[ai], b->genome[bi]);
				count++;
			}
		}
	return sumDistance/count;
}

float genomeDistanceLoger::calcDistance(Gene a,  Gene b)
{
	int id = 0;
	int geneLength = (int)a.seq.size();
	if((int)b.seq.size() != geneLength) { fprintf( stderr, "gene length differ! returning -1"); }
	for(int i = 0; i < geneLength; i++)
	{
		if(a.seq[i] == b.seq[i]) id++;
	}
	return id/geneLength;
}
