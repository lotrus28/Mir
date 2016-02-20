#ifndef DLOGER_H
#define DLOGER_H


class Mir;
class Org;
class Gene;

class genomeDistanceLoger
{
public:
	genomeDistanceLoger(Mir* mir) { this->mir = mir; }
	~genomeDistanceLoger() {} 

private:
public:
	float	calcDistance(Org* a, Org* b);
	float	calcDistance(Gene a, Gene b);

private:
	Mir* mir;

};



#endif // DLOGER_H
