#ifndef DLOGER_H
#define DLOGER_H

class Mir;

class genomeDistanceLoger
{
public:
	genomeDistanceLoger(Mir* mir) { this->mir = mir; }
	~genomeDistanceLoger() {} 

private:
	Mir* mir;

};



#endif // DLOGER_H
