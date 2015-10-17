
#include "mir.h"

#include "stdio.h"
#include <iostream>



int main(int argc, char **argv)
{
	int rand;
	Mir m1(argc, argv);
	m1.init();
	m1.main();
	m1.degrade();
 	return 0;
}
