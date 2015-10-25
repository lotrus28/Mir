
#include "mir.h"

#include "stdio.h"
#include <iostream>



int main(int argc, char **argv)
{
	int rand;
	Mir m1(argc, argv);
	m1.init();
/*	Org* a = new Org(&m1);
	Org* c1 = a->divide();
	Org* c3 = a->divide();
	Org* c2 = c1->divide();
	a->name = "a";	
	m1.giveNames(a);
	c3->alive = false;
	c2->alive = false;
	c1->alive=false;

	c3->maybeDelete();
	c2->maybeDelete();
//	c1->maybeDelete();

	//delete a;
	a->alive = false;
	a->maybeDelete();
*/
	m1.main();
	m1.degrade();
 	return 0;
}
