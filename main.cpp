#include "mir.h"
#include "MirRenderer.h"

#include "stdio.h"
#include "litearray.hpp"
#include <iostream>
#include <string>
#include <time.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL2_gfxPrimitives.h>
//#include "SDL2/SDL.h"
//#include "SDL2/SDL2_gfxPrimitives.h"
#include <random>

int main(int argc, char **argv)
{
    int rand;
	Mir m1(argc, argv);//hello*/
    MirRenderer mr1(&m1);
    mr1.main();
    return 0;
}
