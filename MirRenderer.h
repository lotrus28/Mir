#ifndef MIR_RENDERER_H
#define MIR_RENDERER_H

#include <vector>
#include <queue>
#include "mir.h"
#include "SDL2/SDL.h"
#include "SDL2/SDL2_gfxPrimitives.h"

using namespace std;
using namespace lite;

//////////////////////////////////////////////////////////////////////////
// class MirRenderer
class MirRenderer
{
public:
	MirRenderer(Mir* mir);
	~MirRenderer(void);
public:
	void main();
private:
	void drawSubstances();
	void drawOrgs();
	void drawStat();
	void populateColors();
private:
	Mir* mir;
	int w,h;
	SDL_Window *window, *windowStat;
	SDL_Renderer *renderer, *rendererStat;
	int renderMode;
	float contrast, contrastSubstance;
	bool showSubstances;
	Org* chosen;
	int stepX, stepY;
	lite::array<float[1][3]> substanceColors;
	deque<float> fitVector;
	int statI;
};

#endif
