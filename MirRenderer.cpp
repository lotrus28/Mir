#include "MirRenderer.h"


MirRenderer::MirRenderer(Mir* mir)
{
	this->mir = mir;
	chosen = NULL;
	w = 1000; h = 1000;
	stepX = w/mir->w; stepY = h/mir->h;
	renderMode = 1;
	showSubstances = true;

	// stat window
	windowStat = SDL_CreateWindow("Stat", SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED, 400, 400, SDL_WINDOW_SHOWN);
	rendererStat  =  SDL_CreateRenderer(windowStat, -1, SDL_RENDERER_ACCELERATED);
//	SDL_SetRenderDrawColor(rendererStat, 245, 245, 220, 255);
	SDL_SetRenderDrawColor(rendererStat, 20, 20, 20, 255);
	SDL_RenderClear(rendererStat);


	// mir window
	window = NULL;
	window = SDL_CreateWindow("Mir", SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED, w, h, SDL_WINDOW_SHOWN);
	renderer = NULL;
	renderer  =  SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);

	// try populate colors
	substanceColors.resize(mir->NSubstances,3);
}

MirRenderer::~MirRenderer()
{
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}

//////////////////////////////////////////////////////////////////////////

void MirRenderer::main()
{
    srand(mir->randomSeed);
    populateColors();
	mir->init();
	bool quit = false;
	SDL_Event e;
    float i = 0;
	int statTimer = 0;
	while(!quit)
	{
		mir->tic();
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);
		if(renderMode!= 3)
		{
			if(showSubstances) drawSubstances();
			drawOrgs();
			drawStat();
			SDL_RenderPresent(rendererStat);
			SDL_RenderPresent(renderer);
		}
		int mx, my;
		while( SDL_PollEvent( &e ) != 0 )
		{
			if( e.type == SDL_QUIT ) { quit = true; }
			if( e.type == SDL_KEYUP)
			{
				if (e.key.keysym.sym == SDLK_s) {
					showSubstances = !showSubstances;
				}

				if (e.key.keysym.sym == SDLK_SPACE) {
					mir->degrade();
					mir->loadConfig();
					mir->init();
					populateColors();
				}
				if (e.key.keysym.sym == SDLK_m) {
					renderMode++;
					if(renderMode == 4) renderMode = 0;
				}

				if (e.key.keysym.sym == SDLK_q) {
					quit = true;
				}
			}
			if(SDL_GetMouseState(&mx, &my)&SDL_BUTTON(1))
			{
				mx /= stepX;
				my /= stepY;
				chosen = mir->org(mx, my);
			}
		}
		// stat
		if(++statTimer >= 500)
		{
			printf("mean fit: %f\torgs: %d\tmirAge:%d\n", mir->meanEnzymeFit(), mir->orgsVector.size(), mir->age);
			statTimer = 0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////

void MirRenderer::drawSubstances()
{
	for(int i = 0; i < mir->w; i++)
		for(int j = 0; j < mir->h; j++)
		{
			int x1, x2, y1, y2;
			x1 = stepX*i; x2 = x1 + stepX;
			y1 = stepY*j; y2 = y1 + stepY;
			//renderer
			lite::array<float[3]> color;
			color = 0;
			for(int s = 0; s < mir->NSubstances; s++)
			{
				float intensity = mir->substances(i,j,s)/Mir::maxSubstance;
				lite::array<float[3]> scolor = substanceColors[row(s)];
				color = color + scolor*intensity;
				for(int c = 0; c < 3; c++) if(color(c) > 255) {color(c) = 255;}

			}
			boxRGBA(renderer, x1, y1, x2, y2, color(0), color(1), color(2), 255);

			char str[5];
		}
}

void MirRenderer::drawOrgs()
{
	for(int o = 0; o < mir->orgsVector.size(); o++)
	{
		Org *org = mir->orgsVector[o];
		int x = (org->x+0.5)*stepX;
		int y = (org->y+0.5)*stepY;

		int r = 255, g = 255, b = 255;

		if(renderMode == 1 || renderMode == 2)
		{
			float fit = org->genome[0].fit;
			if(renderMode == 2) fit = 1;
			r = substanceColors(org->genome[0].in,0)*fit;
			g = substanceColors(org->genome[0].in,1)*fit;
			b = substanceColors(org->genome[0].in,2)*fit;
		}


		//filledCircleRGBA(renderer, x, y, min(stepX, stepY)*0.5, 0, 0, 0, 255);
		filledCircleRGBA(renderer, x, y, min(stepX, stepY)*0.4, r, g, b, 255);
		//filledCircleRGBA(renderer, x, y, 1, 255, 255, 255, 255);

	}
}

void MirRenderer::drawStat()
{
	try{
		if(!chosen) return;
		if(chosen->energy < mir->expressionCost*3) { chosen = NULL; return; }
		SDL_RenderClear(rendererStat);
		Uint32 color = 0xFFFFF0FF;
		char str1[1024], str2[1024];
		sprintf(str1,"Energy: %f    Age: %d", chosen->energy, chosen->age);
		sprintf(str2,"Gene0 fit = %f", chosen->genome[0].fit);

		stringColor(rendererStat, 10, 10, str1, color );
		stringColor(rendererStat, 10, 30, str2, color );
	} catch(...) {}
}

//////////////////////////////////////////////////////////////////////////

void MirRenderer::populateColors()
{
	for(int s = 0; s < mir->NSubstances; s++)
	{
		substanceColors(s,0) = rand()%255;
		substanceColors(s,1) = rand()%255;
		substanceColors(s,2) = rand()%255;
	}
}
