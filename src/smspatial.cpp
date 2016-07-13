/* Copyright (C) 2009-2016 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#define USE_SDL 1



#include <stdlib.h>
#include <stdio.h>

//Used for wait time...
#include <time.h>

#ifdef USE_SDL
#ifdef _MSC_VER
#include <SDL.h>
#else
#include <SDL/SDL.h>
#endif /* _MSC_VER */
#endif


#include "randutil.h"
#include "params.h"
#include "hsort.h"
//}

//stringmol
#include "alignment.h"


//metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"
#include "setupSM.h"




//TODO: This was just a test - see smspatial above for the new working version
int oldmain(int argc, char *argv[]) {


	printf("Hello spatial stringmol world\n");

	/*TODO: These must be read from a file:*/

	unsigned int GRID_X = 300;
	unsigned int GRID_Y = 300;


/* Set up SDL if we're using it */
#ifdef USE_SDL
	SDL_Surface *screen;
	SDL_Event sdlEvent;
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
	fprintf(stderr,"*** Unable to init SDL: %s ***\n",SDL_GetError());
	exit(1);
	}
	atexit(SDL_Quit);
	SDL_WM_SetCaption("Spatial Stringmol","nanopond");
	screen = SDL_SetVideoMode(GRID_X,GRID_Y,8,SDL_SWSURFACE);
	if (!screen) {
	fprintf(stderr, "*** Unable to create SDL window: %s ***\n", SDL_GetError());
	exit(1);
	}
	const uintptr_t sdlPitch = screen->pitch;
#endif /* USE_SDL */

	int grid[300][300];
	for(int i=0;i<300;i++){
		for(int j=0;j<300;j++){
			grid[i][j]=0;
		}
	}

	int x,y;

	//Clear the screen
    for (y=0;y<300;++y) {
      for (x=0;x<300;++x){

        ((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = 0;//getColor(grid[x][y]);
      }
      //((uint8_t *)screen->pixels)[i + (20 * sdlPitch)] = i;//getColor(grid[x][y]);
    }



    /*Let's try to init the grid:*/

    int nag = 150;
    for(int n=0;n<nag;n++){
    	int found = 0;
    	while(!found){
			int pos = GRID_X*GRID_Y*rand0to1();

			x = pos/GRID_X;
			y = pos%GRID_X;

			if(grid[x][y]==0){

		        ((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = 0;

		        //Add the partner in the Moore neighborhood
		        int ffound=0;

		        while(!ffound){
					int xx,yy;
					//randy_Moore(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout){
					randy_Moore(x,y,GRID_X,GRID_Y,&xx,&yy);
					if(grid[xx][yy]==0){
						ffound=found=1;
						grid[xx][yy] = grid[x][y] = 1;
					}
		        }
			}
    	}
    }



	/*Let's try a dummy run now*/
	for(int i=0;i<60;i++){


		for(x=0;x<300;++x){
			for (y=0;y<300;++y) {
				if(grid[x][y])
						((uint8_t *)screen->pixels)[y + (x * sdlPitch)] = 201;//getColor(grid[x][y]);
			}
		}
        SDL_UpdateRect(screen,0,0,300,300);

	    unsigned int retTime = time(0) + 1;   // Get finishing time.
	    while (time(0) < retTime);
	    printf("Waited %d seconds\n",i+1);
	}


}



int main(int argc, char *argv[]) {
//	int smspatial(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;
	run = NULL;

	smspatial_init(argv[2],&A,&run);

	int bt=0,ct=0;
	ct = A.nagents(A.nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	int x,y,iteration = 0;


	/* Set up SDL if we're using it */
#ifdef USE_SDL


	int **grid;
	grid = (int **)malloc(run->gridx*sizeof(int *));
	for(int i=0;i<x;i++){
		grid[i] = (int *) malloc(run->gridy*sizeof(int));
		memset(grid[i],0,run->gridy*sizeof(int));
	}

	SDL_Surface *screen;
	SDL_Event sdlEvent;
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
	fprintf(stderr,"*** Unable to init SDL: %s ***\n",SDL_GetError());
	exit(1);
	}
	atexit(SDL_Quit);
	SDL_WM_SetCaption("Spatial Stringmol","nanopond");
	screen = SDL_SetVideoMode(run->gridx,run->gridy,8,SDL_SWSURFACE);
	if (!screen) {
	fprintf(stderr, "*** Unable to create SDL window: %s ***\n", SDL_GetError());
	exit(1);
	}
	const uintptr_t sdlPitch = screen->pitch;
#endif /* USE_SDL */


//	while(A.nagents(A.nowhead,-1)){
	while(iteration < 10000){
		smspatial_step(&A,run);
		ct = A.nagents(A.nowhead,-1);
		bt = ct - A.nagents(A.nowhead,B_UNBOUND);
#ifdef DODEBUG
		printf("Nowhead is %d, Nexthead is %d\n",A.nowhead,A.nexthead);
		s_ag *p;
		p=A.nowhead;
		int mno=0;
		while(p!=NULL){
			int x,y;
			find_ag_gridpos(p,run,&x,&y);
			printf("%d, %d, [%d,%d]  status: %d, bound to %d / %d, prev = %d, next = %d\n",++mno,p,x,y,p->status,p->exec,p->pass,p->prev,p->next);
			p = p->next;
		}
#endif

//		if(iteration == 66396)
//			printf("Pauuuuse\n!");
		iteration++;
		if(!(iteration%100))
				printf("Step %d done, number of molecules is %d, nbound = %d\n",iteration,ct,bt);
		if((!(iteration%10000)) ||     iteration == 66396    ){
			printf("Step %d done, number of molecules is %d, nbound = %d\n",iteration,ct,bt);
			FILE *fp;char fn[128];
			sprintf(fn,"splist%d.dat",iteration);
			fp = fopen(fn,"w");
			SP.print_spp_list(fp);
			fclose(fp);
		}



		int x,y,val;
		for(x=0;x<run->gridx;++x){
			for (y=0;y<run->gridy;++y) {
				switch(run->status[x][y]){
				case G_EMPTY:
					val=0;
					break;
				default:
					switch(run->grid[x][y]->status){
					case B_UNBOUND:
						val = 100;
						break;
					case B_PASSIVE:
						val = 150;
						break;
					case B_ACTIVE:
						val = 200;
						break;
					default:
						val = 50;
					}
				}

				((uint8_t *)screen->pixels)[y + (x * sdlPitch)] = val;//getColor(grid[x][y]);
			}
		}
		while (SDL_PollEvent(&sdlEvent)) {
			if (sdlEvent.type == SDL_QUIT) {
				fprintf(stderr,"[QUIT] Quit signal received!\n");
				exit(0);
			}
		}


        SDL_UpdateRect(screen,0,0,300,300);

	}

	printf("FINISHED smspatial\n");
	fflush(stdout);
	return 0;
}


