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

//#define USE_SDL



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

// Stringmol
#include "alignment.h"

// Metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

// Writing PNGs
#include "lodepng.h"
#include <iostream>
#include "setupSM.h"

typedef enum td_pic{
	pic_spp,
	pic_len,
	pic_bound
}pictype;









int main(int argc, char *argv[]) {
//	int smspatial(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	int rtype = atoi(argv[1]);

	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;
	run = NULL;

	smspatial_init(argv[2],&A,&run,1);

	int bt=0,ct=0;
	ct = A.nagents(A.nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	/* Set up SDL if we're using it */
#ifdef USE_SDL

	int **grid;
	grid = (int **)malloc(run->gridx*sizeof(int *));
	for(int i=0;i<run->gridx;i++){
		grid[i] = (int *) malloc(run->gridy*sizeof(int));
		memset(grid[i],0,run->gridy*sizeof(int));
	}

	printf("grid allocated\n");fflush(stdout);

	SDL_Surface *screen;
	SDL_Event sdlEvent;
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		fprintf(stderr,"*** Unable to init SDL: %s ***\n",SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	SDL_WM_SetCaption("Spatial Stringmol","nanopond");
	screen = SDL_SetVideoMode(run->gridy,run->gridx,8,SDL_SWSURFACE);
	if (!screen) {
		fprintf(stderr, "*** Unable to create SDL window: %s ***\n", SDL_GetError());
		exit(1);
	}
	const uintptr_t sdlPitch = screen->pitch;
	printf("SDL window allocated, pitch is %d \n",screen->pitch);fflush(stdout);
#endif /* USE_SDL */


//	while(A.nagents(A.nowhead,-1)){
	while(A.extit < 1000000){// && A.nagents(A.nowhead,-1)>0){

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

//		if(A.extit == 66396)
//			printf("Pauuuuse\n!");
		A.extit++;
		//if(!(A.extit%100))
		//		printf("Step %d done, number of molecules is %d, nbound = %d\n",(int) A.extit,ct,bt);
		if((!(A.extit%100))  ){// ||     A.extit == 66396    ){
			printf("Step %d done, number of molecules is %d, nbound = %d\n",(int) A.extit,ct,bt);
			FILE *fp;char fn[128];
			sprintf(fn,"splist%d.dat",A.extit);
			fp = fopen(fn,"w");
			SP.print_spp_list(fp);
			fclose(fp);

			printsppct(&A,A.extit);

			//printf("Step %ld done, number of molecules is %d, nbound = %d\n",A.extit,ct,bt);

			FILE *fpp;
			sprintf(fn,"out1_%05ld.conf",A.extit);
			fpp = fopen(fn,"w");
			A.print_conf(fpp);
			fclose(fpp);

			//Now let's reload the file so we can check:
			/*
			SMspp		SPB;
			stringPM	B(&SPB);
			B.load(fn,NULL,0,1);
			B.print_agents(stdout,"NOW",0);
			sprintf(fn,"outB_%05d.conf",A.extit);
			fpp = fopen(fn,"w");
			B.print_conf(fpp);
			fclose(fpp);
			*/
		}

		//if(!(A.extit%10)){
		//	printf("completed iteration %u of simulation\n",A.extit);fflush(stdout);
		//	//printf("grid dimenstions are %d,%d\n",run->gridx,run->gridy);fflush(stdout);
		//}

		//Create the PNG
		std::vector<unsigned char> image;
		image.resize(run->gridx * run->gridy * 4);


	    pictype tp = pic_len;
		int x,y,val;
		for(x=0;x<run->gridx;++x){
			for (y=0;y<run->gridy;++y) {
				switch(run->status[x][y]){
				case G_EMPTY:
					val=0;
					break;
				default:
					switch(tp){
					case pic_spp:
						val=run->grid[x][y]->spp->spp * 50;
						break;
					case pic_len:
						val=strlen(run->grid[x][y]->spp->S) * 50;
						break;
					}
					/*
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
					}*/
				}

#ifdef USE_SDL
				((uint8_t *)screen->pixels)[y + (x * sdlPitch)] = val;//getColor(grid[x][y]);

				if(!(A.extit%10)){
					/*Get the rgb values for writing to PNG*/
					uint8_t r ;
					uint8_t g ;
					uint8_t b ;
					SDL_GetRGB(((uint8_t *)screen->pixels)[y + (x * sdlPitch)], screen->format ,  &r, &g, &b );
					//printf("R is %d G is %d B is %d\n",r,g,b);
					image[4 * run->gridx * y + 4 * x + 0] = r;//255 * !(x & y);
					image[4 * run->gridx * y + 4 * x + 1] = g;//x ^ y;
					image[4 * run->gridx * y + 4 * x + 2] = b;//x | y;
					image[4 * run->gridx * y + 4 * x + 3] = 255;
				}
#endif
			}
		}


#ifdef USE_SDL
		if(!(A.extit%10)){
			char filename[128];
			sprintf(filename,"frame%07d.png",A.extit);
			encodeOneStep(filename, image, run->gridx, run->gridy);
		}



		while (SDL_PollEvent(&sdlEvent)) {
			if (sdlEvent.type == SDL_QUIT) {
				fprintf(stderr,"[QUIT] Quit signal received!\n");

				//TODO: Write the final splist and popdy files here!

				exit(0);
			}
		}

		//printf("sdl grid updated\n",A.extit);fflush(stdout);

        SDL_UpdateRect(screen,0,0,0,0);//run->gridx,run->gridy);

		//printf("sdl rect updated\n",A.extit);fflush(stdout);
#endif

	}

	printf("FINISHED smspatial\n");
//	fflush(stdout);
	return 0;
}


