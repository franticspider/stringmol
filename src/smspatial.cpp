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





int main(int argc, char *argv[]) {


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


    for (int y=0;y<300;++y) {
      for (int x=1;x<300;++x){

        ((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = 0;//getColor(grid[x][y]);
      }
      //((uint8_t *)screen->pixels)[i + (20 * sdlPitch)] = i;//getColor(grid[x][y]);
    }




	/*Let's try a dummy run now*/
	for(int i=0;i<60;i++){

        for (int y=0;y<300;++y) {
          ((uint8_t *)screen->pixels)[y + (20 * sdlPitch)] = i;//getColor(grid[x][y]);
        }




        SDL_UpdateRect(screen,0,0,300,300);

	    unsigned int retTime = time(0) + 1;   // Get finishing time.
	    while (time(0) < retTime);
	    printf("Waited %d seconds\n",i+1);
	}


}
