# Intro

Let's be honest, stringmol needs refactoring. There's too much redundant code in there, it's hard to configure, and there are holes in the functionality. But beyond that, we need a tight code base if we are going to build on what's gone before. 

Although there are many 'standard' ways to refactor, I've not found any of them to be easy to understand and to get started. So this document is going to record what I've tried, what works and what doesn't. 


The first thing I did was ask twitter! 

# Current activity

- 200811: see github.com/franticspider/TryingTravis for how I'm developing a continuous integration environment. 
- Decision: Clean up the Stringmol "reaction" module first. This is used stand-alone in Rstringmol, the web app and as the basis for Centurion

# Roadmap

Plan is to have iterations through the following:

- Identify a module
- Document dependencies
- Write tests

## Features to add/fix

- stringmol 1 on 1 needs to work for any pair of molcules, not just 2 copies of the same one. 
- Resume via logfiles needs to be thoroughly tested
- zero-probability bind issue needs to be addressed


## File by file

### `stringmol.cpp`

There's a *huge* amount of redundancy in this - many functions are minor tweaks of other functions and it isn't clear which ones were ever used. As an example, look for repetition of the variable `divct` in various functions... so a big refactor of the major code blocks within these functions is in order. 


# Background

## Raw Twitter comments

Thanks here to David White (DW), Jerry Swan(JS), John Tuffen(JT) and Emily Dolson (ED) for these suggestions

- DW: Renaming things so they explain what they are
- DW: Cleanup build process if applies and write a http://README.md ... get someone to try following it.
- DW: Organising into modules / packages
- DW: Remove code that is commented out (obv everything is in git anyway right)
- JS: What's the language? The support one can get from the compiler makes some difference to best strategies...
- JS: Language-independent answer - the 'facade' pattern is one of the best general strategies:
https://en.wikipedia.org/wiki/Facade_pattern
- DW: Replace bespoke code with standard libraries as much as possible. Easy way to remove code.
- DW: Remove all generated, temporary, unrelated, files from the repo.
- DW: Use an IDE eg jetbrains to move stuff around and rename (refactor) without breaking anything
- JS: Approximate ordering for large-scale refactoring: 
	1. Delete/replace code.
	2. Modularize. Perhaps first by simply moving source code around so that related functionality is grouped, then  start to channel communication via (preferably contractual) APIs as per the Facade pattern.
- JT: Unit testing of interface APIs first and foremost — subsequent bugs introduced by refactoring will then be found more quickly; DW: I’m assuming any notion of an “API” is a pipe dream here, but if you’ve got something more polished then agree you should definitely write unit tests and write a wee tutorial on how to use it before changing the code; JT: I supppse it also depends on where do you want to win? If it’s  a performance issue then some analysis is probably worth doing to find hotspots, rather than spending time optimising the ‘wrong things; SH: This project is going in several directions at once (R package; Web app; dedicated Hardware; GPU) and if I don't sort this now it's going to bite me. Speed is a huge issue but I know where that particular problem is and I'll fight that dragon another day; JS: Ouch. Encapsulation and layering behind a succession of facades will help prevent the codebase fragmenting in a horribly combinatorial way.
- ED: Low effort: Document the dependencies! Write down any assumptions you're making about how stuff is set up. Make sure there's enough info that you could get it running again a year from now.
- ED: Medium effort: Document the expected output of the code under various circumstances. Write a script that runs the code and confirms that it produces that output. Ta-da you've got basic regression tests!
- ED: Low effort (once you've done the above step and assuming your code is already shared on github or similar): Set up continuous integration (e.g. @travisci) to automatically run  the tests when you make changes. That way you'll know if you do anything that changes behavior.

## Jerry Swan's email

0. Have some regression tests, so that you can be sure that refactoring doesn't break anything. Run them very often (hence prefer that refactorings are small).

1. Most importantly: the integrity of a program is inversely proportional to the amount of mutable state in scope. Hence, see what can readily be done to reduce the visibility of variables, starting with globals, then file-scope statics etc. There's a section on this in the attached.

2. Following on from 1.: prefer to encapsulate data within classes - don't let outside entities tweak class attributes individually via 'setter' methods. Instead, the class should instead be providing higher level services via its methods, and not just be 'a dumb bucket of bits' to be pushed around by the rest of the program. Use const on methods and parameters wherever possible.

3. Start trying to get more help from the compiler, by 'programming in the vocabulary of the problem domain' (again, see attached). Essentially, this comes down using the type system to distinguish between incompatible types: if Speed and Length are incomensurate but are both currently represented as int, then at least use typedefs to explicitly distinguish them, and maybe even make them different types.  

4. When passing arguments by reference in C++, prefer references to pointers wherever possible: barring 'fraud' (deliberately dereferencing a potentially null pointer to make a reference) there is no such thing as a 'null reference'.

5. Start to build integrity into the program via 'design by contract' (see attached), particularly class invariants. This a) affords greater code coverage than unit tests and b) forces one to be very clear about what concept a class actually represents.

## Status before writing this document

- I'd got stringmol on github
- We had a 'lost' web app at stringmol.york.ac.uk
- I'd developed Rstringmol for a paper
- I was working on deploying stringmol to Centurion hardware
- I'd written some basic tests but hadn't maintained them
- A tech. report from 2012, which describes functionality but needs updating
- some runs through valgrind to spot any memory leaks
- an iniitial attempt at documenting the code using doxygen was made
