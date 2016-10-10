
More recently appearing species appear first, (unless they are 'bellybuttonless' - see below). 

Each unique 'parenting' event is listed for each species, to attempt to trace multiple inheritance paths. 


An entry for a species in a splist might look like this: 

```
159722,134932,154542,0.000000,1,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,154752,0.000000,2,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,132361,0.000000,1,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,159722,0.000000,3,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,134932,0.000000,3,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,150129,0.000000,4,999516,12,>$$B^>C$=?>$$BFB%}P
159722,159722,151674,0.000000,2,999516,12,>$$B^>C$=?>$$BFB%}P
```

This is comma-delimited data, where each column indicates the following:

1. The species number
- The active parent species number
- The passive parent species number
- **(UNUSED)** The 'signal' score for a species
- The number of times parentage in this way has happened
- The **first** time this reaction happened
- **(UNUSED)** The 'biomass' of the sequence
- The seqeunce itself

###Caveats

Species that are 'seeded' from a config file are indicated thus:

```
146064,-1,-1,1,$=?>$O
146064,114271,146064,0.000000,63,0,0,$=?>$O
146064,117948,146064,0.000000,145,0,0,$=?>$O
146064,122296,146064,0.000000,83,0,0,$=?>$O
146064,146080,146064,0.000000,17,0,0,$=?>$O
146064,155683,146064,0.000000,2,0,0,$=?>$O
146064,144770,146064,0.000000,8,0,0,$=?>$O
146064,159519,146064,0.000000,2,0,0,$=?>$O
```

The first line is unique - 

Ancestry is *assumed* to be 2-parents-one-offspring, with parents unchanged in the reaction. In stringmol, this isn't always the case! Care needs to be taken if we are interested in what happens to parents *after* a cleave event. 
