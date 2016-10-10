#Rgraphviz is at bioconductor
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")

require("Rgraphviz")

#Following this at: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/rgraphviz/
#> library(Rgraphviz)
#> test.matrix<-matrix(rep(c(0,1,0,0), 9), ncol=6, nrow=6)
#> rownames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
#> colnames(test.matrix)<-c("a", "b", "c", "d", "e", "f")
#> test.matrix
#  a b c d e f
#a 0 0 0 0 0 0
#b 1 0 1 0 1 0
#c 0 0 0 0 0 0
#d 0 1 0 1 0 1
#e 0 0 0 0 0 0
#f 1 0 1 0 1 0
#> am.graph<-new("graphAM", adjMat=test.matrix, edgemode="directed")
#> am.graph
#A graphAM graph with directed edges
#Number of Nodes = 6 
#Number of Edges = 9 
#> plot(am.graph, attrs = list(node = list(fillcolor = "lightblue"),
#                                edge = list(arrowsize=0.5)))






fn<-"~/Desktop/ancestry121.txt"


#######make_tree <- function(fn){


	#x<- read.table(fn)
	x <- read.table(fn,header=T,sep=",")

	spnames <- unique(x$spp)
	#spnames <- unique(x$seq)
	

	x$depth <- as.integer(x$depth)
	x$time <- as.integer(x$time)
	x$spp <- as.integer(x$spp)
	x$act<- as.integer(x$act)
	x$pass<-as.integer(x$pass)
	

	#get a list of unique species
	spp <- unique(x$spp)
	
	
	nspp<- length(spp)

	#Create an adjacency table....
	
	adj.matrix<-matrix(rep(0,nspp*nspp), ncol = nspp, nrow = nspp)
	
	rownames(adj.matrix)<-spnames
	colnames(adj.matrix)<-spnames

	for(i in 1:nrow(x)){
		message(sprintf("Linking %d to active %d and passive %d",x$spp[i],x$act[i],x$pass[i]))
		
		rowno = match(x$spp[i],spp)
		colac = match(x$act[i],spp)
		colpa = match(x$pass[i],spp)
		
		adj.matrix[colac,rowno] = 1
		adj.matrix[colpa,rowno] = 1
	}
	
	gr <- new("graphAM", adjMat=adj.matrix, edgemode="directed")
	plot(gr, attrs = list(node = list(fillcolor = "lightblue"),edge = list(arrowsize=0.5)))

##}

	#Try again with an order: 




