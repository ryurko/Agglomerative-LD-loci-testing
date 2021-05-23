# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0045202","synapse",0.487,9.834,0.639,0.000,"synapse"),
c("GO:0099092","postsynaptic density, intracellular component",0.001,3.258,0.596,0.553,"synapse"),
c("GO:0060076","excitatory synapse",0.007,2.972,0.682,0.633,"synapse"),
c("GO:0099091","postsynaptic specialization, intracellular component",0.001,3.337,0.592,0.650,"synapse"),
c("GO:0098982","GABA-ergic synapse",0.011,4.147,0.675,0.652,"synapse"),
c("GO:1902494","catalytic complex",3.457,5.732,0.839,0.000,"catalytic complex"),
c("GO:0042611","MHC protein complex",0.015,4.566,0.789,0.307,"catalytic complex"),
c("GO:0016342","catenin complex",0.006,2.985,0.799,0.445,"catalytic complex"),
c("GO:0098878","neurotransmitter receptor complex",0.009,2.630,0.795,0.454,"catalytic complex"),
c("GO:0042613","MHC class II protein complex",0.012,2.634,0.792,0.461,"catalytic complex"),
c("GO:1990234","transferase complex",1.442,5.546,0.816,0.491,"catalytic complex"),
c("GO:0098797","plasma membrane protein complex",1.538,3.289,0.725,0.504,"catalytic complex"),
c("GO:0098796","membrane protein complex",3.375,4.628,0.824,0.553,"catalytic complex"),
c("GO:0034708","methyltransferase complex",0.142,4.290,0.784,0.602,"catalytic complex"),
c("GO:0034703","cation channel complex",0.171,3.184,0.827,0.621,"catalytic complex"),
c("GO:1905368","peptidase complex",0.561,2.658,0.829,0.690,"catalytic complex"),
c("GO:0036477","somatodendritic compartment",0.153,5.183,1.000,0.000,"somatodendritic compartment"),
c("GO:0016604","nuclear body",0.231,6.619,0.794,0.000,"nuclear body"),
c("GO:0098576","lumenal side of membrane",0.006,6.163,0.840,0.169,"nuclear body"),
c("GO:0015630","microtubule cytoskeleton",1.054,6.032,0.771,0.197,"nuclear body"),
c("GO:0012507","ER to Golgi transport vesicle membrane",0.074,4.994,0.690,0.355,"nuclear body"),
c("GO:0000785","chromatin",0.550,2.801,0.811,0.458,"nuclear body"),
c("GO:0005794","Golgi apparatus",1.491,3.348,0.769,0.534,"nuclear body"),
c("GO:0005802","trans-Golgi network",0.112,2.582,0.807,0.554,"nuclear body"),
c("GO:0000124","SAGA complex",0.060,2.563,0.692,0.618,"nuclear body"),
c("GO:0005798","Golgi-associated vesicle",0.115,4.090,0.741,0.621,"nuclear body"),
c("GO:0035097","histone methyltransferase complex",0.090,5.020,0.684,0.637,"nuclear body"),
c("GO:0016607","nuclear speck",0.127,3.547,0.802,0.654,"nuclear body"),
c("GO:0030135","coated vesicle",0.276,3.693,0.725,0.686,"nuclear body"),
c("GO:0043005","neuron projection",0.279,9.306,0.876,0.000,"neuron projection"),
c("GO:0032589","neuron projection membrane",0.010,2.571,0.684,0.644,"neuron projection"),
c("GO:1904115","axon cytoplasm",0.011,3.059,0.809,0.648,"neuron projection"),
c("GO:0098590","plasma membrane region",0.523,2.725,0.877,0.000,"plasma membrane region"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_cc_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

