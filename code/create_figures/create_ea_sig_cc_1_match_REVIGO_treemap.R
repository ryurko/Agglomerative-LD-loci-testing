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
revigo.data <- rbind(c("GO:0043005","neuron projection",0.274,9.825,0.871,0.000,"neuron projection"),
c("GO:0032589","neuron projection membrane",0.011,2.611,0.673,0.649,"neuron projection"),
c("GO:1904115","axon cytoplasm",0.011,2.647,0.809,0.651,"neuron projection"),
c("GO:0098576","lumenal side of membrane",0.006,5.883,0.843,0.000,"lumenal side of membrane"),
c("GO:0016604","nuclear body",0.233,5.410,0.812,0.175,"lumenal side of membrane"),
c("GO:0015630","microtubule cytoskeleton",1.020,5.099,0.802,0.206,"lumenal side of membrane"),
c("GO:0031984","organelle subcompartment",2.118,2.561,0.813,0.278,"lumenal side of membrane"),
c("GO:0012507","ER to Golgi transport vesicle membrane",0.073,5.146,0.710,0.362,"lumenal side of membrane"),
c("GO:0005794","Golgi apparatus",1.462,2.707,0.780,0.539,"lumenal side of membrane"),
c("GO:0005798","Golgi-associated vesicle",0.112,4.324,0.757,0.623,"lumenal side of membrane"),
c("GO:0030660","Golgi-associated vesicle membrane",0.102,3.304,0.728,0.639,"lumenal side of membrane"),
c("GO:0035097","histone methyltransferase complex",0.092,5.256,0.684,0.639,"lumenal side of membrane"),
c("GO:0016607","nuclear speck",0.125,2.885,0.820,0.654,"lumenal side of membrane"),
c("GO:0000139","Golgi membrane",0.723,2.581,0.738,0.674,"lumenal side of membrane"),
c("GO:0030135","coated vesicle",0.272,4.048,0.742,0.688,"lumenal side of membrane"),
c("GO:0036477","somatodendritic compartment",0.149,5.262,1.000,0.000,"somatodendritic compartment"),
c("GO:0045202","synapse",0.471,9.756,0.630,0.000,"synapse"),
c("GO:0099092","postsynaptic density, intracellular component",0.001,2.966,0.598,0.555,"synapse"),
c("GO:0060076","excitatory synapse",0.008,2.966,0.671,0.637,"synapse"),
c("GO:0099091","postsynaptic specialization, intracellular component",0.001,2.626,0.595,0.652,"synapse"),
c("GO:0098982","GABA-ergic synapse",0.012,5.446,0.665,0.655,"synapse"),
c("GO:0048786","presynaptic active zone",0.019,2.558,0.652,0.677,"synapse"),
c("GO:0098796","membrane protein complex",3.372,5.216,0.813,0.025,"membrane protein complex"),
c("GO:0042611","MHC protein complex",0.013,4.602,0.790,0.305,"membrane protein complex"),
c("GO:0034708","methyltransferase complex",0.144,4.507,0.788,0.377,"membrane protein complex"),
c("GO:0044815","DNA packaging complex",0.223,2.764,0.864,0.395,"membrane protein complex"),
c("GO:0016342","catenin complex",0.006,3.018,0.799,0.441,"membrane protein complex"),
c("GO:0042613","MHC class II protein complex",0.011,2.654,0.793,0.454,"membrane protein complex"),
c("GO:1990351","transporter complex",1.292,2.692,0.843,0.485,"membrane protein complex"),
c("GO:0098797","plasma membrane protein complex",1.554,3.676,0.719,0.501,"membrane protein complex"),
c("GO:0031588","nucleotide-activated protein kinase complex",0.010,2.581,0.864,0.533,"membrane protein complex"),
c("GO:1990234","transferase complex",1.421,3.889,0.822,0.548,"membrane protein complex"),
c("GO:1902494","catalytic complex",4.505,3.942,0.824,0.578,"membrane protein complex"),
c("GO:0034703","cation channel complex",0.163,2.981,0.826,0.619,"membrane protein complex"),
c("GO:0098590","plasma membrane region",0.437,3.299,0.883,0.037,"plasma membrane region"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_sig_cc_1_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO Gene Ontology treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

