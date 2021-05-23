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
revigo.data <- rbind(c("GO:0003713","transcription coactivator activity",0.085,3.563,0.961,0.000,"transcription coactivator activity"),
c("GO:0060589","nucleoside-triphosphatase regulator activity",0.600,3.249,0.968,0.187,"transcription coactivator activity"),
c("GO:0090079","translation regulator activity, nucleic acid binding",0.849,3.311,0.866,0.325,"transcription coactivator activity"),
c("GO:0003723","RNA binding",5.649,7.092,0.888,0.349,"transcription coactivator activity"),
c("GO:0004672","protein kinase activity",3.471,4.175,0.809,0.000,"protein kinase activity"),
c("GO:0008170","N-methyltransferase activity",0.445,3.288,0.915,0.272,"protein kinase activity"),
c("GO:0008239","dipeptidyl-peptidase activity",0.047,4.100,0.940,0.367,"protein kinase activity"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.428,3.322,0.885,0.409,"protein kinase activity"),
c("GO:0008092","cytoskeletal protein binding",0.859,4.889,0.685,0.000,"cytoskeletal protein binding"),
c("GO:0045296","cadherin binding",0.018,3.689,0.753,0.445,"cytoskeletal protein binding"),
c("GO:0008022","protein C-terminus binding",0.026,3.170,0.747,0.459,"cytoskeletal protein binding"),
c("GO:0050839","cell adhesion molecule binding",0.052,3.229,0.737,0.484,"cytoskeletal protein binding"),
c("GO:0042393","histone binding",0.084,3.280,0.729,0.504,"cytoskeletal protein binding"),
c("GO:0051020","GTPase binding",0.168,4.439,0.677,0.535,"cytoskeletal protein binding"),
c("GO:0003779","actin binding",0.451,3.272,0.698,0.587,"cytoskeletal protein binding"),
c("GO:0098918","structural constituent of synapse",0.002,3.917,0.998,0.007,"structural constituent of synapse"),
c("GO:0005244","voltage-gated ion channel activity",0.253,3.583,0.925,0.010,"voltage-gated ion channel activity"),
c("GO:0042605","peptide antigen binding",0.001,3.159,0.962,0.034,"peptide antigen binding"),
c("GO:0008144","drug binding",0.154,4.737,0.948,0.049,"drug binding"),
c("GO:0003682","chromatin binding",0.187,3.416,0.947,0.050,"chromatin binding"),
c("GO:0032553","ribonucleotide binding",16.281,7.034,0.862,0.083,"ribonucleotide binding"),
c("GO:0030554","adenyl nucleotide binding",14.085,6.821,0.864,0.683,"ribonucleotide binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_mf_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

