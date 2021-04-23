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
revigo.data <- rbind(c("GO:0008144","drug binding",0.158,4.209,0.950,0.000,"drug binding"),
c("GO:0060589","nucleoside-triphosphatase regulator activity",0.557,3.625,1.000,0.000,"nucleoside-triphosphatase regulator activity"),
c("GO:0098918","structural constituent of synapse",0.002,3.943,0.997,0.000,"structural constituent of synapse"),
c("GO:0004672","protein kinase activity",3.474,3.836,0.701,0.005,"protein kinase activity"),
c("GO:0008170","N-methyltransferase activity",0.437,3.680,0.861,0.269,"protein kinase activity"),
c("GO:0008092","cytoskeletal protein binding",0.824,4.141,0.864,0.049,"cytoskeletal protein binding"),
c("GO:0051020","GTPase binding",0.058,4.125,0.877,0.487,"cytoskeletal protein binding"),
c("GO:0030554","adenyl nucleotide binding",14.174,5.767,0.805,0.082,"adenyl nucleotide binding"),
c("GO:0003723","RNA binding",5.619,5.199,0.883,0.235,"adenyl nucleotide binding"),
c("GO:0032553","ribonucleotide binding",16.376,5.740,0.803,0.689,"adenyl nucleotide binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_sig_mf_1_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

