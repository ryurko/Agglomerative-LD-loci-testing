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
revigo.data <- rbind(c("GO:0005215","transporter activity",9.860,4.134,1.000,0.000,"transporter activity"),
c("GO:0008047","enzyme activator activity",0.640,3.770,1.000,0.000,"enzyme activator activity"),
c("GO:0008144","drug binding",0.154,10.606,0.961,0.000,"drug binding"),
c("GO:0015149","hexose transmembrane transporter activity",0.030,3.427,0.928,0.000,"hexose transmembrane transporter activity"),
c("GO:0015075","ion transmembrane transporter activity",3.173,3.273,0.904,0.364,"hexose transmembrane transporter activity"),
c("GO:0005385","zinc ion transmembrane transporter activity",0.018,2.567,0.924,0.506,"hexose transmembrane transporter activity"),
c("GO:0016887","ATPase",4.624,4.266,1.000,0.000,"ATPase"),
c("GO:0030234","enzyme regulator activity",0.977,6.073,0.952,0.000,"enzyme regulator activity"),
c("GO:0003712","transcription coregulator activity",0.292,5.651,0.970,0.218,"enzyme regulator activity"),
c("GO:0003713","transcription coactivator activity",0.085,4.275,0.963,0.452,"enzyme regulator activity"),
c("GO:0140098","catalytic activity, acting on RNA",3.383,6.405,0.972,0.000,"catalytic activity, acting on RNA"),
c("GO:0098918","structural constituent of synapse",0.002,2.511,0.997,0.005,"structural constituent of synapse"),
c("GO:0042605","peptide antigen binding",0.001,3.816,0.972,0.031,"peptide antigen binding"),
c("GO:0031491","nucleosome binding",0.026,3.597,0.931,0.037,"nucleosome binding"),
c("GO:0023026","MHC class II protein complex binding",0.001,3.040,0.944,0.493,"nucleosome binding"),
c("GO:0003954","NADH dehydrogenase activity",0.353,2.692,0.978,0.037,"NADH dehydrogenase activity"),
c("GO:0140101","catalytic activity, acting on a tRNA",1.721,3.611,0.938,0.045,"catalytic activity, acting on a tRNA"),
c("GO:0003724","RNA helicase activity",0.183,3.132,0.947,0.631,"catalytic activity, acting on a tRNA"),
c("GO:0042802","identical protein binding",0.404,8.567,0.750,0.046,"identical protein binding"),
c("GO:0036435","K48-linked polyubiquitin modification-dependent protein binding",0.004,2.728,0.815,0.378,"identical protein binding"),
c("GO:0051787","misfolded protein binding",0.012,2.525,0.802,0.409,"identical protein binding"),
c("GO:0045296","cadherin binding",0.018,7.070,0.789,0.421,"identical protein binding"),
c("GO:0030544","Hsp70 protein binding",0.033,2.692,0.789,0.441,"identical protein binding"),
c("GO:0050839","cell adhesion molecule binding",0.052,8.089,0.783,0.456,"identical protein binding"),
c("GO:0140030","modification-dependent protein binding",0.055,3.067,0.782,0.458,"identical protein binding"),
c("GO:0042393","histone binding",0.084,2.992,0.776,0.474,"identical protein binding"),
c("GO:0044389","ubiquitin-like protein ligase binding",0.094,6.724,0.739,0.478,"identical protein binding"),
c("GO:0051087","chaperone binding",0.104,3.324,0.773,0.481,"identical protein binding"),
c("GO:0031072","heat shock protein binding",0.109,3.245,0.772,0.483,"identical protein binding"),
c("GO:0019904","protein domain specific binding",0.114,6.944,0.771,0.485,"identical protein binding"),
c("GO:0046982","protein heterodimerization activity",0.275,5.145,0.753,0.523,"identical protein binding"),
c("GO:0015631","tubulin binding",0.357,5.063,0.747,0.535,"identical protein binding"),
c("GO:0008092","cytoskeletal protein binding",0.859,5.571,0.735,0.581,"identical protein binding"),
c("GO:0005102","signaling receptor binding",0.505,4.318,0.746,0.594,"identical protein binding"),
c("GO:0046983","protein dimerization activity",1.180,5.412,0.729,0.648,"identical protein binding"),
c("GO:0043175","RNA polymerase core enzyme binding",0.030,3.009,0.757,0.668,"identical protein binding"),
c("GO:0019903","protein phosphatase binding",0.035,2.773,0.754,0.675,"identical protein binding"),
c("GO:0019902","phosphatase binding",0.044,2.611,0.751,0.686,"identical protein binding"),
c("GO:0016817","hydrolase activity, acting on acid anhydrides",2.073,5.135,0.908,0.046,"hydrolase activity, acting on acid anhydrides"),
c("GO:0016788","hydrolase activity, acting on ester bonds",4.977,4.003,0.899,0.384,"hydrolase activity, acting on acid anhydrides"),
c("GO:0003682","chromatin binding",0.187,6.607,0.960,0.047,"chromatin binding"),
c("GO:0001067","transcription regulatory region nucleic acid binding",0.371,5.363,0.932,0.049,"transcription regulatory region nucleic acid binding"),
c("GO:0030554","adenyl nucleotide binding",14.085,11.360,0.918,0.151,"transcription regulatory region nucleic acid binding"),
c("GO:0061980","regulatory RNA binding",0.010,3.825,0.941,0.180,"transcription regulatory region nucleic acid binding"),
c("GO:0043565","sequence-specific DNA binding",1.642,4.174,0.919,0.275,"transcription regulatory region nucleic acid binding"),
c("GO:0035497","cAMP response element binding",0.006,2.550,0.946,0.278,"transcription regulatory region nucleic acid binding"),
c("GO:0035198","miRNA binding",0.007,3.816,0.942,0.290,"transcription regulatory region nucleic acid binding"),
c("GO:0003727","single-stranded RNA binding",0.051,2.739,0.935,0.323,"transcription regulatory region nucleic acid binding"),
c("GO:0003723","RNA binding",5.649,18.245,0.913,0.382,"transcription regulatory region nucleic acid binding"),
c("GO:0032553","ribonucleotide binding",16.281,11.146,0.917,0.683,"transcription regulatory region nucleic acid binding"),
c("GO:0140097","catalytic activity, acting on DNA",3.292,2.680,0.972,0.050,"catalytic activity, acting on DNA"),
c("GO:0004672","protein kinase activity",3.471,4.827,0.819,0.050,"protein kinase activity"),
c("GO:0106019","phosphatidylinositol-4,5-bisphosphate phosphatase activity",0.015,3.819,0.859,0.278,"protein kinase activity"),
c("GO:0008757","S-adenosylmethionine-dependent methyltransferase activity",1.002,2.653,0.918,0.300,"protein kinase activity"),
c("GO:0016407","acetyltransferase activity",1.205,3.137,0.887,0.307,"protein kinase activity"),
c("GO:0016746","acyltransferase activity",3.068,4.223,0.907,0.349,"protein kinase activity"),
c("GO:0016741","transferase activity, transferring one-carbon groups",3.101,2.830,0.907,0.349,"protein kinase activity"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.428,6.953,0.895,0.409,"protein kinase activity"),
c("GO:0004527","exonuclease activity",0.704,3.027,0.849,0.415,"protein kinase activity"),
c("GO:0061659","ubiquitin-like protein ligase activity",0.225,2.660,0.899,0.428,"protein kinase activity"),
c("GO:0019787","ubiquitin-like protein transferase activity",0.498,2.805,0.892,0.467,"protein kinase activity"),
c("GO:0008238","exopeptidase activity",0.956,2.531,0.867,0.505,"protein kinase activity"),
c("GO:0004177","aminopeptidase activity",0.402,2.676,0.872,0.598,"protein kinase activity"),
c("GO:0016791","phosphatase activity",1.224,2.871,0.825,0.617,"protein kinase activity"),
c("GO:0004713","protein tyrosine kinase activity",0.165,2.523,0.863,0.620,"protein kinase activity"),
c("GO:0034212","peptide N-acetyltransferase activity",0.083,3.274,0.908,0.635,"protein kinase activity"),
c("GO:0004532","exoribonuclease activity",0.136,3.015,0.838,0.652,"protein kinase activity"),
c("GO:0042578","phosphoric ester hydrolase activity",1.589,3.031,0.861,0.667,"protein kinase activity"),
c("GO:0008239","dipeptidyl-peptidase activity",0.047,2.511,0.891,0.699,"protein kinase activity"),
c("GO:0016874","ligase activity",3.484,2.561,0.972,0.050,"ligase activity"),
c("GO:0044877","protein-containing complex binding",0.691,7.033,0.956,0.052,"protein-containing complex binding"),
c("GO:0008270","zinc ion binding",3.633,3.553,0.935,0.066,"zinc ion binding"),
c("GO:0005509","calcium ion binding",1.076,3.511,0.942,0.366,"zinc ion binding"),
c("GO:0046914","transition metal ion binding",5.860,3.571,0.932,0.467,"zinc ion binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_mf_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

