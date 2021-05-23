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
revigo.data <- rbind(c("GO:0022008","neurogenesis",0.356,12.137,0.710,0.000,"neurogenesis"),
c("GO:0001964","startle response",0.005,3.568,0.855,0.456,"neurogenesis"),
c("GO:0060322","head development",0.148,5.311,0.797,0.628,"neurogenesis"),
c("GO:0021872","forebrain generation of neurons",0.011,5.303,0.736,0.700,"neurogenesis"),
c("GO:0048002","antigen processing and presentation of peptide antigen",0.011,7.159,0.942,0.000,"antigen processing and presentation of peptide antigen"),
c("GO:0019882","antigen processing and presentation",0.027,5.884,0.961,0.542,"antigen processing and presentation of peptide antigen"),
c("GO:0099177","regulation of trans-synaptic signaling",0.079,9.633,0.849,0.000,"regulation of trans-synaptic signaling"),
c("GO:0098693","regulation of synaptic vesicle cycle",0.002,4.111,0.906,0.143,"regulation of trans-synaptic signaling"),
c("GO:0050803","regulation of synapse structure or activity",0.044,4.807,0.882,0.153,"regulation of trans-synaptic signaling"),
c("GO:0043087","regulation of GTPase activity",0.157,3.876,0.891,0.168,"regulation of trans-synaptic signaling"),
c("GO:0051960","regulation of nervous system development",0.075,6.924,0.819,0.171,"regulation of trans-synaptic signaling"),
c("GO:0010975","regulation of neuron projection development",0.066,6.029,0.815,0.180,"regulation of trans-synaptic signaling"),
c("GO:1902680","positive regulation of RNA biosynthetic process",0.701,7.254,0.775,0.215,"regulation of trans-synaptic signaling"),
c("GO:0051726","regulation of cell cycle",0.598,3.666,0.873,0.258,"regulation of trans-synaptic signaling"),
c("GO:0042176","regulation of protein catabolic process",0.160,3.364,0.881,0.320,"regulation of trans-synaptic signaling"),
c("GO:0009894","regulation of catabolic process",0.358,4.378,0.875,0.332,"regulation of trans-synaptic signaling"),
c("GO:0034248","regulation of cellular amide metabolic process",1.527,4.090,0.854,0.399,"regulation of trans-synaptic signaling"),
c("GO:2000310","regulation of NMDA receptor activity",0.002,3.434,0.865,0.508,"regulation of trans-synaptic signaling"),
c("GO:0010769","regulation of cell morphogenesis involved in differentiation",0.015,4.653,0.843,0.537,"regulation of trans-synaptic signaling"),
c("GO:0045664","regulation of neuron differentiation",0.039,6.043,0.820,0.571,"regulation of trans-synaptic signaling"),
c("GO:0048167","regulation of synaptic plasticity",0.028,6.128,0.829,0.598,"regulation of trans-synaptic signaling"),
c("GO:0031344","regulation of cell projection organization",0.124,4.931,0.833,0.650,"regulation of trans-synaptic signaling"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.144,11.755,0.960,0.008,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0099536","synaptic signaling",0.150,8.459,0.939,0.008,"synaptic signaling"),
c("GO:0007267","cell-cell signaling",0.336,6.321,0.963,0.400,"synaptic signaling"),
c("GO:0030030","cell projection organization",0.664,10.391,0.854,0.009,"cell projection organization"),
c("GO:0000422","autophagy of mitochondrion",0.023,4.756,0.826,0.359,"cell projection organization"),
c("GO:0050808","synapse organization",0.059,8.780,0.849,0.385,"cell projection organization"),
c("GO:0009057","macromolecule catabolic process",2.208,3.691,0.919,0.438,"cell projection organization"),
c("GO:0030031","cell projection assembly",0.328,3.757,0.825,0.447,"cell projection organization"),
c("GO:0051276","chromosome organization",2.083,4.760,0.827,0.540,"cell projection organization"),
c("GO:0043632","modification-dependent macromolecule catabolic process",0.715,3.915,0.900,0.629,"cell projection organization"),
c("GO:0070925","organelle assembly",0.650,3.787,0.832,0.677,"cell projection organization"),
c("GO:0006325","chromatin organization",0.761,3.725,0.838,0.690,"cell projection organization"),
c("GO:0007049","cell cycle",1.745,5.565,0.989,0.011,"cell cycle"),
c("GO:0007017","microtubule-based process",0.707,4.221,0.990,0.011,"microtubule-based process"),
c("GO:0022402","cell cycle process",0.856,4.003,0.990,0.011,"cell cycle process"),
c("GO:0070647","protein modification by small protein conjugation or removal",1.191,4.780,0.920,0.012,"protein modification by small protein conjugation or removal"),
c("GO:0006367","transcription initiation from RNA polymerase II promoter",0.095,3.531,0.943,0.179,"protein modification by small protein conjugation or removal"),
c("GO:0006281","DNA repair",2.380,3.378,0.898,0.319,"protein modification by small protein conjugation or removal"),
c("GO:0018022","peptidyl-lysine methylation",0.141,3.729,0.929,0.401,"protein modification by small protein conjugation or removal"),
c("GO:0032446","protein modification by small protein conjugation",0.925,3.770,0.915,0.483,"protein modification by small protein conjugation or removal"),
c("GO:0098781","ncRNA transcription",0.023,3.446,0.948,0.500,"protein modification by small protein conjugation or removal"),
c("GO:0006259","DNA metabolic process",5.371,4.295,0.925,0.512,"protein modification by small protein conjugation or removal"),
c("GO:0018205","peptidyl-lysine modification",0.449,3.215,0.922,0.584,"protein modification by small protein conjugation or removal"),
c("GO:0006513","protein monoubiquitination",0.034,3.134,0.933,0.691,"protein modification by small protein conjugation or removal"),
c("GO:0046907","intracellular transport",1.822,5.508,0.904,0.012,"intracellular transport"),
c("GO:0048193","Golgi vesicle transport",0.366,4.287,0.948,0.281,"intracellular transport"),
c("GO:0051648","vesicle localization",0.055,4.806,0.920,0.635,"intracellular transport"),
c("GO:0030705","cytoskeleton-dependent intracellular transport",0.071,3.670,0.922,0.684,"intracellular transport"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_bp_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

