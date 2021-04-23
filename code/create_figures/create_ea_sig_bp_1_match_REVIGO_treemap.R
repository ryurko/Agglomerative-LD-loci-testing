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
revigo.data <- rbind(c("GO:0048002","antigen processing and presentation of peptide antigen",0.010,6.097,0.920,0.000,"antigen processing and presentation of peptide antigen"),
c("GO:0019882","antigen processing and presentation",0.025,4.492,0.943,0.579,"antigen processing and presentation of peptide antigen"),
c("GO:0048666","neuron development",0.229,13.960,0.638,0.000,"neuron development"),
c("GO:0007613","memory",0.018,3.460,0.827,0.492,"neuron development"),
c("GO:0060322","head development",0.149,4.433,0.753,0.608,"neuron development"),
c("GO:0099177","regulation of trans-synaptic signaling",0.079,8.867,0.852,0.000,"regulation of trans-synaptic signaling"),
c("GO:0098693","regulation of synaptic vesicle cycle",0.002,3.460,0.886,0.143,"regulation of trans-synaptic signaling"),
c("GO:0050803","regulation of synapse structure or activity",0.043,4.925,0.877,0.153,"regulation of trans-synaptic signaling"),
c("GO:0051960","regulation of nervous system development",0.076,6.972,0.813,0.171,"regulation of trans-synaptic signaling"),
c("GO:0043087","regulation of GTPase activity",0.323,4.362,0.875,0.177,"regulation of trans-synaptic signaling"),
c("GO:0010975","regulation of neuron projection development",0.067,7.049,0.780,0.181,"regulation of trans-synaptic signaling"),
c("GO:1902680","positive regulation of RNA biosynthetic process",0.701,4.347,0.800,0.216,"regulation of trans-synaptic signaling"),
c("GO:0009894","regulation of catabolic process",0.277,3.709,0.870,0.324,"regulation of trans-synaptic signaling"),
c("GO:0034248","regulation of cellular amide metabolic process",1.463,3.321,0.847,0.398,"regulation of trans-synaptic signaling"),
c("GO:0060078","regulation of postsynaptic membrane potential",0.045,3.737,0.877,0.437,"regulation of trans-synaptic signaling"),
c("GO:2000310","regulation of NMDA receptor activity",0.002,3.738,0.854,0.514,"regulation of trans-synaptic signaling"),
c("GO:0010769","regulation of cell morphogenesis involved in differentiation",0.015,4.694,0.832,0.535,"regulation of trans-synaptic signaling"),
c("GO:0045664","regulation of neuron differentiation",0.039,6.532,0.812,0.570,"regulation of trans-synaptic signaling"),
c("GO:0060341","regulation of cellular localization",0.162,3.424,0.855,0.570,"regulation of trans-synaptic signaling"),
c("GO:0048167","regulation of synaptic plasticity",0.029,5.540,0.841,0.604,"regulation of trans-synaptic signaling"),
c("GO:0031344","regulation of cell projection organization",0.127,5.932,0.804,0.654,"regulation of trans-synaptic signaling"),
c("GO:0018022","peptidyl-lysine methylation",0.125,3.772,0.994,0.005,"peptidyl-lysine methylation"),
c("GO:0007156","homophilic cell adhesion via plasma membrane adhesion molecules",0.135,11.329,0.952,0.005,"homophilic cell adhesion via plasma membrane adhesion molecules"),
c("GO:0099536","synaptic signaling",0.146,7.672,0.954,0.005,"synaptic signaling"),
c("GO:0048013","ephrin receptor signaling pathway",0.021,3.518,0.853,0.329,"synaptic signaling"),
c("GO:0007267","cell-cell signaling",0.320,4.459,0.966,0.402,"synaptic signaling"),
c("GO:0030030","cell projection organization",0.670,10.197,0.838,0.006,"cell projection organization"),
c("GO:0000422","autophagy of mitochondrion",0.023,3.953,0.857,0.362,"cell projection organization"),
c("GO:0050808","synapse organization",0.058,9.723,0.827,0.389,"cell projection organization"),
c("GO:0030031","cell projection assembly",0.331,3.871,0.798,0.452,"cell projection organization"),
c("GO:0051276","chromosome organization",2.075,3.558,0.816,0.545,"cell projection organization"),
c("GO:0070925","organelle assembly",0.627,3.292,0.817,0.678,"cell projection organization"),
c("GO:0046907","intracellular transport",1.785,4.967,0.914,0.008,"intracellular transport"),
c("GO:1901998","toxin transport",0.016,3.353,0.972,0.219,"intracellular transport"),
c("GO:0048193","Golgi vesicle transport",0.356,3.867,0.961,0.290,"intracellular transport"),
c("GO:0051668","localization within membrane",0.013,3.734,0.938,0.565,"intracellular transport"),
c("GO:0051648","vesicle localization",0.055,4.732,0.932,0.636,"intracellular transport"),
c("GO:0007017","microtubule-based process",0.697,3.779,0.994,0.008,"microtubule-based process"),
c("GO:0007049","cell cycle",1.774,3.424,0.993,0.008,"cell cycle"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_ea_sig_bp_1_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

