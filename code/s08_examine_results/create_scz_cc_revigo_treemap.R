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
revigo.data <- rbind(c("GO:0043005","neuron projection",0.279,18.779,0.876,0.000,"neuron projection"),
c("GO:0035253","ciliary rootlet",0.003,2.182,0.708,0.596,"neuron projection"),
c("GO:0043198","dendritic shaft",0.004,3.853,0.891,0.603,"neuron projection"),
c("GO:0044327","dendritic spine head",0.001,2.182,0.769,0.614,"neuron projection"),
c("GO:0120111","neuron projection cytoplasm",0.015,7.215,0.765,0.662,"neuron projection"),
c("GO:1902494","catalytic complex",3.457,12.359,0.855,0.000,"catalytic complex"),
c("GO:0035631","CD40 receptor complex",0.002,2.728,0.847,0.266,"catalytic complex"),
c("GO:0000164","protein phosphatase type 1 complex",0.009,2.182,0.806,0.294,"catalytic complex"),
c("GO:0043235","receptor complex",0.111,4.331,0.892,0.368,"catalytic complex"),
c("GO:0044815","DNA packaging complex",0.229,4.924,0.886,0.396,"catalytic complex"),
c("GO:0005681","spliceosomal complex",0.403,3.503,0.685,0.421,"catalytic complex"),
c("GO:0032993","protein-DNA complex",0.414,3.405,0.880,0.422,"catalytic complex"),
c("GO:1904949","ATPase complex",1.228,2.134,0.869,0.481,"catalytic complex"),
c("GO:1990351","transporter complex",1.273,2.951,0.868,0.483,"catalytic complex"),
c("GO:1990234","transferase complex",1.442,6.200,0.850,0.491,"catalytic complex"),
c("GO:0098797","plasma membrane protein complex",1.538,5.118,0.775,0.496,"catalytic complex"),
c("GO:0042613","MHC class II protein complex",0.012,3.389,0.829,0.496,"catalytic complex"),
c("GO:0042611","MHC protein complex",0.015,5.589,0.833,0.504,"catalytic complex"),
c("GO:1990904","ribonucleoprotein complex",2.023,4.337,0.862,0.514,"catalytic complex"),
c("GO:0030123","AP-3 adaptor complex",0.022,2.550,0.781,0.520,"catalytic complex"),
c("GO:0098796","membrane protein complex",3.375,9.999,0.844,0.553,"catalytic complex"),
c("GO:0071011","precatalytic spliceosome",0.036,3.069,0.728,0.569,"catalytic complex"),
c("GO:0034703","cation channel complex",0.171,3.144,0.844,0.621,"catalytic complex"),
c("GO:0000151","ubiquitin ligase complex",0.260,3.572,0.869,0.637,"catalytic complex"),
c("GO:1905369","endopeptidase complex",0.468,2.586,0.863,0.676,"catalytic complex"),
c("GO:1905368","peptidase complex",0.561,4.131,0.861,0.690,"catalytic complex"),
c("GO:0030496","midbody",0.055,2.804,1.000,0.000,"midbody"),
c("GO:0030427","site of polarized growth",0.079,4.234,1.000,0.000,"site of polarized growth"),
c("GO:0044297","cell body",0.083,11.793,1.000,0.000,"cell body"),
c("GO:0031252","cell leading edge",0.121,6.085,1.000,0.000,"cell leading edge"),
c("GO:0036477","somatodendritic compartment",0.153,17.262,1.000,0.000,"somatodendritic compartment"),
c("GO:0032153","cell division site",0.177,3.633,1.000,0.000,"cell division site"),
c("GO:0016604","nuclear body",0.231,15.322,0.766,0.000,"nuclear body"),
c("GO:0015630","microtubule cytoskeleton",1.054,12.418,0.719,0.197,"nuclear body"),
c("GO:0005739","mitochondrion",2.679,9.175,0.751,0.280,"nuclear body"),
c("GO:0042579","microbody",0.279,2.590,0.794,0.285,"nuclear body"),
c("GO:0005774","vacuolar membrane",0.368,5.609,0.728,0.294,"nuclear body"),
c("GO:0005773","vacuole",0.631,5.386,0.780,0.312,"nuclear body"),
c("GO:0005794","Golgi apparatus",1.491,8.072,0.702,0.346,"nuclear body"),
c("GO:0031984","organelle subcompartment",2.183,2.194,0.767,0.359,"nuclear body"),
c("GO:0098576","lumenal side of membrane",0.006,5.100,0.819,0.392,"nuclear body"),
c("GO:0035770","ribonucleoprotein granule",0.123,2.672,0.775,0.397,"nuclear body"),
c("GO:0000785","chromatin",0.550,12.217,0.730,0.458,"nuclear body"),
c("GO:0001650","fibrillar center",0.033,2.885,0.754,0.521,"nuclear body"),
c("GO:0034399","nuclear periphery",0.048,2.923,0.792,0.535,"nuclear body"),
c("GO:0031300","intrinsic component of organelle membrane",0.405,2.300,0.759,0.542,"nuclear body"),
c("GO:0035327","transcriptionally active chromatin",0.009,2.624,0.793,0.554,"nuclear body"),
c("GO:0031227","intrinsic component of endoplasmic reticulum membrane",0.188,3.925,0.711,0.582,"nuclear body"),
c("GO:0005759","mitochondrial matrix",0.366,2.951,0.749,0.588,"nuclear body"),
c("GO:0030133","transport vesicle",0.216,6.387,0.671,0.590,"nuclear body"),
c("GO:0000792","heterochromatin",0.029,3.513,0.779,0.600,"nuclear body"),
c("GO:0030666","endocytic vesicle membrane",0.030,3.085,0.701,0.610,"nuclear body"),
c("GO:0005635","nuclear envelope",0.331,3.317,0.721,0.615,"nuclear body"),
c("GO:0070382","exocytic vesicle",0.066,3.697,0.691,0.645,"nuclear body"),
c("GO:0030139","endocytic vesicle",0.067,4.059,0.722,0.646,"nuclear body"),
c("GO:0000932","P-body",0.072,2.885,0.775,0.646,"nuclear body"),
c("GO:0072686","mitotic spindle",0.070,4.924,0.760,0.652,"nuclear body"),
c("GO:0016607","nuclear speck",0.127,10.698,0.772,0.654,"nuclear body"),
c("GO:0030134","COPII-coated ER to Golgi transport vesicle",0.091,2.787,0.708,0.661,"nuclear body"),
c("GO:0098800","inner mitochondrial membrane protein complex",0.236,2.578,0.632,0.668,"nuclear body"),
c("GO:0030660","Golgi-associated vesicle membrane",0.104,2.365,0.680,0.670,"nuclear body"),
c("GO:0005798","Golgi-associated vesicle",0.115,2.747,0.712,0.673,"nuclear body"),
c("GO:0005730","nucleolus",0.854,5.854,0.690,0.677,"nuclear body"),
c("GO:0000781","chromosome, telomeric region",0.133,5.183,0.747,0.678,"nuclear body"),
c("GO:0099503","secretory vesicle",0.131,3.359,0.709,0.680,"nuclear body"),
c("GO:0005770","late endosome",0.154,3.014,0.676,0.689,"nuclear body"),
c("GO:0009986","cell surface",0.232,2.816,1.000,0.000,"cell surface"),
c("GO:0045202","synapse",0.487,17.186,0.825,0.000,"synapse"),
c("GO:0060076","excitatory synapse",0.007,3.009,0.851,0.633,"synapse"),
c("GO:0098839","postsynaptic density membrane",0.011,2.953,0.675,0.650,"synapse"),
c("GO:0098982","GABA-ergic synapse",0.011,2.392,0.848,0.652,"synapse"),
c("GO:0048786","presynaptic active zone",0.019,2.859,0.840,0.677,"synapse"),
c("GO:0098590","plasma membrane region",0.523,10.830,0.919,0.000,"plasma membrane region"),
c("GO:0030315","T-tubule",0.008,3.854,0.938,0.213,"plasma membrane region"),
c("GO:0032154","cleavage furrow",0.021,2.843,0.886,0.228,"plasma membrane region"),
c("GO:0042383","sarcolemma",0.040,2.773,0.932,0.240,"plasma membrane region"),
c("GO:0005938","cell cortex",0.251,4.348,0.811,0.260,"plasma membrane region"),
c("GO:0031226","intrinsic component of plasma membrane",1.546,5.877,0.911,0.333,"plasma membrane region"),
c("GO:0016323","basolateral plasma membrane",0.047,2.817,0.880,0.647,"plasma membrane region"),
c("GO:0031975","envelope",2.867,6.024,1.000,0.000,"envelope"),
c("GO:0016020","membrane",63.598,8.907,1.000,0.000,"membrane"),
c("GO:0098552","side of membrane",0.189,4.139,0.986,0.029,"side of membrane"),
c("GO:0099568","cytoplasmic region",0.098,8.222,0.877,0.087,"cytoplasmic region"),
c("GO:0048471","perinuclear region of cytoplasm",0.158,8.052,0.873,0.156,"cytoplasmic region"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_cc_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

