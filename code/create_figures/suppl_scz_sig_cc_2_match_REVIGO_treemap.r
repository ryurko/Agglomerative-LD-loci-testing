

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

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0043005","neuron projection",0.190,16.4191,0.755,0.000,"neuron projection"),
c("GO:0097542","ciliary tip",0.003,2.0875,0.688,0.538,"neuron projection"),
c("GO:0071565","nBAF complex",0.003,2.0905,0.562,0.692,"neuron projection"),
c("GO:0032838","cell projection cytoplasm",0.014,5.0729,0.700,0.595,"neuron projection"),
c("GO:0031253","cell projection membrane",0.053,2.9442,0.724,0.653,"neuron projection"),
c("GO:0097447","dendritic tree",0.001,11.1056,0.842,0.627,"neuron projection"),
c("GO:0045202","synapse",0.299,15.5395,0.986,0.000,"synapse"),
c("GO:0070161","anchoring junction",0.151,5.8429,0.974,0.000,"anchoring junction"),
c("GO:0098805","whole membrane",0.888,7.0967,0.970,0.000,"whole membrane"),
c("GO:1902494","catalytic complex",3.734,10.2755,0.824,0.000,"catalytic complex"),
c("GO:0034703","cation channel complex",0.153,2.4746,0.841,0.346,"catalytic complex"),
c("GO:1905348","endonuclease complex",0.129,2.2041,0.839,0.572,"catalytic complex"),
c("GO:1990904","ribonucleoprotein complex",5.291,4.9949,0.818,0.543,"catalytic complex"),
c("GO:0031332","RNAi effector complex",0.017,2.5781,0.783,0.282,"catalytic complex"),
c("GO:0043235","receptor complex",0.115,2.9482,0.864,0.336,"catalytic complex"),
c("GO:1905368","peptidase complex",0.452,2.5201,0.825,0.646,"catalytic complex"),
c("GO:0030123","AP-3 adaptor complex",0.019,2.5781,0.696,0.285,"catalytic complex"),
c("GO:1990351","transporter complex",0.885,2.2591,0.843,0.422,"catalytic complex"),
c("GO:0044815","DNA packaging complex",0.216,5.1042,0.858,0.359,"catalytic complex"),
c("GO:0005839","proteasome core complex",0.183,2.1409,0.720,0.591,"catalytic complex"),
c("GO:1904949","ATPase complex",0.680,2.3964,0.819,0.674,"catalytic complex"),
c("GO:0000164","protein phosphatase type 1 complex",0.009,2.0905,0.713,0.460,"catalytic complex"),
c("GO:0000151","ubiquitin ligase complex",0.232,3.5642,0.731,0.604,"catalytic complex"),
c("GO:0098796","membrane protein complex",2.473,7.9878,0.793,0.484,"catalytic complex"),
c("GO:1990234","transferase complex",1.223,4.1168,0.811,0.439,"catalytic complex"),
c("GO:0032993","protein-DNA complex",0.418,4.8444,0.851,0.386,"catalytic complex"),
c("GO:0098552","side of membrane",0.214,3.3689,0.969,0.034,"side of membrane"),
c("GO:0032154","cleavage furrow",0.012,2.1911,0.924,0.040,"cleavage furrow"),
c("GO:0030496","midbody",0.040,2.1103,0.919,0.044,"midbody"),
c("GO:0031252","cell leading edge",0.086,6.0001,0.915,0.046,"cell leading edge"),
c("GO:0044297","cell body",0.087,8.1296,0.915,0.046,"cell body"),
c("GO:0030427","site of polarized growth",0.091,3.1238,0.915,0.046,"site of polarized growth"),
c("GO:0098794","postsynapse",0.133,10.4981,0.817,0.048,"postsynapse"),
c("GO:0098831","presynaptic active zone cytoplasmic component",0.001,2.1409,0.646,0.663,"postsynapse"),
c("GO:0098978","glutamatergic synapse",0.000,6.2305,0.920,0.510,"postsynapse"),
c("GO:0060076","excitatory synapse",0.004,2.6005,0.901,0.696,"postsynapse"),
c("GO:0016604","nuclear body",0.189,14.6199,0.625,0.049,"nuclear body"),
c("GO:0001650","fibrillar center",0.035,3.0621,0.647,0.519,"nuclear body"),
c("GO:0031300","intrinsic component of organelle membrane",0.221,3.0998,0.656,0.641,"nuclear body"),
c("GO:0035770","ribonucleoprotein granule",0.131,3.0206,0.640,0.411,"nuclear body"),
c("GO:0034399","nuclear periphery",0.066,3.3035,0.652,0.543,"nuclear body"),
c("GO:0072686","mitotic spindle",0.078,4.4404,0.639,0.693,"nuclear body"),
c("GO:0098576","lumenal side of membrane",0.000,6.7111,0.770,0.213,"nuclear body"),
c("GO:0005739","mitochondrion",2.156,7.2802,0.638,0.426,"nuclear body"),
c("GO:0005730","nucleolus",0.664,5.9249,0.580,0.654,"nuclear body"),
c("GO:0071011","precatalytic spliceosome",0.027,2.5432,0.608,0.463,"nuclear body"),
c("GO:0030133","transport vesicle",0.179,6.8698,0.576,0.633,"nuclear body"),
c("GO:0031227","intrinsic component of endoplasmic reticulum membrane",0.142,4.6483,0.591,0.620,"nuclear body"),
c("GO:0005635","nuclear envelope",0.283,3.7189,0.562,0.661,"nuclear body"),
c("GO:0000785","chromatin",0.473,11.7446,0.538,0.463,"nuclear body"),
c("GO:0005681","spliceosomal complex",0.250,3.2590,0.559,0.548,"nuclear body"),
c("GO:0000792","heterochromatin",0.045,2.3138,0.601,0.681,"nuclear body"),
c("GO:0042579","microbody",0.215,2.7010,0.694,0.299,"nuclear body"),
c("GO:0015630","microtubule cytoskeleton",0.900,12.7380,0.641,0.235,"nuclear body"),
c("GO:0031901","early endosome membrane",0.016,2.4384,0.571,0.648,"nuclear body"),
c("GO:0012506","vesicle membrane",0.223,4.4538,0.609,0.682,"nuclear body"),
c("GO:0030904","retromer complex",0.032,2.6056,0.691,0.548,"nuclear body"),
c("GO:0030666","endocytic vesicle membrane",0.019,2.6809,0.587,0.656,"nuclear body"),
c("GO:0099568","cytoplasmic region",0.265,8.4271,0.758,0.109,"nuclear body"),
c("GO:0005774","vacuolar membrane",0.290,4.7969,0.585,0.390,"nuclear body"),
c("GO:0005773","vacuole",0.455,4.3986,0.678,0.360,"nuclear body"),
c("GO:0031984","organelle subcompartment",0.269,2.2046,0.658,0.361,"nuclear body"),
c("GO:0005759","mitochondrial matrix",0.336,2.5140,0.603,0.585,"nuclear body"),
c("GO:0000932","P-body",0.055,3.6043,0.621,0.382,"nuclear body"),
c("GO:0070382","exocytic vesicle",0.048,3.6566,0.602,0.700,"nuclear body"),
c("GO:0048471","perinuclear region of cytoplasm",0.135,7.2541,0.768,0.288,"nuclear body"),
c("GO:0005794","Golgi apparatus",0.969,7.6747,0.587,0.342,"nuclear body"),
c("GO:0098590","plasma membrane region",0.239,8.0825,0.801,0.050,"plasma membrane region"),
c("GO:0042611","MHC protein complex",0.012,6.2284,0.737,0.482,"plasma membrane region"),
c("GO:0031256","leading edge membrane",0.026,3.0733,0.794,0.509,"plasma membrane region"),
c("GO:0042613","MHC class II protein complex",0.008,3.4341,0.737,0.501,"plasma membrane region"),
c("GO:0030315","T-tubule",0.008,3.3365,0.836,0.469,"plasma membrane region"),
c("GO:0032153","cell division site",0.238,2.6999,0.910,0.051,"cell division site"),
c("GO:0031975","envelope",2.324,5.0715,0.895,0.063,"envelope"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_cc_2_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
