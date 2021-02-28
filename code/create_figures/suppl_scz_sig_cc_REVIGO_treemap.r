

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
revigo.data <- rbind(c("GO:0043005","neuron projection",0.190,32.3693,0.731,0.000,"neuron projection"),
c("GO:0098858","actin-based cell projection",0.041,2.7092,0.809,0.641,"neuron projection"),
c("GO:0032838","cell projection cytoplasm",0.014,10.9241,0.691,0.595,"neuron projection"),
c("GO:0031253","cell projection membrane",0.053,5.4059,0.688,0.653,"neuron projection"),
c("GO:0097447","dendritic tree",0.001,19.4587,0.852,0.627,"neuron projection"),
c("GO:0097014","ciliary plasm",0.044,3.6262,0.626,0.644,"neuron projection"),
c("GO:0031143","pseudopodium",0.005,2.1227,0.830,0.554,"neuron projection"),
c("GO:0045202","synapse",0.299,32.8397,0.991,0.000,"synapse"),
c("GO:0070161","anchoring junction",0.151,9.9195,0.961,0.000,"anchoring junction"),
c("GO:0098796","membrane protein complex",2.473,18.1695,0.814,0.000,"membrane protein complex"),
c("GO:0031332","RNAi effector complex",0.017,3.2970,0.809,0.273,"membrane protein complex"),
c("GO:0043235","receptor complex",0.115,4.7675,0.883,0.323,"membrane protein complex"),
c("GO:1990111","spermatoproteasome complex",0.001,3.1723,0.802,0.343,"membrane protein complex"),
c("GO:0098878","neurotransmitter receptor complex",0.012,3.7957,0.897,0.266,"membrane protein complex"),
c("GO:0044815","DNA packaging complex",0.216,4.1037,0.878,0.344,"membrane protein complex"),
c("GO:0005839","proteasome core complex",0.183,2.3816,0.739,0.628,"membrane protein complex"),
c("GO:1904949","ATPase complex",0.680,2.1278,0.839,0.674,"membrane protein complex"),
c("GO:0031461","cullin-RING ubiquitin ligase complex",0.159,2.4695,0.757,0.695,"membrane protein complex"),
c("GO:0032993","protein-DNA complex",0.418,4.8781,0.872,0.368,"membrane protein complex"),
c("GO:0034703","cation channel complex",0.153,6.7014,0.851,0.332,"membrane protein complex"),
c("GO:0034708","methyltransferase complex",0.080,2.2477,0.767,0.659,"membrane protein complex"),
c("GO:1905348","endonuclease complex",0.129,3.6721,0.856,0.572,"membrane protein complex"),
c("GO:1905368","peptidase complex",0.452,2.1602,0.843,0.646,"membrane protein complex"),
c("GO:0030123","AP-3 adaptor complex",0.019,3.5968,0.725,0.419,"membrane protein complex"),
c("GO:1990351","transporter complex",0.885,7.0350,0.865,0.401,"membrane protein complex"),
c("GO:0000164","protein phosphatase type 1 complex",0.009,2.5956,0.737,0.460,"membrane protein complex"),
c("GO:0000151","ubiquitin ligase complex",0.232,4.4610,0.751,0.604,"membrane protein complex"),
c("GO:0070578","RISC-loading complex",0.002,1.9294,0.808,0.420,"membrane protein complex"),
c("GO:1990234","transferase complex",1.223,6.4510,0.832,0.417,"membrane protein complex"),
c("GO:1902494","catalytic complex",3.734,13.7097,0.849,0.484,"membrane protein complex"),
c("GO:0099080","supramolecular complex",0.540,5.1902,0.991,0.000,"supramolecular complex"),
c("GO:0032154","cleavage furrow",0.012,2.2305,0.923,0.040,"cleavage furrow"),
c("GO:0098862","cluster of actin-based cell projections",0.034,3.4900,0.926,0.043,"cluster of actin-based cell projections"),
c("GO:0031225","anchored component of membrane",0.078,1.9414,0.974,0.044,"anchored component of membrane"),
c("GO:0030496","midbody",0.040,2.6696,0.925,0.044,"midbody"),
c("GO:0098805","whole membrane",0.888,16.1855,0.974,0.044,"whole membrane"),
c("GO:0043209","myelin sheath",0.049,2.0127,0.924,0.044,"myelin sheath"),
c("GO:0031252","cell leading edge",0.086,9.5711,0.922,0.046,"cell leading edge"),
c("GO:0044297","cell body",0.087,13.8582,0.922,0.046,"cell body"),
c("GO:0030427","site of polarized growth",0.091,4.5892,0.921,0.046,"site of polarized growth"),
c("GO:0098794","postsynapse",0.133,21.3120,0.805,0.048,"postsynapse"),
c("GO:0098831","presynaptic active zone cytoplasmic component",0.001,3.0402,0.659,0.663,"postsynapse"),
c("GO:0098890","extrinsic component of postsynaptic membrane",0.000,1.9294,0.798,0.568,"postsynapse"),
c("GO:0098978","glutamatergic synapse",0.000,12.6773,0.912,0.510,"postsynapse"),
c("GO:0098982","GABA-ergic synapse",0.000,3.5695,0.901,0.602,"postsynapse"),
c("GO:0098936","intrinsic component of postsynaptic membrane",0.001,3.9063,0.775,0.658,"postsynapse"),
c("GO:0060076","excitatory synapse",0.004,3.8816,0.890,0.696,"postsynapse"),
c("GO:0097470","ribbon synapse",0.001,1.9294,0.895,0.651,"postsynapse"),
c("GO:0044327","dendritic spine head",0.000,1.9294,0.736,0.624,"postsynapse"),
c("GO:0098552","side of membrane",0.214,5.2943,0.973,0.048,"side of membrane"),
c("GO:0019898","extrinsic component of membrane",0.249,1.9569,0.973,0.049,"extrinsic component of membrane"),
c("GO:0015630","microtubule cytoskeleton",0.900,25.5032,0.661,0.056,"microtubule cytoskeleton"),
c("GO:0031298","replication fork protection complex",0.020,2.1453,0.609,0.641,"microtubule cytoskeleton"),
c("GO:0031300","intrinsic component of organelle membrane",0.221,5.2830,0.674,0.641,"microtubule cytoskeleton"),
c("GO:0034399","nuclear periphery",0.066,3.4969,0.678,0.543,"microtubule cytoskeleton"),
c("GO:0005739","mitochondrion",2.156,15.4808,0.668,0.426,"microtubule cytoskeleton"),
c("GO:0005730","nucleolus",0.664,2.3999,0.613,0.654,"microtubule cytoskeleton"),
c("GO:0043240","Fanconi anaemia nuclear complex",0.009,2.0104,0.653,0.600,"microtubule cytoskeleton"),
c("GO:0005720","nuclear heterochromatin",0.024,2.0171,0.614,0.650,"microtubule cytoskeleton"),
c("GO:0031227","intrinsic component of endoplasmic reticulum membrane",0.142,7.1573,0.622,0.620,"microtubule cytoskeleton"),
c("GO:0005635","nuclear envelope",0.283,6.6509,0.595,0.661,"microtubule cytoskeleton"),
c("GO:0000785","chromatin",0.473,10.3255,0.574,0.463,"microtubule cytoskeleton"),
c("GO:0005681","spliceosomal complex",0.250,2.7919,0.598,0.548,"microtubule cytoskeleton"),
c("GO:0000791","euchromatin",0.011,2.1068,0.658,0.615,"microtubule cytoskeleton"),
c("GO:0000792","heterochromatin",0.045,1.9908,0.631,0.681,"microtubule cytoskeleton"),
c("GO:0042579","microbody",0.215,3.8698,0.719,0.335,"microtubule cytoskeleton"),
c("GO:0016363","nuclear matrix",0.019,2.1177,0.700,0.498,"microtubule cytoskeleton"),
c("GO:0031901","early endosome membrane",0.016,3.6554,0.596,0.648,"microtubule cytoskeleton"),
c("GO:0012506","vesicle membrane",0.223,10.2921,0.627,0.682,"microtubule cytoskeleton"),
c("GO:0044545","NSL complex",0.001,1.9479,0.680,0.509,"microtubule cytoskeleton"),
c("GO:0030904","retromer complex",0.032,3.6940,0.724,0.548,"microtubule cytoskeleton"),
c("GO:0099568","cytoplasmic region",0.265,17.2605,0.780,0.342,"microtubule cytoskeleton"),
c("GO:0031965","nuclear membrane",0.100,3.6136,0.599,0.676,"microtubule cytoskeleton"),
c("GO:0005774","vacuolar membrane",0.290,8.0493,0.608,0.390,"microtubule cytoskeleton"),
c("GO:0005773","vacuole",0.455,7.8231,0.704,0.360,"microtubule cytoskeleton"),
c("GO:0031984","organelle subcompartment",0.269,5.4083,0.679,0.361,"microtubule cytoskeleton"),
c("GO:0005759","mitochondrial matrix",0.336,4.6576,0.628,0.585,"microtubule cytoskeleton"),
c("GO:0035861","site of double-strand break",0.032,2.0937,0.686,0.664,"microtubule cytoskeleton"),
c("GO:0000932","P-body",0.055,6.7679,0.660,0.382,"microtubule cytoskeleton"),
c("GO:0070382","exocytic vesicle",0.048,6.6915,0.626,0.700,"microtubule cytoskeleton"),
c("GO:0048471","perinuclear region of cytoplasm",0.135,12.8636,0.790,0.321,"microtubule cytoskeleton"),
c("GO:0005794","Golgi apparatus",0.969,19.2214,0.624,0.314,"microtubule cytoskeleton"),
c("GO:0001650","fibrillar center",0.035,2.2710,0.678,0.519,"microtubule cytoskeleton"),
c("GO:0016528","sarcoplasm",0.019,1.9908,0.862,0.100,"microtubule cytoskeleton"),
c("GO:0035770","ribonucleoprotein granule",0.131,5.9842,0.675,0.411,"microtubule cytoskeleton"),
c("GO:0072686","mitotic spindle",0.078,4.1871,0.656,0.693,"microtubule cytoskeleton"),
c("GO:0098576","lumenal side of membrane",0.000,8.5343,0.773,0.134,"microtubule cytoskeleton"),
c("GO:0098562","cytoplasmic side of membrane",0.147,2.9212,0.960,0.604,"microtubule cytoskeleton"),
c("GO:0071011","precatalytic spliceosome",0.027,2.6432,0.649,0.463,"microtubule cytoskeleton"),
c("GO:0030133","transport vesicle",0.179,11.3408,0.606,0.633,"microtubule cytoskeleton"),
c("GO:0033106","cis-Golgi network membrane",0.001,1.9479,0.730,0.507,"microtubule cytoskeleton"),
c("GO:0016604","nuclear body",0.189,19.9740,0.651,0.235,"microtubule cytoskeleton"),
c("GO:0042824","MHC class I peptide loading complex",0.001,3.1723,0.656,0.439,"microtubule cytoskeleton"),
c("GO:0035097","histone methyltransferase complex",0.030,2.2680,0.618,0.650,"microtubule cytoskeleton"),
c("GO:0030666","endocytic vesicle membrane",0.019,4.5127,0.603,0.656,"microtubule cytoskeleton"),
c("GO:0005881","cytoplasmic microtubule",0.023,2.3866,0.655,0.628,"microtubule cytoskeleton"),
c("GO:0005884","actin filament",0.037,2.3805,0.672,0.655,"microtubule cytoskeleton"),
c("GO:0043292","contractile fiber",0.067,1.9731,0.695,0.683,"microtubule cytoskeleton"),
c("GO:0032153","cell division site",0.238,2.6006,0.917,0.057,"cell division site"),
c("GO:0098590","plasma membrane region",0.239,19.5164,0.806,0.057,"plasma membrane region"),
c("GO:0016011","dystroglycan complex",0.009,2.4068,0.753,0.505,"plasma membrane region"),
c("GO:0042611","MHC protein complex",0.012,10.2763,0.750,0.482,"plasma membrane region"),
c("GO:0031256","leading edge membrane",0.026,5.1761,0.809,0.509,"plasma membrane region"),
c("GO:0019897","extrinsic component of plasma membrane",0.104,3.0675,0.813,0.564,"plasma membrane region"),
c("GO:0042613","MHC class II protein complex",0.008,6.7633,0.751,0.501,"plasma membrane region"),
c("GO:0030315","T-tubule",0.008,5.2351,0.840,0.469,"plasma membrane region"),
c("GO:0042383","sarcolemma",0.024,4.1958,0.871,0.281,"plasma membrane region"),
c("GO:0009986","cell surface",0.241,3.3123,0.917,0.057,"cell surface"),
c("GO:0031975","envelope",2.324,9.4471,0.902,0.073,"envelope"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_cc_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
