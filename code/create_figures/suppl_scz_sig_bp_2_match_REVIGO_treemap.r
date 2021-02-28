

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
revigo.data <- rbind(c("GO:0007049","cell cycle",1.885,18.3842,0.850,0.000,"cell cycle"),
c("GO:0000281","mitotic cytokinesis",0.070,3.8458,0.797,0.659,"cell cycle"),
c("GO:0023061","signal release",0.092,8.5587,0.770,0.402,"cell cycle"),
c("GO:0050890","cognition",0.053,2.9455,0.900,0.561,"cell cycle"),
c("GO:0016049","cell growth",0.153,3.7998,0.865,0.179,"cell cycle"),
c("GO:0051240","positive regulation of multicellular organismal process",0.303,4.9365,0.655,0.697,"cell cycle"),
c("GO:0035567","non-canonical Wnt signaling pathway",0.014,3.9423,0.749,0.691,"cell cycle"),
c("GO:0044839","cell cycle G2/M phase transition",0.054,5.2114,0.803,0.645,"cell cycle"),
c("GO:0051301","cell division",1.230,9.0642,0.856,0.223,"cell cycle"),
c("GO:0002009","morphogenesis of an epithelium",0.153,4.1038,0.786,0.599,"cell cycle"),
c("GO:0061008","hepaticobiliary system development",0.024,3.3035,0.804,0.625,"cell cycle"),
c("GO:0007267","cell-cell signaling",0.407,11.7700,0.838,0.197,"cell cycle"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.064,3.9591,0.951,0.583,"cell cycle"),
c("GO:0061061","muscle structure development",0.142,4.2697,0.804,0.595,"cell cycle"),
c("GO:0046903","secretion",0.810,3.7155,0.839,0.457,"cell cycle"),
c("GO:0043010","camera-type eye development",0.071,3.6533,0.787,0.679,"cell cycle"),
c("GO:0042180","cellular ketone metabolic process",0.423,3.4095,0.857,0.198,"cell cycle"),
c("GO:0022008","neurogenesis",0.423,15.3700,0.693,0.198,"cell cycle"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",0.279,3.1864,0.738,0.441,"cell cycle"),
c("GO:0043043","peptide biosynthetic process",5.770,4.0028,0.914,0.587,"cell cycle"),
c("GO:0007163","establishment or maintenance of cell polarity",0.099,2.9676,0.881,0.172,"cell cycle"),
c("GO:0032990","cell part morphogenesis",0.174,6.7948,0.638,0.658,"cell cycle"),
c("GO:0006520","cellular amino acid metabolic process",5.591,4.2573,0.807,0.384,"cell cycle"),
c("GO:0007059","chromosome segregation",0.476,6.5125,0.866,0.201,"cell cycle"),
c("GO:0048880","sensory system development",0.006,3.1778,0.819,0.566,"cell cycle"),
c("GO:0006082","organic acid metabolic process",9.086,3.3296,0.812,0.608,"cell cycle"),
c("GO:0007017","microtubule-based process",0.658,7.0671,0.863,0.208,"cell cycle"),
c("GO:0022406","membrane docking",0.099,3.7645,0.881,0.172,"cell cycle"),
c("GO:0060322","head development",0.151,8.9738,0.803,0.599,"cell cycle"),
c("GO:0044283","small molecule biosynthetic process",5.677,3.7109,0.896,0.557,"cell cycle"),
c("GO:0022402","cell cycle process",1.053,13.5955,0.774,0.219,"cell cycle"),
c("GO:0030029","actin filament-based process",0.398,4.0944,0.868,0.197,"cell cycle"),
c("GO:0048732","gland development",0.095,3.8084,0.783,0.695,"cell cycle"),
c("GO:0006915","apoptotic process",0.406,5.7756,0.853,0.197,"cell cycle"),
c("GO:0022610","biological adhesion",0.550,4.7930,0.994,0.000,"biological adhesion"),
c("GO:0042592","homeostatic process",1.661,7.8077,0.787,0.000,"homeostatic process"),
c("GO:0060968","regulation of gene silencing",0.029,3.0627,0.761,0.663,"homeostatic process"),
c("GO:0042127","regulation of cell proliferation",0.313,3.2673,0.793,0.255,"homeostatic process"),
c("GO:0010608","posttranscriptional regulation of gene expression",0.719,3.7201,0.798,0.382,"homeostatic process"),
c("GO:0050803","regulation of synapse structure or activity",0.034,3.8509,0.838,0.541,"homeostatic process"),
c("GO:0001505","regulation of neurotransmitter levels",0.055,3.2023,0.833,0.563,"homeostatic process"),
c("GO:0032970","regulation of actin filament-based process",0.179,3.9272,0.780,0.242,"homeostatic process"),
c("GO:0044093","positive regulation of molecular function",0.890,6.4092,0.816,0.286,"homeostatic process"),
c("GO:0033500","carbohydrate homeostasis",0.068,3.5643,0.826,0.574,"homeostatic process"),
c("GO:0016233","telomere capping",0.014,3.0661,0.693,0.504,"homeostatic process"),
c("GO:0032535","regulation of cellular component size",0.179,3.1785,0.663,0.627,"homeostatic process"),
c("GO:0099072","regulation of postsynaptic specialization membrane neurotransmitter receptor levels",0.001,4.7121,0.866,0.424,"homeostatic process"),
c("GO:0060249","anatomical structure homeostasis",0.185,4.9584,0.816,0.629,"homeostatic process"),
c("GO:0040029","regulation of gene expression, epigenetic",0.130,4.1790,0.822,0.325,"homeostatic process"),
c("GO:0009890","negative regulation of biosynthetic process",0.772,6.2723,0.774,0.281,"homeostatic process"),
c("GO:0001678","cellular glucose homeostasis",0.049,3.7998,0.774,0.695,"homeostatic process"),
c("GO:0090630","activation of GTPase activity",0.020,3.0283,0.856,0.604,"homeostatic process"),
c("GO:0009914","hormone transport",0.074,5.6505,0.741,0.578,"homeostatic process"),
c("GO:0046907","intracellular transport",1.564,16.6283,0.827,0.000,"intracellular transport"),
c("GO:0055085","transmembrane transport",8.916,5.4772,0.887,0.440,"intracellular transport"),
c("GO:0032252","secretory granule localization",0.001,3.3034,0.890,0.488,"intracellular transport"),
c("GO:0006811","ion transport",5.344,3.4648,0.893,0.535,"intracellular transport"),
c("GO:0006839","mitochondrial transport",0.182,4.8420,0.919,0.273,"intracellular transport"),
c("GO:0036010","protein localization to endosome",0.005,3.9436,0.874,0.543,"intracellular transport"),
c("GO:1990778","protein localization to cell periphery",0.063,3.6864,0.856,0.661,"intracellular transport"),
c("GO:0042147","retrograde transport, endosome to Golgi",0.055,5.0223,0.858,0.689,"intracellular transport"),
c("GO:0035418","protein localization to synapse",0.006,5.6443,0.905,0.197,"intracellular transport"),
c("GO:0016197","endosomal transport",0.131,5.1597,0.921,0.265,"intracellular transport"),
c("GO:0099003","vesicle-mediated transport in synapse",0.036,3.1962,0.921,0.236,"intracellular transport"),
c("GO:0031503","protein complex localization",0.062,4.7730,0.891,0.413,"intracellular transport"),
c("GO:0099645","neurotransmitter receptor localization to postsynaptic specialization membrane",0.000,3.4341,0.734,0.371,"intracellular transport"),
c("GO:0097120","receptor localization to synapse",0.004,3.3365,0.884,0.535,"intracellular transport"),
c("GO:0030705","cytoskeleton-dependent intracellular transport",0.056,3.4756,0.864,0.690,"intracellular transport"),
c("GO:0072511","divalent inorganic cation transport",0.393,4.1044,0.902,0.689,"intracellular transport"),
c("GO:0070588","calcium ion transmembrane transport",0.157,3.3431,0.906,0.629,"intracellular transport"),
c("GO:0030001","metal ion transport",1.677,4.5686,0.889,0.347,"intracellular transport"),
c("GO:0048193","Golgi vesicle transport",0.297,3.1001,0.909,0.650,"intracellular transport"),
c("GO:0051668","localization within membrane",0.023,4.8649,0.875,0.608,"intracellular transport"),
c("GO:0051648","vesicle localization",0.078,4.8323,0.860,0.674,"intracellular transport"),
c("GO:0048002","antigen processing and presentation of peptide antigen",0.013,7.3559,0.887,0.000,"antigen processing and presentation of peptide antigen"),
c("GO:0002682","regulation of immune system process",0.252,7.1855,0.747,0.683,"antigen processing and presentation of peptide antigen"),
c("GO:0002708","positive regulation of lymphocyte mediated immunity",0.015,3.0035,0.733,0.688,"antigen processing and presentation of peptide antigen"),
c("GO:0051702","interaction with symbiont",0.013,3.8475,0.991,0.000,"interaction with symbiont"),
c("GO:0007276","gamete generation",0.159,3.2301,0.848,0.631,"interaction with symbiont"),
c("GO:0098609","cell-cell adhesion",0.251,7.6907,0.976,0.000,"cell-cell adhesion"),
c("GO:0016071","mRNA metabolic process",0.798,7.8314,0.919,0.037,"mRNA metabolism"),
c("GO:0018193","peptidyl-amino acid modification",1.495,4.4937,0.911,0.525,"mRNA metabolism"),
c("GO:0006310","DNA recombination",1.641,2.9734,0.913,0.659,"mRNA metabolism"),
c("GO:0031329","regulation of cellular catabolic process",0.093,3.2898,0.804,0.529,"mRNA metabolism"),
c("GO:0018205","peptidyl-lysine modification",0.355,4.2896,0.921,0.452,"mRNA metabolism"),
c("GO:0009057","macromolecule catabolic process",1.953,6.1485,0.940,0.128,"mRNA metabolism"),
c("GO:0070647","protein modification by small protein conjugation or removal",0.821,6.8853,0.915,0.146,"mRNA metabolism"),
c("GO:0000375","RNA splicing, via transesterification reactions",0.320,5.9342,0.916,0.593,"mRNA metabolism"),
c("GO:0043603","cellular amide metabolic process",6.879,4.3454,0.929,0.189,"mRNA metabolism"),
c("GO:0043543","protein acylation",0.202,3.4010,0.924,0.429,"mRNA metabolism"),
c("GO:0032446","protein modification by small protein conjugation",0.607,6.1054,0.913,0.477,"mRNA metabolism"),
c("GO:0006259","DNA metabolic process",5.607,5.7242,0.907,0.368,"mRNA metabolism"),
c("GO:0006260","DNA replication",1.577,3.5801,0.909,0.269,"mRNA metabolism"),
c("GO:0009894","regulation of catabolic process",0.146,4.1667,0.818,0.552,"mRNA metabolism"),
c("GO:0008380","RNA splicing",0.413,6.7763,0.917,0.299,"mRNA metabolism"),
c("GO:0006399","tRNA metabolic process",2.495,4.0878,0.910,0.362,"mRNA metabolism"),
c("GO:0051604","protein maturation",0.293,3.1650,0.953,0.259,"mRNA metabolism"),
c("GO:0006470","protein dephosphorylation",0.585,4.6623,0.914,0.475,"mRNA metabolism"),
c("GO:0006479","protein methylation",0.343,3.0314,0.921,0.451,"mRNA metabolism"),
c("GO:0006508","proteolysis",5.223,4.1236,0.948,0.355,"mRNA metabolism"),
c("GO:0034660","ncRNA metabolic process",3.407,3.3329,0.908,0.436,"mRNA metabolism"),
c("GO:0006457","protein folding",0.903,3.1443,0.960,0.038,"protein folding"),
c("GO:0051276","chromosome organization",1.477,16.9281,0.754,0.040,"chromosome organization"),
c("GO:0023056","positive regulation of signaling",0.356,3.0910,0.710,0.686,"chromosome organization"),
c("GO:0034401","chromatin organization involved in regulation of transcription",0.005,4.2716,0.728,0.344,"chromosome organization"),
c("GO:0071824","protein-DNA complex subunit organization",0.238,4.7294,0.793,0.661,"chromosome organization"),
c("GO:0071826","ribonucleoprotein complex subunit organization",0.377,4.6143,0.786,0.690,"chromosome organization"),
c("GO:0016050","vesicle organization",0.130,3.4732,0.795,0.590,"chromosome organization"),
c("GO:0031344","regulation of cell projection organization",0.123,6.2140,0.652,0.679,"chromosome organization"),
c("GO:0022613","ribonucleoprotein complex biogenesis",1.614,3.8887,0.801,0.622,"chromosome organization"),
c("GO:0061024","membrane organization",0.759,6.8037,0.788,0.536,"chromosome organization"),
c("GO:0098930","axonal transport",0.009,5.6329,0.759,0.658,"chromosome organization"),
c("GO:0006325","chromatin organization",0.668,12.0800,0.776,0.528,"chromosome organization"),
c("GO:0050808","synapse organization",0.070,8.1201,0.767,0.425,"chromosome organization"),
c("GO:0097435","supramolecular fiber organization",0.345,3.5334,0.739,0.493,"chromosome organization"),
c("GO:0070841","inclusion body assembly",0.003,4.9562,0.844,0.336,"chromosome organization"),
c("GO:0099563","modification of synaptic structure",0.001,3.4687,0.819,0.592,"chromosome organization"),
c("GO:0044089","positive regulation of cellular component biogenesis",0.193,7.2443,0.661,0.429,"chromosome organization"),
c("GO:0044087","regulation of cellular component biogenesis",0.404,5.1961,0.709,0.510,"chromosome organization"),
c("GO:0099010","modification of postsynaptic structure",0.000,3.5758,0.820,0.295,"chromosome organization"),
c("GO:0030865","cortical cytoskeleton organization",0.035,3.5659,0.802,0.697,"chromosome organization"),
c("GO:0048285","organelle fission",0.463,8.0263,0.776,0.669,"chromosome organization"),
c("GO:0033043","regulation of organelle organization",0.495,12.6425,0.643,0.673,"chromosome organization"),
c("GO:1902903","regulation of supramolecular fiber organization",0.166,5.1760,0.649,0.697,"chromosome organization"),
c("GO:0034723","DNA replication-dependent nucleosome organization",0.003,3.5758,0.834,0.491,"chromosome organization"),
c("GO:0070925","organelle assembly",0.571,4.8459,0.764,0.684,"chromosome organization"),
c("GO:0016569","covalent chromatin modification",0.424,5.0451,0.748,0.698,"chromosome organization"),
c("GO:0007005","mitochondrion organization",0.418,4.5978,0.777,0.662,"chromosome organization"),
c("GO:0051094","positive regulation of developmental process",0.258,5.3665,0.645,0.652,"chromosome organization"),
c("GO:0030030","cell projection organization",0.608,9.7286,0.728,0.523,"chromosome organization"),
c("GO:0010256","endomembrane system organization",0.189,4.1407,0.809,0.465,"chromosome organization"),
c("GO:0000226","microtubule cytoskeleton organization",0.293,6.7205,0.699,0.638,"chromosome organization"),
c("GO:0006974","cellular response to DNA damage stimulus",2.360,7.4112,0.883,0.042,"cellular response to DNA damage stimulus"),
c("GO:1901700","response to oxygen-containing compound",0.503,3.0615,0.921,0.568,"cellular response to DNA damage stimulus"),
c("GO:0035967","cellular response to topologically incorrect protein",0.044,3.7998,0.890,0.652,"cellular response to DNA damage stimulus"),
c("GO:0035966","response to topologically incorrect protein",0.053,4.6853,0.918,0.660,"cellular response to DNA damage stimulus"),
c("GO:0009628","response to abiotic stimulus",0.571,4.0173,0.937,0.414,"cellular response to DNA damage stimulus"),
c("GO:0071322","cellular response to carbohydrate stimulus",0.023,3.0251,0.897,0.652,"cellular response to DNA damage stimulus"),
c("GO:0070482","response to oxygen levels",0.055,3.6779,0.943,0.652,"cellular response to DNA damage stimulus"),
c("GO:0009743","response to carbohydrate",0.041,3.2520,0.923,0.649,"cellular response to DNA damage stimulus"),
c("GO:0034097","response to cytokine",0.136,5.9524,0.918,0.357,"cellular response to DNA damage stimulus"),
c("GO:0071453","cellular response to oxygen levels",0.026,4.8695,0.902,0.455,"cellular response to DNA damage stimulus"),
c("GO:0006979","response to oxidative stress",0.575,3.1060,0.929,0.661,"cellular response to DNA damage stimulus"),
c("GO:0080135","regulation of cellular response to stress",0.182,5.8154,0.768,0.679,"cellular response to DNA damage stimulus"),
c("GO:0070555","response to interleukin-1",0.017,4.5873,0.922,0.611,"cellular response to DNA damage stimulus"),
c("GO:0016311","dephosphorylation",1.250,5.2016,0.948,0.065,"dephosphorylation"),
c("GO:0090407","organophosphate biosynthetic process",4.110,3.3765,0.931,0.457,"dephosphorylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_bp_2_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
