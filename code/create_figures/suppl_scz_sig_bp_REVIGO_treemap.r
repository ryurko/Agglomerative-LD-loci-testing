

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
revigo.data <- rbind(c("GO:0000003","reproduction",0.769,6.8595,1.000,0.000,"reproduction"),
c("GO:0007610","behavior",0.170,8.1008,0.994,0.000,"behavior"),
c("GO:0009057","macromolecule catabolic process",1.953,11.4550,0.941,0.000,"macromolecule catabolism"),
c("GO:0018193","peptidyl-amino acid modification",1.495,10.3195,0.922,0.138,"macromolecule catabolism"),
c("GO:0006310","DNA recombination",1.641,4.8566,0.923,0.169,"macromolecule catabolism"),
c("GO:0031329","regulation of cellular catabolic process",0.093,5.8941,0.799,0.529,"macromolecule catabolism"),
c("GO:0070647","protein modification by small protein conjugation or removal",0.821,9.8018,0.926,0.525,"macromolecule catabolism"),
c("GO:0043603","cellular amide metabolic process",6.879,7.3954,0.934,0.209,"macromolecule catabolism"),
c("GO:0032446","protein modification by small protein conjugation",0.607,9.8406,0.926,0.508,"macromolecule catabolism"),
c("GO:0031399","regulation of protein modification process",0.565,6.1338,0.771,0.504,"macromolecule catabolism"),
c("GO:0035303","regulation of dephosphorylation",0.046,5.2622,0.821,0.674,"macromolecule catabolism"),
c("GO:0009894","regulation of catabolic process",0.146,7.2349,0.811,0.552,"macromolecule catabolism"),
c("GO:1901361","organic cyclic compound catabolic process",1.164,5.1783,0.943,0.700,"macromolecule catabolism"),
c("GO:0006470","protein dephosphorylation",0.585,8.4862,0.916,0.506,"macromolecule catabolism"),
c("GO:0006479","protein methylation",0.343,5.5838,0.931,0.478,"macromolecule catabolism"),
c("GO:0006508","proteolysis",5.223,9.9625,0.956,0.385,"macromolecule catabolism"),
c("GO:0019953","sexual reproduction",0.275,6.3143,0.975,0.000,"sexual reproduction"),
c("GO:0022008","neurogenesis",0.423,31.7721,0.647,0.000,"neurogenesis"),
c("GO:0007389","pattern specification process",0.147,5.4583,0.738,0.683,"neurogenesis"),
c("GO:0050890","cognition",0.053,7.2974,0.863,0.561,"neurogenesis"),
c("GO:0016050","vesicle organization",0.130,5.1806,0.838,0.558,"neurogenesis"),
c("GO:0016049","cell growth",0.153,5.8420,0.865,0.156,"neurogenesis"),
c("GO:0031344","regulation of cell projection organization",0.123,13.1669,0.667,0.679,"neurogenesis"),
c("GO:0021700","developmental maturation",0.074,5.5495,0.784,0.555,"neurogenesis"),
c("GO:0051240","positive regulation of multicellular organismal process",0.303,11.4948,0.626,0.697,"neurogenesis"),
c("GO:0061024","membrane organization",0.759,12.3344,0.828,0.501,"neurogenesis"),
c("GO:0098930","axonal transport",0.009,6.9944,0.763,0.658,"neurogenesis"),
c("GO:0044839","cell cycle G2/M phase transition",0.054,5.7830,0.845,0.645,"neurogenesis"),
c("GO:0050808","synapse organization",0.070,17.0531,0.794,0.394,"neurogenesis"),
c("GO:0048598","embryonic morphogenesis",0.171,5.9820,0.728,0.692,"neurogenesis"),
c("GO:0051301","cell division",1.230,11.8839,0.856,0.223,"neurogenesis"),
c("GO:0002009","morphogenesis of an epithelium",0.153,7.4175,0.757,0.599,"neurogenesis"),
c("GO:0007267","cell-cell signaling",0.407,22.7107,0.833,0.170,"neurogenesis"),
c("GO:0097435","supramolecular fiber organization",0.345,7.0405,0.769,0.452,"neurogenesis"),
c("GO:0099536","synaptic signaling",0.188,15.5137,0.816,0.427,"neurogenesis"),
c("GO:0007224","smoothened signaling pathway",0.038,4.8079,0.760,0.675,"neurogenesis"),
c("GO:0061061","muscle structure development",0.142,7.0012,0.772,0.595,"neurogenesis"),
c("GO:0044089","positive regulation of cellular component biogenesis",0.193,11.9058,0.697,0.562,"neurogenesis"),
c("GO:0051932","synaptic transmission, GABAergic",0.007,4.8896,0.849,0.689,"neurogenesis"),
c("GO:0043010","camera-type eye development",0.071,6.2258,0.739,0.679,"neurogenesis"),
c("GO:0044087","regulation of cellular component biogenesis",0.404,10.4141,0.733,0.604,"neurogenesis"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",0.279,6.7814,0.725,0.441,"neurogenesis"),
c("GO:0007163","establishment or maintenance of cell polarity",0.099,5.3207,0.881,0.151,"neurogenesis"),
c("GO:0032990","cell part morphogenesis",0.174,14.4771,0.639,0.658,"neurogenesis"),
c("GO:0048285","organelle fission",0.463,12.5611,0.822,0.627,"neurogenesis"),
c("GO:0007612","learning",0.029,5.9702,0.819,0.672,"neurogenesis"),
c("GO:0006520","cellular amino acid metabolic process",5.591,5.2829,0.808,0.607,"neurogenesis"),
c("GO:0061640","cytoskeleton-dependent cytokinesis",0.080,5.9189,0.840,0.667,"neurogenesis"),
c("GO:0033043","regulation of organelle organization",0.495,21.0289,0.676,0.631,"neurogenesis"),
c("GO:1902903","regulation of supramolecular fiber organization",0.166,7.4827,0.669,0.697,"neurogenesis"),
c("GO:0007049","cell cycle",1.885,25.9048,0.850,0.198,"neurogenesis"),
c("GO:0007059","chromosome segregation",0.476,7.4864,0.867,0.201,"neurogenesis"),
c("GO:0048880","sensory system development",0.006,6.3383,0.780,0.566,"neurogenesis"),
c("GO:0006082","organic acid metabolic process",9.086,6.9877,0.811,0.293,"neurogenesis"),
c("GO:0001775","cell activation",0.171,6.1916,0.876,0.158,"neurogenesis"),
c("GO:0007017","microtubule-based process",0.658,16.9711,0.863,0.208,"neurogenesis"),
c("GO:0007010","cytoskeleton organization",0.786,21.7319,0.814,0.490,"neurogenesis"),
c("GO:0007005","mitochondrion organization",0.418,7.1413,0.823,0.621,"neurogenesis"),
c("GO:0022406","membrane docking",0.099,5.5178,0.881,0.151,"neurogenesis"),
c("GO:0044283","small molecule biosynthetic process",5.677,5.8890,0.890,0.608,"neurogenesis"),
c("GO:0060322","head development",0.151,19.5118,0.771,0.599,"neurogenesis"),
c("GO:0022402","cell cycle process",1.053,19.4099,0.812,0.219,"neurogenesis"),
c("GO:0051094","positive regulation of developmental process",0.258,11.2897,0.625,0.652,"neurogenesis"),
c("GO:0034622","cellular macromolecular complex assembly",1.211,14.8400,0.803,0.526,"neurogenesis"),
c("GO:0030029","actin filament-based process",0.398,9.7082,0.868,0.170,"neurogenesis"),
c("GO:0030030","cell projection organization",0.608,24.9024,0.759,0.206,"neurogenesis"),
c("GO:0010256","endomembrane system organization",0.189,7.2346,0.846,0.428,"neurogenesis"),
c("GO:0048732","gland development",0.095,6.0566,0.736,0.695,"neurogenesis"),
c("GO:0000226","microtubule cytoskeleton organization",0.293,13.8002,0.742,0.601,"neurogenesis"),
c("GO:0006915","apoptotic process",0.406,11.1730,0.852,0.170,"neurogenesis"),
c("GO:0022610","biological adhesion",0.550,15.1674,0.994,0.000,"biological adhesion"),
c("GO:0040007","growth",0.317,5.0701,0.994,0.000,"growth"),
c("GO:0042592","homeostatic process",1.661,18.7488,0.771,0.000,"homeostatic process"),
c("GO:0002456","T cell mediated immunity",0.014,5.2670,0.889,0.686,"homeostatic process"),
c("GO:0045815","positive regulation of gene expression, epigenetic",0.014,5.3311,0.796,0.194,"homeostatic process"),
c("GO:1902679","negative regulation of RNA biosynthetic process",0.612,10.6486,0.725,0.274,"homeostatic process"),
c("GO:2001257","regulation of cation channel activity",0.020,5.8553,0.742,0.604,"homeostatic process"),
c("GO:0048002","antigen processing and presentation of peptide antigen",0.013,9.8181,0.918,0.683,"homeostatic process"),
c("GO:0016071","mRNA metabolic process",0.798,9.5318,0.927,0.311,"homeostatic process"),
c("GO:0050852","T cell receptor signaling pathway",0.017,5.1086,0.676,0.693,"homeostatic process"),
c("GO:0042127","regulation of cell proliferation",0.313,5.5099,0.790,0.255,"homeostatic process"),
c("GO:0050803","regulation of synapse structure or activity",0.034,8.1538,0.827,0.541,"homeostatic process"),
c("GO:0001505","regulation of neurotransmitter levels",0.055,4.6228,0.822,0.563,"homeostatic process"),
c("GO:0000375","RNA splicing, via transesterification reactions",0.320,5.9673,0.930,0.593,"homeostatic process"),
c("GO:0032970","regulation of actin filament-based process",0.179,6.3402,0.773,0.242,"homeostatic process"),
c("GO:0044093","positive regulation of molecular function",0.890,14.2418,0.804,0.286,"homeostatic process"),
c("GO:0006259","DNA metabolic process",5.607,9.9503,0.915,0.386,"homeostatic process"),
c("GO:0032535","regulation of cellular component size",0.179,5.7543,0.673,0.627,"homeostatic process"),
c("GO:0002682","regulation of immune system process",0.252,12.4346,0.770,0.250,"homeostatic process"),
c("GO:0099072","regulation of postsynaptic specialization membrane neurotransmitter receptor levels",0.001,6.0665,0.857,0.424,"homeostatic process"),
c("GO:0060249","anatomical structure homeostasis",0.185,7.3918,0.801,0.629,"homeostatic process"),
c("GO:0048872","homeostasis of number of cells",0.086,5.4167,0.811,0.586,"homeostatic process"),
c("GO:0042391","regulation of membrane potential",0.135,6.8646,0.811,0.610,"homeostatic process"),
c("GO:0008380","RNA splicing",0.413,7.2119,0.928,0.292,"homeostatic process"),
c("GO:0006399","tRNA metabolic process",2.495,5.0612,0.918,0.436,"homeostatic process"),
c("GO:0060333","interferon-gamma-mediated signaling pathway",0.004,6.3097,0.725,0.635,"homeostatic process"),
c("GO:0051174","regulation of phosphorus metabolic process",0.580,7.2878,0.787,0.371,"homeostatic process"),
c("GO:0034660","ncRNA metabolic process",3.407,4.7559,0.916,0.376,"homeostatic process"),
c("GO:0009914","hormone transport",0.074,7.5061,0.723,0.578,"homeostatic process"),
c("GO:0046907","intracellular transport",1.564,28.6666,0.824,0.000,"intracellular transport"),
c("GO:0055085","transmembrane transport",8.916,15.7199,0.875,0.497,"intracellular transport"),
c("GO:0006811","ion transport",5.344,13.4196,0.881,0.535,"intracellular transport"),
c("GO:0006839","mitochondrial transport",0.182,6.3781,0.910,0.273,"intracellular transport"),
c("GO:1990778","protein localization to cell periphery",0.063,4.9339,0.855,0.661,"intracellular transport"),
c("GO:0042147","retrograde transport, endosome to Golgi",0.055,6.0725,0.853,0.689,"intracellular transport"),
c("GO:0046903","secretion",0.810,9.8324,0.824,0.319,"intracellular transport"),
c("GO:0035418","protein localization to synapse",0.006,7.2646,0.905,0.197,"intracellular transport"),
c("GO:0016197","endosomal transport",0.131,7.9489,0.913,0.265,"intracellular transport"),
c("GO:0099003","vesicle-mediated transport in synapse",0.036,5.1271,0.911,0.650,"intracellular transport"),
c("GO:0031503","protein complex localization",0.062,6.7215,0.890,0.413,"intracellular transport"),
c("GO:0097120","receptor localization to synapse",0.004,5.1019,0.885,0.535,"intracellular transport"),
c("GO:0002790","peptide secretion",0.041,5.2657,0.841,0.428,"intracellular transport"),
c("GO:0030705","cytoskeleton-dependent intracellular transport",0.056,6.1351,0.862,0.690,"intracellular transport"),
c("GO:0072511","divalent inorganic cation transport",0.393,9.2738,0.887,0.624,"intracellular transport"),
c("GO:0070588","calcium ion transmembrane transport",0.157,8.2279,0.889,0.566,"intracellular transport"),
c("GO:0034220","ion transmembrane transport",3.528,12.9034,0.867,0.381,"intracellular transport"),
c("GO:0048193","Golgi vesicle transport",0.297,6.5609,0.897,0.287,"intracellular transport"),
c("GO:0051668","localization within membrane",0.023,6.4911,0.873,0.608,"intracellular transport"),
c("GO:0051648","vesicle localization",0.078,8.0785,0.858,0.674,"intracellular transport"),
c("GO:0098609","cell-cell adhesion",0.251,17.4067,0.976,0.000,"cell-cell adhesion"),
c("GO:0032259","methylation",3.103,5.4657,0.990,0.019,"methylation"),
c("GO:0006974","cellular response to DNA damage stimulus",2.360,11.4806,0.886,0.035,"cellular response to DNA damage stimulus"),
c("GO:1901701","cellular response to oxygen-containing compound",0.345,5.9628,0.892,0.654,"cellular response to DNA damage stimulus"),
c("GO:1901698","response to nitrogen compound",0.178,5.3403,0.929,0.581,"cellular response to DNA damage stimulus"),
c("GO:1901700","response to oxygen-containing compound",0.503,7.5684,0.923,0.568,"cellular response to DNA damage stimulus"),
c("GO:0035966","response to topologically incorrect protein",0.053,5.1677,0.924,0.660,"cellular response to DNA damage stimulus"),
c("GO:0009628","response to abiotic stimulus",0.571,8.7429,0.936,0.414,"cellular response to DNA damage stimulus"),
c("GO:0070482","response to oxygen levels",0.055,5.2785,0.943,0.652,"cellular response to DNA damage stimulus"),
c("GO:0009719","response to endogenous stimulus",0.526,4.9582,0.936,0.411,"cellular response to DNA damage stimulus"),
c("GO:0034097","response to cytokine",0.136,9.6346,0.925,0.357,"cellular response to DNA damage stimulus"),
c("GO:0071453","cellular response to oxygen levels",0.026,6.0746,0.908,0.455,"cellular response to DNA damage stimulus"),
c("GO:1902531","regulation of intracellular signal transduction",0.547,8.9033,0.698,0.474,"cellular response to DNA damage stimulus"),
c("GO:0006979","response to oxidative stress",0.575,4.9640,0.929,0.661,"cellular response to DNA damage stimulus"),
c("GO:0016311","dephosphorylation",1.250,10.0657,0.944,0.040,"dephosphorylation"),
c("GO:0090407","organophosphate biosynthetic process",4.110,5.8236,0.924,0.457,"dephosphorylation"),
c("GO:0019637","organophosphate metabolic process",6.148,4.8283,0.932,0.579,"dephosphorylation"),
c("GO:0006629","lipid metabolic process",3.522,7.4760,0.901,0.093,"lipid metabolism"),
c("GO:1901566","organonitrogen compound biosynthetic process",14.064,10.5164,0.952,0.286,"lipid metabolism"),
c("GO:0008610","lipid biosynthetic process",2.123,5.1452,0.892,0.218,"lipid metabolism"),
c("GO:0043043","peptide biosynthetic process",5.770,5.6017,0.918,0.587,"lipid metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_bp_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
