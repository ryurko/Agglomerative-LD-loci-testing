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
revigo.data <- rbind(c("GO:0000003","reproduction",0.780,3.669,1.000,0.000,"reproduction"),
c("GO:0007049","cell cycle",1.745,24.415,0.992,0.000,"cell cycle"),
c("GO:0007610","behavior",0.124,3.556,1.000,0.000,"behavior"),
c("GO:0019882","antigen processing and presentation",0.027,7.922,0.945,0.000,"antigen processing and presentation"),
c("GO:0002377","immunoglobulin production",0.005,4.590,0.926,0.520,"antigen processing and presentation"),
c("GO:0002440","production of molecular mediator of immune response",0.006,4.755,0.927,0.521,"antigen processing and presentation"),
c("GO:0048002","antigen processing and presentation of peptide antigen",0.011,8.174,0.939,0.542,"antigen processing and presentation"),
c("GO:0046649","lymphocyte activation",0.067,5.816,0.935,0.605,"antigen processing and presentation"),
c("GO:0002449","lymphocyte mediated immunity",0.028,4.063,0.941,0.607,"antigen processing and presentation"),
c("GO:0002460","adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",0.032,5.499,0.900,0.612,"antigen processing and presentation"),
c("GO:0022610","biological adhesion",0.553,6.998,1.000,0.000,"biological adhesion"),
c("GO:0040007","growth",0.122,4.084,1.000,0.000,"growth"),
c("GO:0051702","biological process involved in interaction with symbiont",0.019,3.431,0.997,0.000,"biological process involved in interaction with symbiont"),
c("GO:0051726","regulation of cell cycle",0.598,15.667,0.826,0.000,"regulation of cell cycle"),
c("GO:0051960","regulation of nervous system development",0.075,13.726,0.804,0.199,"regulation of cell cycle"),
c("GO:0010975","regulation of neuron projection development",0.066,8.414,0.797,0.209,"regulation of cell cycle"),
c("GO:0022407","regulation of cell-cell adhesion",0.072,3.782,0.843,0.211,"regulation of cell cycle"),
c("GO:0051983","regulation of chromosome segregation",0.076,3.721,0.850,0.212,"regulation of cell cycle"),
c("GO:0099177","regulation of trans-synaptic signaling",0.079,7.593,0.830,0.212,"regulation of cell cycle"),
c("GO:0050865","regulation of cell activation",0.097,5.486,0.848,0.216,"regulation of cell cycle"),
c("GO:0044093","positive regulation of molecular function",0.460,8.721,0.827,0.217,"regulation of cell cycle"),
c("GO:0002682","regulation of immune system process",0.288,11.631,0.842,0.223,"regulation of cell cycle"),
c("GO:0060341","regulation of cellular localization",0.160,15.212,0.799,0.225,"regulation of cell cycle"),
c("GO:0080135","regulation of cellular response to stress",0.182,6.899,0.815,0.228,"regulation of cell cycle"),
c("GO:0032970","regulation of actin filament-based process",0.224,4.362,0.839,0.232,"regulation of cell cycle"),
c("GO:0042592","homeostatic process",1.029,9.502,0.797,0.235,"regulation of cell cycle"),
c("GO:0042127","regulation of cell population proliferation",0.296,5.123,0.835,0.238,"regulation of cell cycle"),
c("GO:0010941","regulation of cell death",0.357,6.983,0.833,0.242,"regulation of cell cycle"),
c("GO:0010564","regulation of cell cycle process",0.361,9.974,0.798,0.242,"regulation of cell cycle"),
c("GO:0044087","regulation of cellular component biogenesis",0.459,6.233,0.830,0.248,"regulation of cell cycle"),
c("GO:1902680","positive regulation of RNA biosynthetic process",0.701,14.164,0.743,0.258,"regulation of cell cycle"),
c("GO:0040029","regulation of gene expression, epigenetic",0.068,4.556,0.846,0.291,"regulation of cell cycle"),
c("GO:0062012","regulation of small molecule metabolic process",0.110,3.640,0.842,0.299,"regulation of cell cycle"),
c("GO:0048024","regulation of mRNA splicing, via spliceosome",0.033,3.539,0.845,0.318,"regulation of cell cycle"),
c("GO:0043484","regulation of RNA splicing",0.057,4.076,0.841,0.332,"regulation of cell cycle"),
c("GO:0009894","regulation of catabolic process",0.358,5.563,0.828,0.332,"regulation of cell cycle"),
c("GO:0031329","regulation of cellular catabolic process",0.311,4.620,0.817,0.337,"regulation of cell cycle"),
c("GO:0051174","regulation of phosphorus metabolic process",0.660,3.798,0.814,0.362,"regulation of cell cycle"),
c("GO:0031399","regulation of protein modification process",0.673,3.487,0.808,0.367,"regulation of cell cycle"),
c("GO:0009890","negative regulation of biosynthetic process",0.914,7.095,0.765,0.396,"regulation of cell cycle"),
c("GO:0010608","posttranscriptional regulation of gene expression",1.580,4.113,0.806,0.398,"regulation of cell cycle"),
c("GO:0099072","regulation of postsynaptic membrane neurotransmitter receptor levels",0.007,5.577,0.854,0.482,"regulation of cell cycle"),
c("GO:0009914","hormone transport",0.028,6.758,0.802,0.535,"regulation of cell cycle"),
c("GO:0010769","regulation of cell morphogenesis involved in differentiation",0.015,3.432,0.831,0.537,"regulation of cell cycle"),
c("GO:0001678","cellular glucose homeostasis",0.041,3.660,0.827,0.551,"regulation of cell cycle"),
c("GO:0050803","regulation of synapse structure or activity",0.044,6.065,0.838,0.554,"regulation of cell cycle"),
c("GO:0048167","regulation of synaptic plasticity",0.028,4.536,0.814,0.598,"regulation of cell cycle"),
c("GO:2001251","negative regulation of chromosome organization",0.055,4.551,0.760,0.615,"regulation of cell cycle"),
c("GO:1902275","regulation of chromatin organization",0.067,3.683,0.802,0.623,"regulation of cell cycle"),
c("GO:0032535","regulation of cellular component size",0.216,3.547,0.746,0.637,"regulation of cell cycle"),
c("GO:0045814","negative regulation of gene expression, epigenetic",0.043,4.264,0.807,0.639,"regulation of cell cycle"),
c("GO:0043087","regulation of GTPase activity",0.157,3.864,0.839,0.642,"regulation of cell cycle"),
c("GO:0010942","positive regulation of cell death",0.101,4.789,0.782,0.647,"regulation of cell cycle"),
c("GO:0002684","positive regulation of immune system process",0.160,8.911,0.743,0.659,"regulation of cell cycle"),
c("GO:0045595","regulation of cell differentiation",0.295,13.275,0.791,0.660,"regulation of cell cycle"),
c("GO:0060284","regulation of cell development",0.087,11.085,0.797,0.668,"regulation of cell cycle"),
c("GO:2001257","regulation of cation channel activity",0.035,3.647,0.799,0.682,"regulation of cell cycle"),
c("GO:0044089","positive regulation of cellular component biogenesis",0.189,7.951,0.769,0.684,"regulation of cell cycle"),
c("GO:0043085","positive regulation of catalytic activity",0.374,4.991,0.826,0.692,"regulation of cell cycle"),
c("GO:1901976","regulation of cell cycle checkpoint",0.012,3.427,0.800,0.695,"regulation of cell cycle"),
c("GO:0016049","cell growth",0.040,5.489,0.984,0.008,"cell growth"),
c("GO:0022406","membrane docking",0.060,3.593,0.994,0.009,"membrane docking"),
c("GO:0007163","establishment or maintenance of cell polarity",0.075,3.697,0.994,0.009,"establishment or maintenance of cell polarity"),
c("GO:0001775","cell activation",0.110,4.418,0.994,0.009,"cell activation"),
c("GO:0098609","cell-cell adhesion",0.215,9.115,0.985,0.010,"cell-cell adhesion"),
c("GO:0006915","apoptotic process",0.218,8.114,0.981,0.010,"apoptotic process"),
c("GO:0007267","cell-cell signaling",0.336,14.275,0.959,0.010,"cell-cell signaling"),
c("GO:0023061","signal release",0.046,9.920,0.914,0.365,"cell-cell signaling"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",0.223,4.149,0.773,0.413,"cell-cell signaling"),
c("GO:0035567","non-canonical Wnt signaling pathway",0.007,4.280,0.809,0.659,"cell-cell signaling"),
c("GO:0022008","neurogenesis",0.356,18.795,0.843,0.010,"neurogenesis"),
c("GO:0050890","cognition",0.046,3.811,0.916,0.535,"neurogenesis"),
c("GO:0021700","developmental maturation",0.062,4.040,0.900,0.569,"neurogenesis"),
c("GO:0007276","gamete generation",0.171,4.552,0.898,0.594,"neurogenesis"),
c("GO:0061061","muscle structure development",0.113,4.261,0.894,0.614,"neurogenesis"),
c("GO:0002009","morphogenesis of an epithelium",0.131,3.683,0.885,0.621,"neurogenesis"),
c("GO:0007612","learning",0.023,3.543,0.917,0.627,"neurogenesis"),
c("GO:0060322","head development",0.148,11.233,0.892,0.628,"neurogenesis"),
c("GO:0061008","hepaticobiliary system development",0.027,4.170,0.882,0.628,"neurogenesis"),
c("GO:0048732","gland development",0.084,4.725,0.869,0.684,"neurogenesis"),
c("GO:0043010","camera-type eye development",0.076,3.432,0.867,0.685,"neurogenesis"),
c("GO:0030029","actin filament-based process",0.368,4.411,0.993,0.010,"actin filament-based process"),
c("GO:0007059","chromosome segregation",0.440,8.279,0.993,0.010,"chromosome segregation"),
c("GO:0007017","microtubule-based process",0.707,8.221,0.993,0.011,"microtubule-based process"),
c("GO:0022402","cell cycle process",0.856,17.221,0.945,0.011,"cell cycle process"),
c("GO:0044839","cell cycle G2/M phase transition",0.013,5.983,0.951,0.579,"cell cycle process"),
c("GO:0044770","cell cycle phase transition",0.056,9.571,0.948,0.649,"cell cycle process"),
c("GO:0006457","protein folding",0.909,3.593,0.993,0.011,"protein folding"),
c("GO:0016071","mRNA metabolic process",1.124,9.807,0.947,0.011,"mRNA metabolic process"),
c("GO:0070647","protein modification by small protein conjugation or removal",1.191,7.782,0.949,0.156,"mRNA metabolic process"),
c("GO:0009057","macromolecule catabolic process",2.208,7.401,0.958,0.169,"mRNA metabolic process"),
c("GO:0043603","cellular amide metabolic process",6.951,4.552,0.961,0.274,"mRNA metabolic process"),
c("GO:0006508","proteolysis",5.739,5.448,0.959,0.356,"mRNA metabolic process"),
c("GO:0098781","ncRNA transcription",0.023,3.855,0.951,0.359,"mRNA metabolic process"),
c("GO:0006310","DNA recombination",1.722,3.602,0.946,0.399,"mRNA metabolic process"),
c("GO:0006479","protein methylation",0.403,3.741,0.951,0.443,"mRNA metabolic process"),
c("GO:0018205","peptidyl-lysine modification",0.449,4.744,0.953,0.448,"mRNA metabolic process"),
c("GO:0006470","protein dephosphorylation",0.671,4.390,0.947,0.467,"mRNA metabolic process"),
c("GO:0008380","RNA splicing",0.610,8.333,0.944,0.483,"mRNA metabolic process"),
c("GO:0032446","protein modification by small protein conjugation",0.925,6.263,0.948,0.483,"mRNA metabolic process"),
c("GO:0006259","DNA metabolic process",5.371,7.046,0.940,0.488,"mRNA metabolic process"),
c("GO:0006354","DNA-templated transcription, elongation",0.088,3.812,0.946,0.497,"mRNA metabolic process"),
c("GO:0018193","peptidyl-amino acid modification",2.741,5.749,0.945,0.549,"mRNA metabolic process"),
c("GO:0006399","tRNA metabolic process",2.675,5.108,0.940,0.570,"mRNA metabolic process"),
c("GO:0000375","RNA splicing, via transesterification reactions",0.426,6.888,0.943,0.597,"mRNA metabolic process"),
c("GO:0034660","ncRNA metabolic process",3.716,4.195,0.940,0.669,"mRNA metabolic process"),
c("GO:0044265","cellular macromolecule catabolic process",1.526,7.026,0.945,0.689,"mRNA metabolic process"),
c("GO:0051301","cell division",1.474,9.715,0.992,0.012,"cell division"),
c("GO:0051276","chromosome organization",2.083,18.735,0.887,0.012,"chromosome organization"),
c("GO:0070841","inclusion body assembly",0.001,5.457,0.936,0.310,"chromosome organization"),
c("GO:0050808","synapse organization",0.059,11.274,0.915,0.424,"chromosome organization"),
c("GO:0034723","DNA replication-dependent nucleosome organization",0.004,3.537,0.914,0.432,"chromosome organization"),
c("GO:0010256","endomembrane system organization",0.233,3.802,0.916,0.483,"chromosome organization"),
c("GO:0097435","supramolecular fiber organization",0.313,3.826,0.914,0.498,"chromosome organization"),
c("GO:0061024","membrane organization",0.560,5.710,0.910,0.529,"chromosome organization"),
c("GO:0030030","cell projection organization",0.664,11.084,0.909,0.540,"chromosome organization"),
c("GO:0032200","telomere organization",0.126,4.185,0.904,0.572,"chromosome organization"),
c("GO:0034622","cellular protein-containing complex assembly",1.206,9.769,0.892,0.578,"chromosome organization"),
c("GO:0006323","DNA packaging",0.208,4.115,0.900,0.601,"chromosome organization"),
c("GO:0048285","organelle fission",0.277,9.818,0.904,0.618,"chromosome organization"),
c("GO:0000226","microtubule cytoskeleton organization",0.298,7.870,0.886,0.623,"chromosome organization"),
c("GO:0007005","mitochondrion organization",0.368,4.005,0.902,0.637,"chromosome organization"),
c("GO:0098930","axonal transport",0.009,5.929,0.905,0.650,"chromosome organization"),
c("GO:0022613","ribonucleoprotein complex biogenesis",1.879,4.528,0.905,0.680,"chromosome organization"),
c("GO:0006325","chromatin organization",0.761,13.224,0.890,0.690,"chromosome organization"),
c("GO:0007010","cytoskeleton organization",0.790,9.656,0.896,0.693,"chromosome organization"),
c("GO:0030865","cortical cytoskeleton organization",0.035,3.508,0.910,0.697,"chromosome organization"),
c("GO:0046907","intracellular transport",1.822,17.946,0.895,0.012,"intracellular transport"),
c("GO:0035418","protein localization to synapse",0.007,6.470,0.946,0.192,"intracellular transport"),
c("GO:0097120","receptor localization to synapse",0.008,3.513,0.963,0.194,"intracellular transport"),
c("GO:0031503","protein-containing complex localization",0.175,5.191,0.958,0.251,"intracellular transport"),
c("GO:0006839","mitochondrial transport",0.212,4.582,0.957,0.266,"intracellular transport"),
c("GO:0048870","cell motility",0.572,4.039,0.937,0.284,"intracellular transport"),
c("GO:0046903","secretion",0.631,6.067,0.952,0.298,"intracellular transport"),
c("GO:0006812","cation transport",3.432,5.520,0.936,0.366,"intracellular transport"),
c("GO:0006403","RNA localization",0.203,3.796,0.939,0.408,"intracellular transport"),
c("GO:0060401","cytosolic calcium ion transport",0.031,3.517,0.953,0.494,"intracellular transport"),
c("GO:0055085","transmembrane transport",14.451,7.472,0.927,0.509,"intracellular transport"),
c("GO:0036010","protein localization to endosome",0.003,3.897,0.928,0.510,"intracellular transport"),
c("GO:0006811","ion transport",5.089,4.540,0.942,0.549,"intracellular transport"),
c("GO:1990778","protein localization to cell periphery",0.056,4.373,0.915,0.636,"intracellular transport"),
c("GO:0030705","cytoskeleton-dependent intracellular transport",0.071,3.460,0.918,0.684,"intracellular transport"),
c("GO:0042147","retrograde transport, endosome to Golgi",0.079,3.721,0.913,0.691,"intracellular transport"),
c("GO:0051656","establishment of organelle localization",0.157,6.170,0.909,0.699,"intracellular transport"),
c("GO:0006974","cellular response to DNA damage stimulus",2.542,8.391,0.914,0.013,"cellular response to DNA damage stimulus"),
c("GO:0071453","cellular response to oxygen levels",0.022,5.606,0.924,0.326,"cellular response to DNA damage stimulus"),
c("GO:0009719","response to endogenous stimulus",0.539,3.652,0.942,0.401,"cellular response to DNA damage stimulus"),
c("GO:0009628","response to abiotic stimulus",0.598,6.102,0.942,0.406,"cellular response to DNA damage stimulus"),
c("GO:1901700","response to oxygen-containing compound",0.566,4.942,0.921,0.484,"cellular response to DNA damage stimulus"),
c("GO:0035966","response to topologically incorrect protein",0.066,6.895,0.920,0.537,"cellular response to DNA damage stimulus"),
c("GO:0043620","regulation of DNA-templated transcription in response to stress",0.021,3.700,0.808,0.555,"cellular response to DNA damage stimulus"),
c("GO:0070555","response to interleukin-1",0.017,4.921,0.925,0.570,"cellular response to DNA damage stimulus"),
c("GO:0035967","cellular response to topologically incorrect protein",0.053,5.669,0.909,0.614,"cellular response to DNA damage stimulus"),
c("GO:1901701","cellular response to oxygen-containing compound",0.406,4.222,0.907,0.615,"cellular response to DNA damage stimulus"),
c("GO:0070482","response to oxygen levels",0.045,4.809,0.948,0.632,"cellular response to DNA damage stimulus"),
c("GO:0034097","response to cytokine",0.171,6.882,0.918,0.667,"cellular response to DNA damage stimulus"),
c("GO:0006979","response to oxidative stress",0.590,3.773,0.934,0.671,"cellular response to DNA damage stimulus"),
c("GO:0070498","interleukin-1-mediated signaling pathway",0.006,4.263,0.796,0.687,"cellular response to DNA damage stimulus"),
c("GO:0042180","cellular ketone metabolic process",0.447,4.020,0.978,0.065,"cellular ketone metabolic process"),
c("GO:0006520","cellular amino acid metabolic process",5.778,5.426,0.961,0.375,"cellular ketone metabolic process"),
c("GO:0006082","organic acid metabolic process",9.322,4.386,0.970,0.592,"cellular ketone metabolic process"),
c("GO:0044283","small molecule biosynthetic process",6.682,4.115,0.975,0.609,"cellular ketone metabolic process"),
c("GO:0016311","dephosphorylation",1.604,5.145,0.975,0.075,"dephosphorylation"),
c("GO:0090407","organophosphate biosynthetic process",4.883,3.460,0.962,0.465,"dephosphorylation"),
c("GO:0019637","organophosphate metabolic process",6.558,3.534,0.970,0.579,"dephosphorylation"),
c("GO:0015980","energy derivation by oxidation of organic compounds",1.365,3.528,0.977,0.076,"energy derivation by oxidation of organic compounds"),
c("GO:1901566","organonitrogen compound biosynthetic process",15.347,3.780,0.961,0.401,"organonitrogen compound biosynthetic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_bp_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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

