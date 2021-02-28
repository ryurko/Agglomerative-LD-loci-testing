

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
revigo.data <- rbind(c("GO:0003700","transcription factor activity, sequence-specific DNA binding",4.217,14.6561,0.985,0.000,"transcription factor activity, sequence-specific DNA binding"),
c("GO:0001216","bacterial-type RNA polymerase transcriptional activator activity, sequence-specific DNA binding",0.000,4.6970,0.984,0.343,"transcription factor activity, sequence-specific DNA binding"),
c("GO:0003712","transcription cofactor activity",0.250,6.7219,0.969,0.000,"transcription cofactor activity"),
c("GO:0005215","transporter activity",8.494,9.3942,0.989,0.000,"transporter activity"),
c("GO:0008144","drug binding",0.178,15.2149,0.953,0.000,"drug binding"),
c("GO:0015075","ion transmembrane transporter activity",3.726,7.5030,0.876,0.000,"ion transmembrane transporter activity"),
c("GO:0022803","passive transmembrane transporter activity",0.699,3.4583,0.900,0.647,"ion transmembrane transporter activity"),
c("GO:0022836","gated channel activity",0.424,4.3434,0.890,0.610,"ion transmembrane transporter activity"),
c("GO:0005385","zinc ion transmembrane transporter activity",0.025,4.7699,0.899,0.510,"ion transmembrane transporter activity"),
c("GO:0015144","carbohydrate transmembrane transporter activity",0.325,3.0154,0.888,0.655,"ion transmembrane transporter activity"),
c("GO:0015136","sialic acid transmembrane transporter activity",0.001,3.1661,0.926,0.384,"ion transmembrane transporter activity"),
c("GO:0072509","divalent inorganic cation transmembrane transporter activity",0.274,3.7481,0.886,0.643,"ion transmembrane transporter activity"),
c("GO:0046915","transition metal ion transmembrane transporter activity",0.206,3.7532,0.886,0.630,"ion transmembrane transporter activity"),
c("GO:0016788","hydrolase activity, acting on ester bonds",4.643,10.0922,0.897,0.000,"hydrolase activity, acting on ester bonds"),
c("GO:0016791","phosphatase activity",1.017,5.3920,0.833,0.659,"hydrolase activity, acting on ester bonds"),
c("GO:0016817","hydrolase activity, acting on acid anhydrides",7.223,7.5015,0.893,0.409,"hydrolase activity, acting on ester bonds"),
c("GO:0047499","calcium-independent phospholipase A2 activity",0.001,2.7351,0.901,0.372,"hydrolase activity, acting on ester bonds"),
c("GO:0008233","peptidase activity",4.049,3.9010,0.898,0.372,"hydrolase activity, acting on ester bonds"),
c("GO:0008238","exopeptidase activity",0.865,2.9272,0.895,0.625,"hydrolase activity, acting on ester bonds"),
c("GO:0004177","aminopeptidase activity",0.403,4.7507,0.895,0.272,"hydrolase activity, acting on ester bonds"),
c("GO:0042578","phosphoric ester hydrolase activity",1.364,6.8356,0.846,0.317,"hydrolase activity, acting on ester bonds"),
c("GO:0070003","threonine-type peptidase activity",0.127,2.2483,0.908,0.527,"hydrolase activity, acting on ester bonds"),
c("GO:0004724","magnesium-dependent protein serine/threonine phosphatase activity",0.002,2.8529,0.882,0.457,"hydrolase activity, acting on ester bonds"),
c("GO:0003724","RNA helicase activity",0.087,2.7891,0.921,0.495,"hydrolase activity, acting on ester bonds"),
c("GO:0016887","ATPase activity",4.560,4.9301,0.894,0.379,"hydrolase activity, acting on ester bonds"),
c("GO:0034593","phosphatidylinositol bisphosphate phosphatase activity",0.010,3.5101,0.864,0.434,"hydrolase activity, acting on ester bonds"),
c("GO:0008408","3'-5' exonuclease activity",0.342,3.3558,0.847,0.698,"hydrolase activity, acting on ester bonds"),
c("GO:0004535","poly(A)-specific ribonuclease activity",0.009,2.8121,0.876,0.430,"hydrolase activity, acting on ester bonds"),
c("GO:0004540","ribonuclease activity",0.612,5.0109,0.846,0.623,"hydrolase activity, acting on ester bonds"),
c("GO:0052866","phosphatidylinositol phosphate phosphatase activity",0.016,2.7580,0.870,0.506,"hydrolase activity, acting on ester bonds"),
c("GO:0046030","inositol trisphosphate phosphatase activity",0.004,2.7351,0.880,0.468,"hydrolase activity, acting on ester bonds"),
c("GO:0030234","enzyme regulator activity",0.856,10.0134,0.910,0.000,"enzyme regulator activity"),
c("GO:0072542","protein phosphatase activator activity",0.003,2.8121,0.929,0.582,"enzyme regulator activity"),
c("GO:0032395","MHC class II receptor activity",0.000,3.9520,0.973,0.000,"MHC class II receptor activity"),
c("GO:0098918","structural constituent of synapse",0.000,3.0402,0.988,0.000,"structural constituent of synapse"),
c("GO:0042605","peptide antigen binding",0.002,5.5822,0.937,0.031,"peptide antigen binding"),
c("GO:0042277","peptide binding",0.089,3.2758,0.940,0.618,"peptide antigen binding"),
c("GO:0003954","NADH dehydrogenase activity",0.368,2.4967,0.972,0.033,"NADH dehydrogenase activity"),
c("GO:0019787","ubiquitin-like protein transferase activity",0.378,6.6414,0.920,0.033,"ubiquitin-like protein transferase activity"),
c("GO:0008276","protein methyltransferase activity",0.257,3.9169,0.902,0.213,"ubiquitin-like protein transferase activity"),
c("GO:0000828","inositol hexakisphosphate kinase activity",0.002,3.1723,0.927,0.376,"ubiquitin-like protein transferase activity"),
c("GO:0061659","ubiquitin-like protein ligase activity",0.094,4.9712,0.928,0.196,"ubiquitin-like protein transferase activity"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.504,10.1843,0.895,0.421,"ubiquitin-like protein transferase activity"),
c("GO:0016278","lysine N-methyltransferase activity",0.087,2.8272,0.906,0.535,"ubiquitin-like protein transferase activity"),
c("GO:0008757","S-adenosylmethionine-dependent methyltransferase activity",0.940,3.8137,0.892,0.655,"ubiquitin-like protein transferase activity"),
c("GO:0016779","nucleotidyltransferase activity",1.954,2.3588,0.890,0.616,"ubiquitin-like protein transferase activity"),
c("GO:0042054","histone methyltransferase activity",0.076,2.9926,0.906,0.529,"ubiquitin-like protein transferase activity"),
c("GO:0019199","transmembrane receptor protein kinase activity",0.071,3.5512,0.898,0.615,"ubiquitin-like protein transferase activity"),
c("GO:0004672","protein kinase activity",3.390,6.3511,0.877,0.275,"ubiquitin-like protein transferase activity"),
c("GO:0004713","protein tyrosine kinase activity",0.138,3.7928,0.902,0.584,"ubiquitin-like protein transferase activity"),
c("GO:0004714","transmembrane receptor protein tyrosine kinase activity",0.045,4.0680,0.897,0.526,"ubiquitin-like protein transferase activity"),
c("GO:0016410","N-acyltransferase activity",1.120,3.0808,0.896,0.313,"ubiquitin-like protein transferase activity"),
c("GO:0008170","N-methyltransferase activity",0.412,3.3256,0.899,0.686,"ubiquitin-like protein transferase activity"),
c("GO:0034212","peptide N-acetyltransferase activity",0.091,2.5320,0.912,0.642,"ubiquitin-like protein transferase activity"),
c("GO:0016746","transferase activity, transferring acyl groups",2.893,5.0543,0.906,0.355,"ubiquitin-like protein transferase activity"),
c("GO:0016741","transferase activity, transferring one-carbon groups",3.043,4.6307,0.905,0.358,"ubiquitin-like protein transferase activity"),
c("GO:0060090","binding, bridging",0.052,2.4891,0.956,0.037,"binding, bridging"),
c("GO:0016875","ligase activity, forming carbon-oxygen bonds",1.067,2.9339,0.971,0.038,"ligase activity, forming carbon-oxygen bonds"),
c("GO:0042802","identical protein binding",0.400,13.2188,0.732,0.044,"identical protein binding"),
c("GO:0008022","protein C-terminus binding",0.035,4.1804,0.769,0.442,"identical protein binding"),
c("GO:0030544","Hsp70 protein binding",0.006,2.8331,0.789,0.393,"identical protein binding"),
c("GO:0061629","RNA polymerase II sequence-specific DNA binding transcription factor binding",0.001,2.2582,0.800,0.354,"identical protein binding"),
c("GO:0042974","retinoic acid receptor binding",0.004,3.7481,0.772,0.535,"identical protein binding"),
c("GO:0031072","heat shock protein binding",0.091,3.6390,0.756,0.476,"identical protein binding"),
c("GO:0005167","neurotrophin TRK receptor binding",0.001,2.7351,0.787,0.512,"identical protein binding"),
c("GO:0001042","RNA polymerase I core binding",0.001,2.7351,0.772,0.563,"identical protein binding"),
c("GO:0005154","epidermal growth factor receptor binding",0.006,2.6546,0.772,0.567,"identical protein binding"),
c("GO:0050839","cell adhesion molecule binding",0.088,8.8974,0.756,0.475,"identical protein binding"),
c("GO:0033130","acetylcholine receptor binding",0.002,2.3735,0.782,0.525,"identical protein binding"),
c("GO:0046983","protein dimerization activity",1.157,8.4193,0.712,0.634,"identical protein binding"),
c("GO:0046982","protein heterodimerization activity",0.280,5.8419,0.736,0.523,"identical protein binding"),
c("GO:0030165","PDZ domain binding",0.013,3.6167,0.781,0.412,"identical protein binding"),
c("GO:0035254","glutamate receptor binding",0.004,4.4873,0.777,0.381,"identical protein binding"),
c("GO:0051787","misfolded protein binding",0.007,2.6370,0.788,0.396,"identical protein binding"),
c("GO:0050811","GABA receptor binding",0.002,3.6782,0.785,0.518,"identical protein binding"),
c("GO:0002039","p53 binding",0.013,2.8944,0.781,0.413,"identical protein binding"),
c("GO:0015631","tubulin binding",0.327,8.1229,0.729,0.531,"identical protein binding"),
c("GO:0044325","ion channel binding",0.017,3.4556,0.778,0.421,"identical protein binding"),
c("GO:0042393","histone binding",0.087,3.3323,0.757,0.475,"identical protein binding"),
c("GO:0045296","cadherin binding",0.062,6.8792,0.761,0.462,"identical protein binding"),
c("GO:0051087","chaperone binding",0.101,3.1423,0.754,0.481,"identical protein binding"),
c("GO:0005102","receptor binding",0.420,9.2648,0.731,0.573,"identical protein binding"),
c("GO:0008134","transcription factor binding",0.320,3.4984,0.736,0.530,"identical protein binding"),
c("GO:0019902","phosphatase binding",0.032,7.7423,0.736,0.680,"identical protein binding"),
c("GO:0036435","K48-linked polyubiquitin binding",0.002,2.4068,0.801,0.364,"identical protein binding"),
c("GO:0019904","protein domain specific binding",0.147,10.8978,0.749,0.495,"identical protein binding"),
c("GO:0019903","protein phosphatase binding",0.022,6.7764,0.735,0.664,"identical protein binding"),
c("GO:0044389","ubiquitin-like protein ligase binding",0.108,7.8967,0.717,0.483,"identical protein binding"),
c("GO:0030881","beta-2-microglobulin binding",0.001,3.5505,0.809,0.345,"identical protein binding"),
c("GO:0043548","phosphatidylinositol 3-kinase binding",0.005,2.3008,0.791,0.388,"identical protein binding"),
c("GO:0005516","calmodulin binding",0.072,3.3674,0.759,0.468,"identical protein binding"),
c("GO:0046875","ephrin receptor binding",0.007,3.4263,0.770,0.553,"identical protein binding"),
c("GO:0086080","protein binding involved in heterotypic cell-cell adhesion",0.001,3.7338,0.808,0.349,"identical protein binding"),
c("GO:0008092","cytoskeletal protein binding",0.708,11.0724,0.722,0.570,"identical protein binding"),
c("GO:0034185","apolipoprotein binding",0.002,2.8529,0.800,0.368,"identical protein binding"),
c("GO:0016874","ligase activity",3.540,3.9648,0.969,0.045,"ligase activity"),
c("GO:0003682","chromatin binding",0.220,7.9103,0.902,0.045,"chromatin binding"),
c("GO:0031491","nucleosome binding",0.027,2.8972,0.909,0.679,"chromatin binding"),
c("GO:0030169","low-density lipoprotein particle binding",0.003,2.2698,0.921,0.584,"chromatin binding"),
c("GO:0001067","regulatory region nucleic acid binding",0.282,10.8009,0.933,0.046,"regulatory region nucleic acid binding"),
c("GO:0097322","7SK snRNA binding",0.001,2.4068,0.941,0.261,"regulatory region nucleic acid binding"),
c("GO:0030554","adenyl nucleotide binding",14.356,16.3342,0.922,0.685,"regulatory region nucleic acid binding"),
c("GO:0048027","mRNA 5'-UTR binding",0.008,2.3806,0.936,0.408,"regulatory region nucleic acid binding"),
c("GO:0070180","large ribosomal subunit rRNA binding",0.068,2.7351,0.929,0.482,"regulatory region nucleic acid binding"),
c("GO:0035198","miRNA binding",0.004,4.5687,0.939,0.152,"regulatory region nucleic acid binding"),
c("GO:0032553","ribonucleotide binding",16.747,16.4466,0.921,0.137,"regulatory region nucleic acid binding"),
c("GO:0003723","RNA binding",5.283,14.8516,0.917,0.362,"regulatory region nucleic acid binding"),
c("GO:0042162","telomeric DNA binding",0.035,2.6904,0.938,0.308,"regulatory region nucleic acid binding"),
c("GO:0043565","sequence-specific DNA binding",2.222,9.9582,0.919,0.252,"regulatory region nucleic acid binding"),
c("GO:0019843","rRNA binding",1.409,2.5972,0.915,0.385,"regulatory region nucleic acid binding"),
c("GO:0044877","macromolecular complex binding",0.740,10.0790,0.948,0.050,"macromolecular complex binding"),
c("GO:0033218","amide binding",0.413,2.9756,0.950,0.050,"amide binding"),
c("GO:0005509","calcium ion binding",0.967,9.3326,0.933,0.054,"calcium ion binding"),
c("GO:0008270","zinc ion binding",4.511,6.3039,0.923,0.409,"calcium ion binding"),
c("GO:0000287","magnesium ion binding",1.785,2.2749,0.930,0.444,"calcium ion binding"),
c("GO:0030145","manganese ion binding",0.256,3.5448,0.937,0.548,"calcium ion binding"),
c("GO:0046914","transition metal ion binding",6.942,7.1799,0.922,0.548,"calcium ion binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_mf_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
