

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
revigo.data <- rbind(c("GO:0003712","transcription cofactor activity",0.250,5.8518,0.949,0.000,"transcription cofactor activity"),
c("GO:0005215","transporter activity",8.494,3.2178,0.975,0.000,"transporter activity"),
c("GO:0015075","ion transmembrane transporter activity",3.726,2.7769,0.921,0.000,"ion transmembrane transporter activity"),
c("GO:0046915","transition metal ion transmembrane transporter activity",0.206,2.6351,0.922,0.630,"ion transmembrane transporter activity"),
c("GO:0005385","zinc ion transmembrane transporter activity",0.025,3.0169,0.928,0.510,"ion transmembrane transporter activity"),
c("GO:0016788","hydrolase activity, acting on ester bonds",4.643,4.3454,0.873,0.000,"hydrolase activity, acting on ester bonds"),
c("GO:0016791","phosphatase activity",1.017,2.7990,0.812,0.659,"hydrolase activity, acting on ester bonds"),
c("GO:0016796","exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters",0.244,2.8228,0.813,0.677,"hydrolase activity, acting on ester bonds"),
c("GO:0016817","hydrolase activity, acting on acid anhydrides",7.223,2.9813,0.869,0.409,"hydrolase activity, acting on ester bonds"),
c("GO:0016887","ATPase activity",4.560,3.0356,0.873,0.379,"hydrolase activity, acting on ester bonds"),
c("GO:0034593","phosphatidylinositol bisphosphate phosphatase activity",0.010,2.6056,0.845,0.410,"hydrolase activity, acting on ester bonds"),
c("GO:0004540","ribonuclease activity",0.612,3.2037,0.812,0.286,"hydrolase activity, acting on ester bonds"),
c("GO:0004177","aminopeptidase activity",0.403,2.9178,0.894,0.272,"hydrolase activity, acting on ester bonds"),
c("GO:0042578","phosphoric ester hydrolase activity",1.364,3.0256,0.817,0.623,"hydrolase activity, acting on ester bonds"),
c("GO:0030234","enzyme regulator activity",0.856,5.0762,0.907,0.000,"enzyme regulator activity"),
c("GO:0050839","cell adhesion molecule binding",0.088,6.9848,0.709,0.000,"cell adhesion molecule binding"),
c("GO:0015631","tubulin binding",0.327,4.3494,0.677,0.531,"cell adhesion molecule binding"),
c("GO:0042393","histone binding",0.087,3.1238,0.709,0.424,"cell adhesion molecule binding"),
c("GO:0045296","cadherin binding",0.062,5.9502,0.715,0.414,"cell adhesion molecule binding"),
c("GO:0051087","chaperone binding",0.101,3.5924,0.706,0.445,"cell adhesion molecule binding"),
c("GO:0005102","receptor binding",0.420,3.3989,0.679,0.573,"cell adhesion molecule binding"),
c("GO:0008022","protein C-terminus binding",0.035,2.6357,0.724,0.398,"cell adhesion molecule binding"),
c("GO:0019902","phosphatase binding",0.032,3.9497,0.675,0.680,"cell adhesion molecule binding"),
c("GO:0019904","protein domain specific binding",0.147,6.0528,0.700,0.440,"cell adhesion molecule binding"),
c("GO:0031072","heat shock protein binding",0.091,2.8071,0.708,0.441,"cell adhesion molecule binding"),
c("GO:0019903","protein phosphatase binding",0.022,3.9423,0.681,0.664,"cell adhesion molecule binding"),
c("GO:0044389","ubiquitin-like protein ligase binding",0.108,5.6849,0.652,0.447,"cell adhesion molecule binding"),
c("GO:0030881","beta-2-microglobulin binding",0.001,2.7526,0.771,0.318,"cell adhesion molecule binding"),
c("GO:0046983","protein dimerization activity",1.157,3.4899,0.656,0.634,"cell adhesion molecule binding"),
c("GO:0046982","protein heterodimerization activity",0.280,3.8738,0.687,0.523,"cell adhesion molecule binding"),
c("GO:0042802","identical protein binding",0.400,5.8037,0.680,0.495,"cell adhesion molecule binding"),
c("GO:0008092","cytoskeletal protein binding",0.708,4.5549,0.667,0.570,"cell adhesion molecule binding"),
c("GO:0043175","RNA polymerase core enzyme binding",0.022,2.9791,0.682,0.663,"cell adhesion molecule binding"),
c("GO:0042605","peptide antigen binding",0.002,4.8058,0.942,0.029,"peptide antigen binding"),
c("GO:0003954","NADH dehydrogenase activity",0.368,3.0613,0.955,0.033,"NADH dehydrogenase activity"),
c("GO:0016746","transferase activity, transferring acyl groups",2.893,3.9344,0.882,0.043,"transferase activity, transferring acyl groups"),
c("GO:0016772","transferase activity, transferring phosphorus-containing groups",8.504,5.3939,0.871,0.421,"transferase activity, transferring acyl groups"),
c("GO:0019787","ubiquitin-like protein transferase activity",0.378,3.1831,0.898,0.270,"transferase activity, transferring acyl groups"),
c("GO:0016410","N-acyltransferase activity",1.120,4.0268,0.851,0.306,"transferase activity, transferring acyl groups"),
c("GO:0034212","peptide N-acetyltransferase activity",0.091,3.6974,0.870,0.642,"transferase activity, transferring acyl groups"),
c("GO:0019199","transmembrane receptor protein kinase activity",0.071,2.7245,0.883,0.615,"transferase activity, transferring acyl groups"),
c("GO:0004672","protein kinase activity",3.390,3.0920,0.856,0.355,"transferase activity, transferring acyl groups"),
c("GO:0004714","transmembrane receptor protein tyrosine kinase activity",0.045,3.1093,0.886,0.526,"transferase activity, transferring acyl groups"),
c("GO:0044877","macromolecular complex binding",0.740,6.1060,0.927,0.044,"macromolecular complex binding"),
c("GO:0008144","drug binding",0.178,5.7932,0.932,0.046,"drug binding"),
c("GO:0003682","chromatin binding",0.220,5.8395,0.879,0.047,"chromatin binding"),
c("GO:0031492","nucleosomal DNA binding",0.008,2.9791,0.873,0.693,"chromatin binding"),
c("GO:0031491","nucleosome binding",0.027,3.4287,0.886,0.679,"chromatin binding"),
c("GO:0030145","manganese ion binding",0.256,2.7992,0.931,0.048,"manganese ion binding"),
c("GO:0001067","regulatory region nucleic acid binding",0.282,4.4758,0.909,0.048,"regulatory region nucleic acid binding"),
c("GO:0030554","adenyl nucleotide binding",14.356,6.9262,0.890,0.134,"regulatory region nucleic acid binding"),
c("GO:0035198","miRNA binding",0.004,3.1185,0.926,0.152,"regulatory region nucleic acid binding"),
c("GO:0032553","ribonucleotide binding",16.747,6.3053,0.889,0.685,"regulatory region nucleic acid binding"),
c("GO:0043565","sequence-specific DNA binding",2.222,3.1554,0.895,0.252,"regulatory region nucleic acid binding"),
c("GO:0003723","RNA binding",5.283,17.5329,0.891,0.362,"regulatory region nucleic acid binding"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf(file="figures/suppl/f_scz_sig_mf_2_match_revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
