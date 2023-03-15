# Reproduce analyses for Gunn et al. [DATE]
<font size="+1">Follow the steps listed below in the <b><i>Analyses</i></b> section to reproduce analyses for Gunn et al. (DATE). Each step below gives a summary of the analysis and directs you to a general code file which then works through the analysis step-by-step. This general file will usually point you to other Rmd code, bash shell scripts, or python scripts. Each analysis is contained within subdirectories of the same name in the main R project directory.</font>

## Project: Effects of admixture on fitness in Neosho Bass populations 
We assessed the effect of admixture on fitness in two stream populations within the native range of the Neosho Bass (<i>M. velox</i>) which are known to have extensively hybridized with Smallmouth Bass (<i>Micropterus dolomieu</i>). Specifically, we used 14 microsatellite loci in a Bayesian analysis of population structure to estimate proportions of interspecific ancestry in individuals collected from Big Sugar Creek and the Elk River in southwestern Missouri (Central Interior Highlands ecoregion (CIH), North America). We used ancestry inference to identify fish as "Pure Neosho Bass", "Pure Smallmouth Bass", or "Admixed". For each group, we measured age and total length and projected individual growth using the standard parameterization of the von Berlanffy growth model, comparing average theoretical maximum length among groups. Finally, we used calculated a body condition as a proxy of fitness and generated heterozygosity-fitness correlations of body condition across the global dataset, within stream populations, and within ancestry groups. We ultimately sought to understand the short-term genetic consequences of admixture for Neosho Bass populations in order to better inform management and long-term viability of distinct, economically and ecologically important sportfish species in the CIH.

## General information on repository structure
This is a publicly visible GitHub repository storing code (and a small amount of data, although we have done our best to avoid uploading large amounts of data due to the limited storage ing GitHub) for Gunn et al. (DATE). In the home directory of the repository (SMB_Fitness), you will find a README.md file (the source script for this information), the R Project file (SMB_Fitness.Rproj), a project info folder (project_info, which includes all important information on data procurement for this project), a .gitignore file, and "analysis" directories, each of which corresponds with a specific analysis conducted in our study:

1) map_analysis
2) filtering_analysis
3) ancestry_analysis
4) growth_analysis
5) hfc_analysis

Within each analysis directory, you will find an R markdown script (.Rmd) with the name of the analysis, which contains all of the code needed to run the full analysis. Additionally, you will find:

1) code

The code directory will store all source code, shell scripts, lists of bash commands, and software packages needed for analysis. 

Once you have downloaded the repository and located the code directory, you should create two additional sub-directories within each analysis (on the same level as the code directory):

2) data
3) figures

The data directory will store all processed data and metadata needed for analysis. The figures folder will contain any raw figures generated in ggplot for each analysis. Ideally, the Rmd script should have paths set up so that the code reads all data and scripts and generates figures seamlessly.

## Using the code
To reproduce all analyses in Gunn et al. (DATE), download this data repository and place in a desired home directory. This may be done on your local machine, but we recommend downloading to a high-performance computing cluster so that all code will run seamlessly in one environment, as long as Rstudio is installed and the GUI can be called on the cluster.

Once all directories are downloaded, create a new sub-directory within the home directory (same level as the seven analysis directories, .Rproj, README.md, etc.) called "raw_data". This is where you will store the raw genomic data and associated sample metadata (see <i><b>Data</i></b> section below).

## Data

Raw genotype data and accompanying metadata are available at Zenodo.org: [LINK]

Download these data into to your `/raw_data` directory within the home working directory.

You should have 2 new items in the directory: <br>

1.  <br>
2.  <br>

If you have any questions or issues with data and/or code, please don't hesitate to contact me: jcgunn@uvm.edu

## Analyses

### Analysis 1: Generating species native range maps
In this analysis, we generated easily readable maps displaying the native distributions of the two species of interest (Smallmouth Bass and Neosho Bass) and the two streams sampled for age and growth analysis (Big Sugar Creek and the Elk River). We generated two types of maps: 1) a map of the Central Interior Highlands (CIH), showing the distributions of Smallmouth Bass and Neosho Bass within the ecoregion and labeling Smallmouth Bass reference populations, and 2) a close-up map of the convergence of Big Sugar Creek and the Elk River, showing approximate sample collection sites. In R, we generated only georeferenced outlines of these maps. Shapes representing stream sites and/or populations were superimposed <i>a posteriori</i> on the maps in PowerPoint.

#### Run the code: `map_analysis/smb_fitness_map_analysis.Rmd`

### Analysis 2: Data filtering and summarization
For this aim, we cleaned, filtered, and summarized the full genotype data (14 microsatellite loci) and phenotypic data (consensus age and total length) for 135 raw samples, which comprised fish collected from two streams (Big Sugar Creek and Elk River) within the recognized native range of the Neosho Bass (<i>M. velox</i>; NB) and a reference set of samples representing Smallmouth Bass (<i>M. dolomieu</i>; SMB) collected from streams in the White River drainage in the Smallmouth Bass native range. We did not collect phenotypic data on the Smallmouth Bass reference samples.

#### Run the code: `filtering_analysis/smb_fitness_filtering_analysis.Rmd`

### Analysis 3: Ancestry inference with Bayesian clustering in STRUCTURE


#### Run the code: `ancestry_analysis/smb_fitness_ancestry_analysis.Rmd`

### Analysis 4: von Bertalanffy individual growth modeling


#### Run the code: `growth_analysis/smb_fitness_growth_analysis.Rmd`

### Analysis 5: Heterozygosity fitness correlation analysis


#### Run the code: `hfc_analysis/smb_fitness_hfc_analysis.Rmd`



