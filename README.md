# Reproduce analyses for Gunn et al. 2024
<font size="+1">Follow the steps listed below in the <b><i>Analyses</i></b> section to reproduce analyses for Gunn et al. 2024. Each step below gives a summary of the analysis and directs you to a general code file which then works through the analysis step-by-step. This general file will usually point you to other Rmd code, bash shell scripts, or python scripts. Each analysis is contained within subdirectories of the same name in the main R project directory.</font>

<b>Find the article here</b>: <a href="url">...</a> 

<b>Citation here</b>: Gunn, J. C., Clements, S. J., Adams, G., Sterling, E. M., Moore, M. J., Volkers, T. N., & Eggert, L. S. (2024). Phenotypic homogenization and potential fitness constraints following non-native introgression in an endemic sportfish. <i>Journal of Evolutionary Biology</i>,..., ....

## Project: Effects of admixture on fitness in Neosho Bass populations 
We assessed the effect of admixture on fitness in two stream populations within the native range of the Neosho Bass (<i>M. velox</i>; NB) which are known to have extensively hybridized with Smallmouth Bass (<i>Micropterus dolomieu</i>; SMB). Specifically, we used 14 microsatellite loci in a Bayesian analysis of population structure to estimate proportions of interspecific ancestry in individuals collected from Big Sugar Creek and the Elk River in southwestern Missouri (Central Interior Highlands ecoregion (CIH), North America). We used ancestry inference to estimate the proportion of ancestry derived from SMB and NB. For each individual, we measured age and total length and projected individual growth using the standard parameterization of the von Bertalanffy growth model. Finally, we used body condition as a proxy for fitness and generated an ancestry-condition correlation across the global dataset. We ultimately sought to understand the short-term genetic consequences of admixture for NB populations in order to better inform management and long-term viability of distinct, economically and ecologically important sportfish species in the CIH.

## General information on repository structure
This is a publicly visible GitHub repository storing code (and a small amount of data, although we have done our best to avoid uploading large amounts of data due to the limited storage in GitHub) for Gunn et al. 2024. In the home directory of the repository (SMB_Fitness), you will find a README.md file (the source script for this information), the R Project file (SMB_Fitness.Rproj), a project info folder (project_info, which includes all important information on data procurement for this project), a .gitignore file, and "analysis" directories, each of which corresponds with a specific analysis conducted in our study:

1) 01_map_analysis (generation of map images)
2) 02_filtering_analysis (raw data manipulation, cleaning, and filtering)
3) 03_ancestry_analysis (estimation of genetic ancestry proportions)
4) 04_growth_analysis (von Bertalanffy growth models)
5) 05_ac_analysis (correlation of body condition and ancestry)

Within each analysis directory, you will find an R markdown script (.Rmd) with the name of the analysis, which contains all of the code needed to run the full analysis. Additionally, you will find:

1) code

The code directory will store all source code, shell scripts, lists of bash commands, and software packages needed for analysis. 

Once you have downloaded the repository and located the code directory, you should create two additional sub-directories within each analysis (on the same level as the code directory):

2) data
3) figures

The data directory will store all processed data and metadata needed for analysis. The figures folder will contain any raw figures generated in ggplot for each analysis. Ideally, the Rmd script should have paths set up so that the code reads all data and scripts and generates figures seamlessly.

## Using the code
To reproduce all analyses in Gunn et al. 2024, download this data repository and place in a desired home directory. This may be done on your local machine, but we recommend downloading to a high-performance computing cluster so that all code will run seamlessly in one environment, as long as Rstudio is installed and the GUI can be called on the cluster.

Once all directories are downloaded, create a new sub-directory within the home directory (same level as the five analysis directories, .Rproj, README.md, etc.) called "raw_data". This is where you will store the raw genomic data and associated sample metadata (see <i><b>Data</i></b> section below).

## Data

Raw genotype data, accompanying metadata, and data descriptions are available on Dryad: https://doi.org/10.5061/dryad.xksn02vr6. Data descriptions are also provided in the Rmarkdown files associated with each analysis in this Rproject. Detailed instructions are given for 

Download these data into to your `/raw_data` directory within the home working directory.

You should have 4 new items in the directory: <br>

1.  genotype_data.xlsx <br>
2.  phenotype_data.xlsx <br>
3.  metadata.xlsx <br>
4.  otolith_images.zip <br>

Un-compressed the "otolith_images.zip" folder, which contains all raw otolith images collected for individual fish in this study. You should now have a folder called "otolith_images" in your `/raw_data` directory

If you have any questions or issues with data and/or code, please don't hesitate to contact me: jcgunn@uvm.edu

## Analyses

### Analysis 1: Generating species native range maps
In this analysis, we generated easily readable maps displaying the native distributions of the two species of interest (Smallmouth Bass and Neosho Bass) and the two streams sampled for age and growth analysis (Big Sugar Creek and the Elk River). We generated two types of maps: 1) a map of the Central Interior Highlands (CIH), showing the distributions of Smallmouth Bass and Neosho Bass within the ecoregion and labeling Smallmouth Bass reference populations, and 2) a close-up map of the convergence of Big Sugar Creek and the Elk River, showing approximate sample collection sites. In R, we generated only georeferenced outlines of these maps. Shapes representing stream sites and/or populations were superimposed <i>a posteriori</i> on the maps in PowerPoint.

#### Run the code: `01_map_analysis/smb_fitness_map_analysis.Rmd`

### Analysis 2: Data filtering and summarization
For this analysis, we cleaned, filtered, and summarized the full genotype data (14 microsatellite loci) and phenotypic data (consensus age and total length) for 135 raw samples, which comprised fish collected from two streams (Big Sugar Creek and Elk River) within the recognized native range of the Neosho Bass (<i>M. velox</i>; NB) and a reference set of samples representing Smallmouth Bass (<i>M. dolomieu</i>; SMB) collected from streams in the White River drainage in the Smallmouth Bass native range. We did not collect phenotypic data on the Smallmouth Bass reference samples.

#### Run the code: `02_filtering_analysis/smb_fitness_filtering_analysis.Rmd`

### Analysis 3: Ancestry inference analysis
For this analysis, we use Bayesian clustering analysis in the program STRUCTURE (see citation below in under "Programs Needed") to assess individual proportions of ancestry derived from Smallmouth Bass and Neosho Bass. We analyze all sample fish (obtained from Big Sugar Creek and Elk River within the Neosho Bass native range) together with reference fish (obtained from Crooked Creek, White River, and Tablerock Lake within the Smallmouth Bass native range). Knowing that Smallmouth Bass in the White River drainage are of pure genomic origin (Gunn et al. 2022), we use the minimum ancestry proportion of Smallmouth Bass derived from our microsatellites as a lower bound to identify "pure" vs. "admixed" fish in the Neosho Bass range. Ancestry groups are then used to assess growth and body condition in subsequent analyses.

#### Run the code: `03_ancestry_analysis/smb_fitness_ancestry_analysis.Rmd`

### Analysis 4: von Bertalanffy individual growth modeling
For this analysis, we use a Bayesian hierarchical framework to paramaterize the von Bertalanffy growth model using back-calculated total length and consensus age of samples in the Elk River and Big Sugar Creek to assess the contribution of non-native SMB ancestry to growth. We estimate the linear relationship between population-level growth parameters of the von Bertalanffy model (maximum theoretical total length-at-age, the Brody growth coefficient, and theoretical age at length-0) and ancestry proportion by quantifying average deviations from the global parameters due to SMB vs. NB ancestry, and we account for any potential differences in growth due to sex (male or female) or stream of origin (Big Sugar Creek or Elk River). For all growth models, we include individual ID as a random effect to account for individual variation in linear back-calculation estimates.

#### Run the code: `04_growth_analysis/smb_fitness_growth_analysis.Rmd`

### Analysis 5: Ancestry condition correlation analysis
For this analysis, we assessed the linear relationship (correlation) between interspecific ancestry (as measured by proportion of non-native SMB ancestry due to introgression) and individual physiological condition (a proxy of individual fitness) to determine the effect of admixture on fitness in fish in Big Sugar Creek and Elk River in the NB native range. Ancestry-fitness correlation may be positive (fitness increases with ancestry; heterosis by relief from inbreeding depression, hybrid vigor, or adaptive introgression), negative (fitness decreases with ancestry; outbreeding depression by breakdown of coadapted gene complexes, or introduction of deleterious recessive alleles from a large population), or neutral (no effect of ancestry on fitness; formation of a stable hybrid zone). 

#### Run the code: `05_ac_analysis/smb_fitness_ac_analysis.Rmd`
