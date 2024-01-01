#Most of this script has been adapted from a script used for analysing differential 
#splicing in Singh et al. 2023 (MBE)

# This script is written in both R and Linux/Bash commands 
# The Linux and R portions are clearly labelled as being separate 

#The script uses a GTF file that can obtained by running the following command on Linux:
# wget https://ftp.ensembl.org/pub/release-105/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.105.gtf.gz

#The QoRTs java utility runs on aligned BAM files that were obtained by aligning the RNA-seq fastqs using bwa-mem,
#followed by some processing (e.g., removal of duplicates, marking reads, etc) using GATK 

## This script requires the R package 'JunctionSeq'. As of Sep 2022, the package has not been updated since
## 2018, and requires older versions of R and DESeq2 to run. Detailed instructions for installing JunctionSeq
## are in the file named "JunctionSeq.Installation.Steps.txt"

##############
#  LINUX    #
##############
#run QoRTs java utility QC on BAM files (need to be sorted by position or query; here, by query)
#QoRTs QC is highly memory-intense; change java specifications to increase the quota of memory; here it's set to 100 GB memory
#runs in parallel over all 48 BAM files, makes an output folder for each BAM 
ls *.bam | parallel -j 24 \
"mkdir {.} && java -Xmx100G -jar ~/apps/hartleys-QoRTs-099881f/QoRTs.jar QC \
                    --stranded --verbose \
                    {} ~/scratch/2nd.Project/Alternative.Splicing/References/Drosophila_melanogaster.BDGP6.32.105.gtf {.}/"

#Make a 'decoder' file that will use the sample names as the directories to find QC outputs
ls ~/scratch/2nd.Project/Alternative.Splicing/QoRTs > QoRTs.decoder.file.txt


############
#   R      #
############
#run QorTs in R to obtain 'size factors'
require(QoRTs)
qorts.results <- read.qc.results.data("~/scratch/2nd.Project/Alternative.Splicing/QoRTs/", 
                 decoder.files = "~/scratch/2nd.Project/Alternative.Splicing/Text.Files/QoRTs.decoder.file.txt", 
                 calc.DESeq2 = TRUE)

get.size.factors(qorts.results, outfile = "~/scratch/2nd.Project/Alternative.Splicing/QoRTs/Normalisation.Size.Factors.txt")


##############
#  LINUX     #
##############
## Create novel junction splice sites using QoRTs MergeNovelSplices
## Inputs: QC outputs folder, decoder file, GTF files and folder name for outputs

java -Xmx10G -jar ~/apps/hartleys-QoRTs-099881f/QoRTs.jar mergeNovelSplices --minCount 20 --stranded --verbose \
~/scratch/2nd.Project/Alternative.Splicing/QoRTs/ \
~/scratch/2nd.Project/Alternative.Splicing/Text.Files/QoRTs.decoder.file.txt \
~/scratch/2nd.Project/Alternative.Splicing/References/Drosophila_melanogaster.BDGP6.32.105.gtf \
~/scratch/2nd.Project/Alternative.Splicing/MergeNovelSplices/ &


############
#   R      #
############

require(DESeq2)
require(JunctionSeq)
require(BiocParallel)

#Set global variables
numCores <- 5
FDRThreshold <- 0.01
mappedReadsThreshold <- 50

#Load in decoder files and add fields for conditions

#Run this only once to edit the decoder file
decoder.file <-  read.delim("~/scratch/2nd.Project/Alternative.Splicing/Text.Files/QoRTs.decoder.file.txt", 
                            header = T, stringsAsFactors = F)
#colnames(decoder.file) <- "sample.ID"
sample.sex <- rep(c("Female", "Female", "Male", "Male"), 6)
decoder.file$sex <-  sample.sex
write.table(decoder.file, "~/scratch/2nd.Project/Alternative.Splicing/Text.Files/QoRTs.decoder.file.for.JunctionSeq.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#Load and subset decoder file to only consider the "body" samples
decoder <- read.table("~/scratch/2nd.Project/Alternative.Splicing/Text.Files/QoRTs.decoder.file.for.JunctionSeq.txt", header = T, stringsAsFactors = F)
decoder.for.junctionseq <- decoder[!(grepl("head", decoder$sample.ID)),]


#Providing the directory for the count files:
countFiles <- paste0("~/scratch/2nd.Project/Alternative.Splicing/MergeNovelSplices/", decoder.for.junctionseq$sample.ID, "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

#Run the differential exon usage (DEU) analysis
jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = decoder.for.junctionseq$sample.ID,
  condition = decoder.for.junctionseq$sex,
  flat.gff.file = "~/scratch/2nd.Project/Alternative.Splicing/MergeNovelSplices/withNovel.forJunctionSeq.gff.gz", 
  nCores = 5,
  verbose=TRUE,
  debug.mode = TRUE
)

# Save output to file
writeCompleteResults(count.set.object, "scratch/2nd.Project/Alternative.Splicing/Differential.Splicing.Outputs/Test.Run.Body/",
                     gzip.output = FALSE,
                     FDR.threshold = FDRThreshold,
                     save.allGenes = TRUE, save.sigGenes = TRUE,
                     save.fit = FALSE, save.VST = FALSE,
                     save.bedTracks = TRUE,
                     save.jscs = TRUE,
                     bedtrack.format = c("BED", "GTF", "GFF3"),
                     verbose = TRUE)