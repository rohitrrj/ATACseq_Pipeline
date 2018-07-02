#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=10G
#$ -pe shm 6
#$ -l h_rt=6:00:00
#$ -A projects_goronzy
#$ -m ea
#$ -M jadhav@stanford.edu

#print the time and date
echo Starting: $(date)

#-------------------------------#
#initialize job context         #
#-------------------------------#

module load python
module load fastqc/0.11.2
module load bowtie/2.2.4
module load MACS2/2.1.1
module load java/latest  python samtools/1.3  picard-tools/1.92  preseq/1.0.2 bedtools/2.25.0 r/3.2.2 igvtools/2.3.3  preseq/1.0.2 ngsplot/2.47

myDATADIR="/srv/gsfs0/projects/goronzy/Huimin/raw_data"
myPROJDIR="/srv/gsfs0/projects/goronzy/Huimin"
mySCRIPTSDIR="/srv/gsfs0/projects/goronzy/RA/scripts"
mySampleFile="${myPROJDIR}/Run1.txt"

##-----------------------------##
##create FolderLayout          ##
##-----------------------------##
function createFolderLayout()
{
###be aware of the input parameter
###${1} --- means the ${jobSample}

myAnalyDIR="${myPROJDIR}/Analysis"
if [ ! -d "${myAnalyDIR}" ]
then
echo "${myAnalyDIR} does not exist, creating it now."
tmpCMDSTR="mkdir ${myAnalyDIR}"
eval "${tmpCMDSTR}"
else
echo "${myAnalyDIR} already exist, just use it."
fi


mytrimDIR="${myAnalyDIR}/trim"
if [ ! -d "${mytrimDIR}" ]
then
echo "${mytrimDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${mytrimDIR}; mkdir ${mytrimDIR}/${1}"
tmpCMDSTR="mkdir ${mytrimDIR}; mkdir ${mytrimDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${mytrimDIR} already exist, just use it."
#####for each sample#####
#sampletrimDir=${mytrimDIR}/${1}
sampletrimDir=${mytrimDIR}/${jobSample}
if [ ! -d "${sampletrimDir}" ]
then
echo "${sampletrimDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${sampletrimDir}"
eval "${tmpCMDSTR}"
else
echo "${sampletrimDir} already exist, just use it."
fi
fi

myfastqcDIR="${myAnalyDIR}/fastqc"
if [ ! -d "${myfastqcDIR}" ]
then
echo "${myfastqcDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myfastqcDIR}; mkdir ${myfastqcDIR}/${1}"
tmpCMDSTR="mkdir ${myfastqcDIR}; mkdir ${myfastqcDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myfastqcDIR} already exist, just use it."
#####for each sample#####
#samplefastqcDir=${myfastqcDIR}/${1}
samplefastqcDir=${myfastqcDIR}/${jobSample}
if [ ! -d "${samplefastqcDir}" ]
then
echo "${samplefastqcDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplefastqcDir}"
eval "${tmpCMDSTR}"
else
echo "${samplefastqcDir} already exist, just use it."
fi
fi

mybowtieDIR="${myAnalyDIR}/bowtie"
if [ ! -d "${mybowtieDIR}" ]
then
echo "${mybowtieDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${mybowtieDIR}; mkdir ${mybowtieDIR}/${1}"
tmpCMDSTR="mkdir ${mybowtieDIR}; mkdir ${mybowtieDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${mybowtieDIR} already exist, just use it."
#####for each sample#####
#samplebowtieDir=${mybowtieDIR}/${1}
samplebowtieDir=${mybowtieDIR}/${jobSample}
if [ ! -d "${samplebowtieDir}" ]
then
echo "${samplebowtieDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplebowtieDir}"
eval "${tmpCMDSTR}"
else
echo "${samplebowtieDir} already exist, just use it."
fi
fi

myQCDIR="${myAnalyDIR}/QC"
if [ ! -d "${myQCDIR}" ]
then
echo "${myQCDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myQCDIR}; mkdir ${myQCDIR}/${1}"
tmpCMDSTR="mkdir ${myQCDIR}; mkdir ${myQCDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myQCDIR} already exist, just use it."
#####for each sample#####
#sampleQCDir=${myQCDIR}/${1}
sampleQCDir=${myQCDIR}/${jobSample}
if [ ! -d "${sampleQCDir}" ]
then
echo "${sampleQCDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${sampleQCDir}"
eval "${tmpCMDSTR}"
else
echo "${sampleQCDir} already exist, just use it."
fi
fi

myPeaksDIR="${myAnalyDIR}/Peaks"
if [ ! -d "${myPeaksDIR}" ]
then
echo "${myPeaksDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myPeaksDIR}; mkdir ${myPeaksDIR}/${1}"
tmpCMDSTR="mkdir ${myPeaksDIR}; mkdir ${myPeaksDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myPeaksDIR} already exist, just use it."
#####for each sample#####
#samplePeaksDir=${myPeaksDIR}/${1}
samplePeaksDir=${myPeaksDIR}/${jobSample}
if [ ! -d "${samplePeaksDir}" ]
then
echo "${samplePeaksDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplePeaksDir}"
eval "${tmpCMDSTR}"
else
echo "${samplePeaksDir} already exist, just use it."
fi
fi

}

##-----------------------------##
##run Trimming                 ##
##-----------------------------##
function runTrim()
{

tempRead1="${myDATADIR}/${jobSample}_R1.fastq"
tempRead2="${myDATADIR}/${jobSample}_R2.fastq"

trimCMD="python ${mySCRIPTSDIR}/pyadapter_trim.py -a ${tempRead1} -b ${tempRead2}"
echo ${trimCMD}
eval ${trimCMD}
}

##-----------------------------##
##move Trim Files              ##
##-----------------------------##
function mvTrimFile()
{
myTrimDIR="${myPROJDIR}/qsub_scripts"

tempTrimFile="${myTrimDIR}/${jobSample}*.trim.fastq"

mvTrimFileCMD="mv ${tempTrimFile} ${myPROJDIR}/Analysis/trim/${jobSample}/"
echo ${mvTrimFileCMD}
eval ${mvTrimFileCMD}

}

##-----------------------------##
##run Fastqc                   ##
##-----------------------------##
function runfastqc()
{

tempRead1="${myDATADIR}/${jobSample}_R1.fastq"
tempRead2="${myDATADIR}/${jobSample}_R2.fastq"

fastqcCMD="fastqc -o ${myPROJDIR}/Analysis/fastqc/${jobSample} -f fastq ${tempRead1} ${tempRead2}"
echo ${fastqcCMD}
eval ${fastqcCMD}

}

##-----------------------------##
##run Bowtie                   ##
##-----------------------------##
function runbowtie()
{
myTrimDIR="${myPROJDIR}/Analysis/trim/${jobSample}"

tempRead1="${myTrimDIR}/${jobSample}_R1.trim.fastq"
tempRead2="${myTrimDIR}/${jobSample}_R2.trim.fastq"

bowtieCMD="bowtie2 --no-discordant  --no-dovetail  --no-mixed  --no-unal  --very-sensitive  -X 2000 -p 25 -x /srv/gsfs0/shared_data/RefGenomes/H_sapiens/hg19/indexes/bowtie/2.2.4/hg19 --rg-id ${jobSample} -1 ${tempRead1} -2 ${tempRead2} -S ${myPROJDIR}/Analysis/bowtie/${jobSample}/${jobSample}.sam"
echo ${bowtieCMD}
eval ${bowtieCMD}
}

##-----------------------------##
##run QualityControl           ##
##-----------------------------##
function runQC()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#converts a .sam to .bam file (more compressed) along with removing reads with mapping quality < 20
samtobamCMD="samtools view -q 20 -uT /srv/gsfs0/shared_data/RefGenomes/H_sapiens/hg19/hg19.fa ${myBowtieDir}/${jobSample}.sam -o ${myQCDIR}/${jobSample}.bam"
echo ${samtobamCMD}
eval ${samtobamCMD}

#sorts file to be in order; samples need to be sorted before they can be merged
bamsortCMD="samtools sort -@ 6  -T ${myQCDIR}/tmp  -o ${myQCDIR}/${jobSample}.sorted.bam ${myQCDIR}/${jobSample}.bam"
echo ${bamsortCMD}
eval ${bamsortCMD}

#generates an index file for later processing step
bamindexCMD="samtools index ${myQCDIR}/${jobSample}.sorted.bam"
echo ${bamindexCMD}
eval ${bamindexCMD}

#saves reads from chrom 1-22 in a file (Remove chrM, chrX, chrY)
rmchrCMD="samtools view -u ${myQCDIR}/${jobSample}.sorted.bam `for c in $(seq 1 22); do echo -n "chr${c} "; done` | samtools rmdup - ${myQCDIR}/${jobSample}.rmdup.bam"
echo ${rmchrCMD}
eval ${rmchrCMD}

#removes duplicate reads
#rmdupCMD="samtools rmdup ${myQCDIR}/${jobSample}.sorted.bam ${myQCDIR}/${jobSample}.rmdup.bam"
#echo ${rmdupCMD}
#eval ${rmdupCMD}

#Sort reads again
sortrmdup="samtools sort -@ 6 -T ${myQCDIR}/tmp -o ${myQCDIR}/${jobSample}.rmdup.sorted.bam ${myQCDIR}/${jobSample}.rmdup.bam"
echo ${sortrmdup}
eval ${sortrmdup}

#generates index for later processing steps
rmdupindexCMD="samtools index ${myQCDIR}/${jobSample}.rmdup.sorted.bam"
echo ${rmdupindexCMD}
eval ${rmdupindexCMD}

#converts the bam file to a bedgraph file
bamtobedgraphCMD="genomeCoverageBed -bg -split -ibam ${myQCDIR}/${jobSample}.rmdup.sorted.bam -trackopts 'type=bedGraph name='${jobSample} -g /srv/gsfs0/shared_data/RefGenomes/H_sapiens/hg19/hg19.fa.fai > ${myQCDIR}/${jobSample}.rmdup.bedGraph"
echo ${bamtobedgraphCMD}
eval ${bamtobedgraphCMD}

# convert to tdf
bamtotdf="igvtools count ${myQCDIR}/${jobSample}.rmdup.sorted.bam ${myQCDIR}/${jobSample}.rmdup.tdf /srv/gsfs0/shared_data/RefGenomes/H_sapiens/hg19/hg19.fa.fai"
echo ${bamtotdf}
eval ${bamtotdf}

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
#vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/parsed_hg19_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
#echo ${vplotCMD}
#eval ${vplotCMD}

#generates fragement plot, also a measurement of openness
#insertsizeCMD="java -Xmx8g -jar /srv/gsfs0/software/picard-tools/1.92/CollectInsertSizeMetrics.jar INPUT=${myQCDIR}/${jobSample}.rmdup.sorted.bam OUTPUT=${myQCDIR}/${jobSample}.InsertSizeMetrics.txt HISTOGRAM_FILE=${myQCDIR}/${jobSample}.InsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT"
#echo ${insertsizeCMD}
#eval ${insertsizeCMD}

#estimates library complexity
#preseqcurveCMD="preseq c_curve -o ${myQCDIR}/${jobSample}.c.curve.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
#echo ${preseqcurveCMD}
#eval ${preseqcurveCMD}
#preseqextrapCMD="preseq lc_extrap -o ${myQCDIR}/${jobSample}.lc.extrap.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
#echo ${preseqextrapCMD}
#eval ${preseqextrapCMD}
}

##-----------------------------##
##run QC Metrics Seprately     ##
##-----------------------------##
function runQCMetrics()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/parsed_hg19_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
echo ${vplotCMD}
eval ${vplotCMD}

#generates fragement plot, also a measurement of openness
insertsizeCMD="java -Xmx8g -jar /srv/gsfs0/software/picard-tools/1.92/CollectInsertSizeMetrics.jar INPUT=${myQCDIR}/${jobSample}.rmdup.sorted.bam OUTPUT=${myQCDIR}/${jobSample}.InsertSizeMetrics.txt HISTOGRAM_FILE=${myQCDIR}/${jobSample}.InsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT"
echo ${insertsizeCMD}
eval ${insertsizeCMD}

#estimates library complexity
preseqcurveCMD="preseq c_curve -o ${myQCDIR}/${jobSample}.c.curve.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
echo ${preseqcurveCMD}
eval ${preseqcurveCMD}
preseqextrapCMD="preseq lc_extrap -o ${myQCDIR}/${jobSample}.lc.extrap.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
echo ${preseqextrapCMD}
eval ${preseqextrapCMD}
}

##-----------------------------##
##run V-Plot Seprately     ##
##-----------------------------##
function runVPlot()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/parsed_hg19_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
echo ${vplotCMD}
eval ${vplotCMD}
}

##-----------------------------##
##run MACS (Peak calling)      ##
##-----------------------------##
function runMACS()
{
myPeakDIR="${myPROJDIR}/Analysis/Peaks/${jobSample}"

MACSCMD="macs2 callpeak -t ${myPROJDIR}/Analysis/QC/${jobSample}/${jobSample}.rmdup.sorted.bam -f BAM --keep-dup all -n ${jobSample} --nomodel --shift -100 --extsize 200 -q .01 --nolambda --call-summits --outdir ${myPeakDIR}"
echo ${MACSCMD}
eval ${MACSCMD}
}

##-----------------------------##
##run Gunzip      ##
##-----------------------------##
function runGzip()
{
GZIPCMD1="gunzip ${myDATADIR}/${jobSample}_R1.fastq.gz"
GZIPCMD2="gunzip ${myDATADIR}/${jobSample}_R2.fastq.gz"
echo ${GZIPCMD1}
eval ${GZIPCMD1}
echo ${GZIPCMD2}
eval ${GZIPCMD2}
}


#-------------------------------#
#starting job specification     #
#-------------------------------#
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray=($(cat ${mySampleFile}))' #sample-index mapping

for i in "${mySampleArray[@]}"
do
jobSample=$i
echo ${jobSample}
createFolderLayout ${jobSample}
#runGzip ${jobSample}
#runTrim ${jobSample}
mvTrimFile ${jobSample}
runfastqc ${jobSample}
runbowtie ${jobSample}
runQC ${jobSample}
runMACS ${jobSample}
runQCMetrics ${jobSample}
#runVPlot ${jobSample}

done

#print the time and date again
echo Ending: $(date)
