task getSampleId {
  File inBam

  command {
     samtools view -H ${inBam} | grep "SM:" | sed 's/.*SM:\(.*\)\t.*/\1/g'
  }

  output {
    String SM = read_string(stdout())
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task compareGenotype {
  File controlBam
  File tumorBam
  String baseName
  String outputDir

  command {
    compareBamGenotypes.pl \
    -o ${outputDir + "/genotype"} \
    -nb ${controlBam} \
    -j ${outputDir + "/genotype/" + baseName + "_summary.json"} \
    -tb ${tumorBam}
  }

  output {
    File genotype = "${outputDir}/genotype/${baseName}_summary.json"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


task analyzeContamination {
  File bamFile
  Array[File]? ascatSegmentFiles
  Int contamDownsampOneIn = 25
  String process
  String baseName
  String outputDir

  command <<<
    if [ ${process} == "tumor" ]; then
      verifyBamHomChk.pl \
      -o ${outputDir + "/contamination"} \
      -b ${bamFile} \
      -d ${contamDownsampOneIn} \
      -j ${outputDir + "/contamination/" + baseName + "_" + process + "_summary.json"}
      -a ${ascatSegmentFiles}
    else
      verifyBamHomChk.pl \
      -o ${outputDir + "/contamination"} \
      -b ${bamFile} \
      -d ${contamDownsampOneIn} \
      -j ${outputDir + "/contamination/" + baseName + "_" + process + "_summary.json"}
    fi
  >>>

  output {
    File contamFile = "${outputDir}/contamination/${baseName}_${process}_summary.json"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task bam_stats {
  File bamFile
  String baseName
  String outputDir

  command {
    bam_stats -i ${bamFile} \
              -o ${outputDir + "/" + baseName + ".bas"}
  }

  output {
    File basFile = "${outputDir}/${baseName}.bas"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task bbAlleleCount {
  File bamFile
  File bbRefLoci
  String baseName
  String outputDir

  command {
    execute_with_sample ${bamFile} alleleCounter \
    -l ${bbRefLoci} \ 
    -o ${outputDir + "/bbCounts/" + baseName + "_" + bbRefLoci + ".tsv"} \
    -b ${bamFile};
  }

  output {
    Array[File] alleleCounts = glob("${outputDir}/${baseName}.*.tsv")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task qc_metrics {
  File controlBam
  File tumorBam
  String outputDir

  command {
    qc_and_metrics.pl ${outputDir} ${controlBam} ${tumorBam}
  }

  output {
    File qc_metrics = "${outputDir}/qc_metrics.json"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task ascat {
  File tumorBam
  File controlBam
  File genomeFa
  File genomeFai
  File snpPosFile
  File snpLociFile
  File snpGcCorrectionsFile
  String process
  Int index
  String seqType
  String assembly
  String species
  String gender
  String outputDir

  command {
    ascat.pl
    -p ${process} \
    -i ${index} \
    -r ${genomeFa} \
    -pr ${seqType} \
    -ra ${assembly}
    -rs ${species} \
    -g ${gender} \
    -pl "ILLUMINA" \    
    -s ${snpLociFile} \
    -sp ${snpPosFile} \
    -sg ${snpGcCorrectionsFile} \
    -o ${outputDir + "/ascat"} \
    -t ${tumorBam} \
    -n ${controlBam} \
    -f
  }

  output {
    Array[File] ascatOutput = glob("${outputDir}/ascat/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task pindel {
  File tumorBam
  File controlBam
  File genomeFa
  File genomeFai
  File simpleRepeatsFile
  File vcfFilterRulesFile
  File vcfFilterSoftRulesFile
  File codingGeneFootprintsFile
  File unmatchedNormalPanelGff3
  File badAnchorLociFile
  String refExclude = "MT,GL%,hs37d5,NC_007605"
  String seqType
  String assembly
  String species
  String outputDir
  String process
  Int? pindelInputThreads
  Int? pindelNormalisedThreads
  Int? index

  command {
      pindel.pl \
      -p ${process} \
      -r ${genomeFa} \
      -e ${refExclude} \
      -st ${seqType} \
      -as ${assembly}
      -sp ${species} \
      -s ${simpleRepeatsFile} \
      -f ${vcfFilterRulesFile} \
      -g ${codingGeneFootprintsFile} \
      -u ${unmatchedNormalPanelGff3} \
      -sf ${vcfFilterSoftRulesFile} \
      -b ${badAnchorLociFile} \
      -o ${outputDir + "/pindel"} \
      -t ${tumorBam} \
      -n ${controlBam}
      ${"-i " + index}
      ${"-c " + pindelInputThreads}
      ${"-l " + pindelNormalisedThreads}
  }

  output {
    Array[File] pindelOutput = glob("${outputDir}/pindel/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


task brass {
  File tumorBam
  File controlBam
  File genomeFa
  File genomeFai
  File refExclude
  File ignoredRegionsFile
  File normalPanelGroupsFile
  File genomeCacheFile
  File virusSeqsFile
  String microbeSeqsFilesPrefix
  Array[File] microbeSeqsFiles
  File bedCoordFile
  File? cnPath
  File? cnStats
  String process
  String refBase
  String seqType
  String assembly
  String species
  String outputDir
  Int? index
  Int? threads  

  command {
    brass.pl \
    -j 4 \
    -k 4 \
    -p ${process} \
    -g ${genomeFa} \
    -e ${refExclude} \
    -pr ${seqType} \
    -as ${assembly} \
    -s ${species} \
    -pl "ILLUMINA" \
    -d  ${ignoredRegionsFile} \
    -f  ${normalPanelGroupsFile} \
    -g_cache  ${genomeCacheFile} \
    -o ${outputDir + "/brass"} \
    -t ${tumorBam} \
    -n ${controlBam} \
    -vi ${virusSeqsFile} \
    -mi ${microbeSeqsFilesPrefix} \
    -b ${bedCoordFile} \
    ${"-i " + index} \
    ${"-a" + cnPath} \
    ${"-ss" + cnStats} \
    ${"-c " + threads} \
    ${"-l " + threads}
  }

  output {
    Array[File] brassOut = glob("${outputDir}/brass/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


task caveCnPrep {
  File cnPath
  String type
  String outputDir
  
  command <<<
  if [[${type} -e "tumor"]]; then 
    export OFFSET=6 ;
  else
    export OFFSET=4 ;
  fi ;
  perl -ne '@F=(split q{,}, $_)[1,2,3," + $OFFSET + "]; $F[1]-1; print join(\"\\t\",@F).\"\\n\";' < ${cnPath} > ${outputDir + "/" + type + ".cn.bed"} ;
  >>>
 
  output {
    File caveCnPrepOut = "${outputDir}/${type}.cn.bed"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


task caveman {
  File tumorBam
  File controlBam
  File genomeFa
  File genomeFai
  File ascatContamFile
  File ignoredRegionsFile
  File? pindelGermlineMutsFile
  String flagBedFilesDir
  Array[File] flagBedFiles
  String unmatchedNormalFilesDir
  Array[File] unmatchedNormalFiles
  String process
  String seqType
  String assembly
  String species
  String seqProtocol
  String outputDir
  Int? index
  Int? threads
  
  command {
    caveman.pl \
    -p ${process} \
    -ig ${ignoredRegionsFile} \
    -b ${flagBedFilesDir} \
    -np ${seqType} \
    -tp ${seqType} \
    -sa ${assembly} \
    -s ${species} \
    -st ${seqProtocol} \
    -o ${outputDir + "/caveman"} \
    -tc ${outputDir + "/tumor.cn.bed"} \
    -nc ${outputDir + "/normal.cn.bed"} \
    -k ${ascatContamFile} \
    -tb ${tumorBam} \
    -nb ${controlBam} \
    -r ${genomeFai} \
    -u ${unmatchedNormalFilesDir} \
    ${"-i " + index} \
    ${"-in " + pindelGermlineMutsFile} \
    ${"-l " + threads} \
    ${"-t " + threads}
  }

  output {
    Array[File] cavemanOut = glob("${outputDir}/caveman/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

workflow sanger_cgp_somatic_vc {
  File controlBam
  File tumorBam
  File genomeFa
  File genomeFai
  String seqType = "WGS"
  String seqProtocol = "genomic"
  String assembly = "GRCh37"
  String species = "human"
  String gender
  String globalOutputDir = "/output/"

  # bbAlleleCount
  Array[File] bbRefLociFiles

  # ASCAT
  File snpPosFile
  File snpLociFile
  File snpGcCorrectionsFile

  # PINDEL
  File simpleRepeatsFile
  File vcfFilterRulesFile
  File vcfFilterSoftRulesFile
  File codingGeneFootprintsFile
  File unmatchedNormalPanelGff3
  File badAnchorLociFile
  Int pindelInputThreads
  Int pindelNormalisedThreads
  Array[Int] pindelScatterIndices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

  # BRASS
  File refExclude
  File ignoredRegionsFile
  File normalPanelGroupsFile
  File genomeCacheFile
  File virusSeqsFile
  String microbeSeqsFilesPrefix
  Array[File] microbeSeqsFiles
  File bedCoordFile
  File cnPath
  File cnStats
  Int brassThreads

  # CAVEMAN
  File ascatContamFile
  File ignoredRegionsFile
  File? pindelGermlineMutsFile
  String flagBedFilesDir
  Array[File] flagBedFiles
  String unmatchedNormalFilesDir
  Array[File] unmatchedNormalFiles
  Array[Int] cavemandSplitIndices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86]
  Int cavemanMstepThreads
  Int cavemanEstepThreads

  ##
  ## QC/Prep steps
  ##
  call getSampleId as tumor_sampleId {
    input: inBam = tumorBam
  }

  call getSampleId as control_sampleId {
    input: inBam = controlBam
  }

  call compareGenotype {
    input: controlBam = controlBam, 
            tumorBam = tumorBam,
            baseName = "genotype",
            outputDir = globalOutputDir
  }

  call analyzeContamination as control_contam {
    input: process = "control",
            bamFile = controlBam,
            baseName = control_sampleId.SM,
            outputDir = globalOutputDir
  }

  call analyzeContamination as tumor_contam {
    input: process = "tumor",
            bamFile = tumorBam,
            baseName = tumor_sampleId.SM,
            outputDir = globalOutputDir
  }
  
  call bam_stats as control_bam_stats {
    input: bamFile = controlBam, 
            baseName = control_sampleId.SM,
            outputDir = globalOutputDir
  }

  call bam_stats as tumor_bam_stats {
    input: bamFile = tumorBam, 
            baseName = tumor_sampleId.SM,
            outputDir = globalOutputDir
  }
  
  scatter(chrLoci in bbRefLociFiles) {
    call bbAlleleCount as control_bbAlleleCount {
      input: bamFile = controlBam, 
              baseName = control_sampleId.SM,
              bbRefLoci = chrLoci,
              outputDir = globalOutputDir
    }
  
    call bbAlleleCount as tumor_bbAlleleCount {
      input: bamFile = tumorBam, 
              baseName = tumor_sampleId.SM,
              bbRefLoci = chrLoci,
              outputDir = globalOutputDir
    }
  }

  call qc_metrics {
    input: controlBam = controlBam,
            tumorBam = tumorBam,
            outputDir = globalOutputDir
  }

  ##
  ## ASCAT - copynumber analysis
  ##
  call ascat as ascat_allele_count {
    input: process = "allele_count",
            controlBam = controlBam,
            tumorBam = tumorBam,
            index = 1, 
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            snpPosFile = snpPosFile,
            snpLociFile = snpLociFile,
            snpGcCorrectionsFile = snpGcCorrectionsFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            gender = gender,
            outputDir = globalOutputDir
  }
  
  call ascat {
    input: process = "ascat",
            controlBam = controlBam,
            tumorBam = tumorBam,
            index = 1, 
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            snpPosFile = snpPosFile,
            snpLociFile = snpLociFile,
            snpGcCorrectionsFile = snpGcCorrectionsFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            gender = gender,
            outputDir = globalOutputDir
  }

  call ascat as ascat_finalise {
    input: process = "finalise",
            controlBam = controlBam,
            tumorBam = tumorBam,
            index = 1, 
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            snpPosFile = snpPosFile,
            snpLociFile = snpLociFile,
            snpGcCorrectionsFile = snpGcCorrectionsFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            gender = gender,
            outputDir = globalOutputDir
  }


  ##
  ## Pindel - InDel calling
  ##
  call pindel as pindel_input1 {
    input: process="input", 
            pindelInputThreads = pindelInputThreads,
            index = 1,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            simpleRepeatsFile = simpleRepeatsFile,
            vcfFilterRulesFile = vcfFilterRulesFile,
            vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
            codingGeneFootprintsFile = codingGeneFootprintsFile,
            unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
            badAnchorLociFile = badAnchorLociFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call pindel as pindel_input2 {
    input: process="input", 
            pindelInputThreads = pindelInputThreads,
            index = 2,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            simpleRepeatsFile = simpleRepeatsFile,
            vcfFilterRulesFile = vcfFilterRulesFile,
            vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
            codingGeneFootprintsFile = codingGeneFootprintsFile,
            unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
            badAnchorLociFile = badAnchorLociFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call pindel {
    input: process="pindel",
            pindelInputThreads = pindelNormalisedThreads,
            pindelNormalisedThreads = pindelNormalisedThreads,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            simpleRepeatsFile = simpleRepeatsFile,
            vcfFilterRulesFile = vcfFilterRulesFile,
            vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
            codingGeneFootprintsFile = codingGeneFootprintsFile,
            unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
            badAnchorLociFile = badAnchorLociFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  # number of refs to process
  # 1-22, X, Y
  scatter(i in pindelScatterIndices) {
    call pindel as pindel_pin2vcf {
      input: process = "pin2vcf", 
              index = i,
              controlBam = controlBam,
              tumorBam = tumorBam,
              genomeFa = genomeFa,
              genomeFai = genomeFai, 
              simpleRepeatsFile = simpleRepeatsFile,
              vcfFilterRulesFile = vcfFilterRulesFile,
              vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
              codingGeneFootprintsFile = codingGeneFootprintsFile,
              unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
              badAnchorLociFile = badAnchorLociFile,
              seqType = seqType,
              assembly = assembly,
              species = species,
              outputDir = globalOutputDir
    }
  }

  call pindel as pindel_merge {
    input: process = "merge",
            index = 1, 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            simpleRepeatsFile = simpleRepeatsFile,
            vcfFilterRulesFile = vcfFilterRulesFile,
            vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
            codingGeneFootprintsFile = codingGeneFootprintsFile,
            unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
            badAnchorLociFile = badAnchorLociFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call pindel as pindel_flag {
    input: process = "flag", 
            index = 1,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            simpleRepeatsFile = simpleRepeatsFile,
            vcfFilterRulesFile = vcfFilterRulesFile,
            vcfFilterSoftRulesFile = vcfFilterSoftRulesFile,
            codingGeneFootprintsFile = codingGeneFootprintsFile,
            unmatchedNormalPanelGff3 = unmatchedNormalPanelGff3,
            badAnchorLociFile = badAnchorLociFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  ##
  ## BRASS - breakpoint analysis
  ##
  call brass as brass_input1 {
    input: process = "input",
            index = 1,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_input2 {
    input: process = "input", 
            index = 2,
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_cover {
    input: process = "cover", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_merge {
    input: process = "merge", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_group {
    input: process = "group", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_isize {
    input: process = "isize", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_normcn {
    input: process = "normcn", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            cnPath = cnPath,
            cnStats = cnStats,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_filter {
    input: process = "filter", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            cnPath = cnPath,
            cnStats = cnStats,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_split {
    input: process = "split", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            refBase = referenceDir, 
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            outputDir = globalOutputDir
  }

  call brass as brass_assemble {
    input: process = "assemble", 
            controlBam = controlBam,
            tumorBam = tumorBam,
            genomeFa = genomeFa,
            genomeFai = genomeFai, 
            refExclude = refExclude,
            ignoredRegionsFile = ignoredRegionsFile,
            normalPanelGroupsFile = normalPanelGroupsFile,
            genomeCacheFile = genomeCacheFile,
            virusSeqsFile = virusSeqsFile,
            microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
            microbeSeqsFiles = microbeSeqsFiles,
            bedCoordFile = bedCoordFile,
            seqType = seqType,
            assembly = assembly,
            species = species,
            threads = brassThreads,
            outputDir = globalOutputDir
  }

  call brass as brass_grass {
    input: process = "grass", 
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai, 
           refExclude = refExclude,
           ignoredRegionsFile = ignoredRegionsFile,
           normalPanelGroupsFile = normalPanelGroupsFile,
           genomeCacheFile = genomeCacheFile,
           virusSeqsFile = virusSeqsFile,
           microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
           microbeSeqsFiles = microbeSeqsFiles,
           bedCoordFile = bedCoordFile,
           cnPath = cnPath,
           cnStats = cnStats,
           seqType = seqType,
           assembly = assembly,
           species = species,
           outputDir = globalOutputDir
  }

  call brass as brass_tabix {
    input: process = "tabix", 
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai, 
           refExclude = refExclude,
           ignoredRegionsFile = ignoredRegionsFile,
           normalPanelGroupsFile = normalPanelGroupsFile,
           genomeCacheFile = genomeCacheFile,
           virusSeqsFile = virusSeqsFile,
           microbeSeqsFilesPrefix = microbeSeqsFilesPrefix,
           microbeSeqsFiles = microbeSeqsFiles,
           bedCoordFile = bedCoordFile,
           cnPath = cnPath,
           cnStats = cnStats,
           seqType = seqType,
           assembly = assembly,
           species = species,
           outputDir = globalOutputDir
  }

  ##
  ## Caveman - SNV analysis
  ##
  # TODO - caveCnPrep
  call caveCnPrep as control_caveCnPrep {
    input: type = "control",
           CnPath = ,
           outputDir = globalOutputDir
  }

  call caveCnPrep as tumor_caveCnPrep {
    input: type = "tumor",
           CnPath = ,
           outputDir = globalOutputDir
  }

  scatter(i in caveManSplitIndices) {
    call caveman as caveman_split {
      input: process = "split", 
             index = i,
             controlBam = controlBam,
             tumorBam = tumorBam,
             genomeFa = genomeFa,
             genomeFai = genomeFai,
             ascatContamFile = ascatContamFile,
             ignoredRegionsFile = ignoredRegionsFile,
             pindelGermlineMutsFile = pindelGermlineMutsFile,
             flagBedFilesDir = flagBedFilesDir,
             flagBedFiles = flagBedFiles,
             unmatchedNormalFilesDir = unmatchedNormalFilesDir,
             unmatchedNormalFiles = unmatchedNormalFiles,
             seqType = seqType,
             assembly = assembly,
             species = species,
             seqProtocol = seqProtocol,
             outputDir = globalOutputDir
    }
  }

  call caveman as caveman_split_concat {
    input: process = "split_concat", 
           index = 1,
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
  }

  call caveman as caveman_mstep {
    input: process = "mstep", 
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
           threads = cavemanMstepThreads
  }

  call caveman as caveman_merge {
    input: process = "merge", 
           index = 1,
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
  }

  call caveman as caveman_estep {
    input: process = "estep",
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
           threads = cavemanEstepThreads
  }

  call caveman as caveman_merge_results {
    input: process = "merge_results", 
           index = 1,
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
  }

  call caveman as caveman_add_ids {
    input: process = "add_ids", 
           index = 1,
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
  }

  call caveman as caveman_flag {
    input: process = "flag", 
           index = 1,
           controlBam = controlBam,
           tumorBam = tumorBam,
           genomeFa = genomeFa,
           genomeFai = genomeFai,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           ascatContamFile = ascatContamFile,
           ignoredRegionsFile = ignoredRegionsFile,
           pindelGermlineMutsFile = pindelGermlineMutsFile,
           flagBedFilesDir = flagBedFilesDir,
           flagBedFiles = flagBedFiles,
           unmatchedNormalFilesDir = unmatchedNormalFilesDir,
           unmatchedNormalFiles = unmatchedNormalFiles,
           seqType = seqType,
           assembly = assembly,
           species = species,
           seqProtocol = seqProtocol,
           outputDir = globalOutputDir
  }
}
