task get_bam_basename {
  File bamFile

  command {
    basename ${bamFile} .bam
    
  }

  output {
    String base = read_string(stdout())
  }

  runtime {
    docker: "bwa-workflow"
  }
}

task bbAlleleCount {
  File bamFile
  String bamName
  String refBase
  String outputDir

  command <<<
    for chr in {1..23}; do
      execute_with_sample ${bamFile} alleleCounter \
      -l ${refBase + "/battenberg/1000genomesloci/1000genomesloci2012_chr" + chr + ".txt"} \ 
      -o ${outputDir + "/" + bamName + "." + chr ".tsv"} \
      -b ${bamFile};
    done
  >>>

  output {
    Array[File] alleleCounts = glob("${outputDir}/${bamName}.*.tsv")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

# TODO
task bbAlleleMerge {
  File controlBam
  String bbDir

  command {
    packageImpute.pl ${controlBam} ${bbDir}
  }

  output {
    # File alleleCounts = "${outputDir}/${bamName}.${chr}.tsv"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task bam_stats {
  File bamFile
  String bamName
  String outputDir

  command {
    bam_stats -i ${bamFile} \
              -o ${outputDir + "/" + bamName + ".bas"}
  }

  output {
    File bamStats = ${outputDir + "/" + bamName + ".bas"}
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task qc_and_metrics {
  File controlBamFile
  Array[File] tumorBamFiles
  String outputDir

  command {
    qc_and_metrics.pl ${outputDir} ${controlBamFile} ${sep=" " tumorBamFiles}
  }

  output {
    File qc_metrics = "${outputDir}/qc_metrics.json"
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


task compareGenotype {
  File normalBamFile
  Array[File]+ tumorBamFiles
  String bamName
  String outputDir

  command {
    compareBamGenotypes.pl \
    -o ${outputDir + "/genotype"} \
    -nb ${normalBamFile} \
    -j ${outputDir + "/genotype/" + bamName + "_summary.json"} \
    -tb ${sep=" -tb " tumorBamFiles}
  }

  output {
    File genotype = ${outputDir + "/genotype/" + bamName + "_summary.json"}
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

# TODO
task analyzeContamination {
  File bamFile
  String bamName
  String outputDir
  # File contamDownsampOneIn

  command {
    verifyBamHomChk.pl \
    -o ${outputDir + "/contamination"} \
    -b ${bamFile} \
    -d ${contamDownsampOneIn} \
    -j ${outputDir + "/contamination/" + bamName + "_summary.json"} 

    # if(process.equals("tumour")) {
    #   thisJob.getCommand().addArgument("-a " + OUTDIR + "/" + tumourCount + "/ascat/*.copynumber.caveman.csv"); // not the best approach but works
    # }

  }

  output {
    File contamFile = ${outputDir + "/contamination/" + bamName + "_summary.json"}
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task ascat {
  File tumorBamFile
  File normalBamFile
  File genomeFa
  String refBase
  String seqType
  String assembly
  String species
  String outputDir
  String process
  String gender
  Int index
  Int tumorCount

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
    -s ${refBase + "/ascat/SnpLocus.tsv"} \
    -sp ${refBase + "/ascat/SnpPositions.tsv"} \
    -sg ${refBase + "/ascat/SnpGcCorrections.tsv"} \
    -o ${outputDir + "/" + tumorCount + "/ascat"} \
    -t ${tumorBamFile} \
    -n ${normalBamFile} \
    -f
}

  output {
    Array[File] ascatOutput = glob("${outputDir}/${tumorCount}/ascat/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }

}

# TODO
task pindel {
  File tumorBamFile
  File normalBamFile
  File genomeFa
  File refExclude
  String refBase
  String seqType
  String assembly
  String species
  String outputDir
  String process
  Int pindelInputThreads
  Int index

  command {
    pindel.pl \
    -p ${process} \
    -r ${genomeFa} \
    -e ${refExclude} \
    -st ${seqType} \
    -as ${assembly}
    -sp ${species} \
    -s ${refBase + "/pindel/simpleRepeats.bed.gz"} \
    -f ${refBase + "/pindel/genomicRules.lst"} \
    -g ${refBase + "/pindel/human.GRCh37.indelCoding.bed.gz"} \
    -u ${refBase + "/pindel/pindel_np.gff3.gz"} \
    -sf ${refBase + "/pindel/softRules.lst"} \
    -b ${refBase + "/brass/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz"} \
    -o ${outputDir + "/pindel"} \
    -t ${tumorBamFile} \
    -n ${normalBamFile} \
    -c ${pindelInputThreads}

    # if(!process.equals("pindel")) {
    #   thisJob.getCommand().addArgument("-i " + index);
    # }

  }

  output {
    Array[File] pindelOutput = glob("${outputDir}/pindel/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task caveCnPrep {
  File cnPath
  Int tumorCount
  String type
  String outputDir
  
  command <<<
  if [[${type} -e "tumor"]]; then 
    export OFFSET=6 ;
  else
    export OFFSET=4 ;
  fi ;
  perl -ne '@F=(split q{,}, $_)[1,2,3," + $OFFSET + "]; $F[1]-1; print join(\"\\t\",@F).\"\\n\";' < ${CnPath} > ${outputDir + "/" + tumorCount + "/" + type + ".cn.bed"} ;
  >>>
 
  output {
    File caveCnPrepOut = ${outputDir + "/" + tumorCount + "/" + type + ".cn.bed"}
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

# TODO
task cavemanBaseJob {
  String process
  String refBase
  String seqType
  String assembly
  String species
  String seqProtocol
  Int tumorCount
  File ascatContamFile
  File tumorBam
  File controlBam
  File genomeFa
  File genomeFai
  String outputDir
  
  command {
    # process logic
    # if [[]] ; then
    # ;
    # else
    # fi ;
    caveman.pl \
    -p ${process} \
    -ig ${refBase + "/caveman/ucscHiDepth_0.01_merge1000_no_exon.tsv"} \
    -b ${refBase + "/caveman/flagging"} \
    -np ${seqType} \
    -tp ${seqType} \
    -sa ${assembly} \
    -s ${species} \
    -st ${seqProtocol} \
    -o ${outputDir + "/" + tumourCount + "/caveman"} \
    -tc ${outputDir + "/" + tumourCount + "/tumour.cn.bed"} \
    -nc ${outputDir + "/" + tumourCount + "/normal.cn.bed"} \
    -k ${ascatContamFile} \
    -tb ${tumorBam} \
    -nb ${normalBam} \
    -r ${genomeFai} \
    -u ${refBase + "/caveman"}
    # ${processFlag}
  }

  output {
    # File cavemanOut = ${outputDir + "/"}
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}

task brass {
  File tumorBam
  File controlBam
  File genomeFa
  File refExclude
  String process
  String refBase
  String seqType
  String assembly
  String species
  String outputDir
  Int tumorCount

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
    -d  ${refBase + "/brass/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz"} \
    -f  ${refBase + "/brass/brass_np.groups.gz"} \
    -g_cache  ${refBase + "/vagrent/e75/Homo_sapiens.GRCh37.75.vagrent.cache.gz"} \
    -o ${outputDir + "/" + tumourCount + "/brass"} \
    -t ${tumourBam} \
    -n ${controlBam} \
    -vi ${refBase + "/brass/viral.1.1.genomic.fa"} \
    -mi ${refBase + "/brass/all_ncbi_bacteria.20150703"} \
    -b ${refBase + "/brass/hs37d5_500bp_windows.gc.bed.gz"}
  }
  
  output {
    Array[File] brassOut = glob("${outputDir}/brass/*")
  }

  runtime {
    docker: "sanger-somatic-vc-workflow"
  }
}


workflow sanger-cgp-somatic-vc {
  Array[Int] cavemandSplitIndices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86]

  ##
  ## QC/Prep steps
  ##
  call compareGenotype {
    inputs:
  }

  call analyzeContamination as control_contam {
    inputs:
  }

  call analyzeContamination as tumor_contam {
    inputs:
  }
  
  call bam_stats as control_bam_stats {
    inputs:
  }

  scatter(tbam in tumorBams) {
    call bam_stats {
      inputs:
    }
    
    call bbAlleleCount {
      inputs:
    }
  }

  call bbAlleleMerge {
   inputs:
  }


  ##
  ## ASCAT - copynumber analysis
  ##
  # scatter ?
  call ascat as ascat_allele_count {
    inputs: process=allele_count,
  }
  
   call ascat {
    inputs: process=ascat,
  }

  call ascat as ascat_finalise {
    inputs: process=finalise,
  }


  ##
  ## Pindel - InDel calling
  ##
  call pindel as pindel_input1 {
    inputs: process=input, 
            index=1,
  }

  call pindel as pindel_input2 {
    inputs: process=input, 
            index=2,
  }

  call pindel {
    inputs: process=pindel, 
  }

  # scatter?
  call pindel as pindel_pin2vcf {
    inputs: process=pin2vcf, 
  }

  call pindel as pindel_merge {
    inputs: process=merge, 
  }

  call pindel as pindel_flag {
    inputs: process=flag, 
  }


  ##
  ## BRASS - breakpoint analysis
  ##
  call brass as brass_input1 {
    inputs: process=input, 
  }

  call brass as brass_input2 {
    inputs: process=input, 
            index=2,
  }

  call brass as brass_cover {
    inputs: process=cover, 
  }

  call brass as brass_merge {
    inputs: process=merge, 
  }

  call brass as brass_group {
    inputs: process=group, 
  }

  call brass as brass_isize {
    inputs: process=isize, 
  }

  call brass as brass_normcn {
    inputs: process=normcn, 
  }

  call brass as brass_filter {
    inputs: process=filter, 
  }

  call brass as brass_split {
    inputs: process=split, 
  }

  call brass as brass_assemble {
    inputs: process=split, 
  }

  call brass as brass_grass {
    inputs: process=grass, 
  }

  call brass as brass_tabix {
    inputs: process=tabix, 
  }

  ##
  ## Caveman - SNV analysis
  ##
  call caveCnPrep as control_caveCnPrep {
    inputs: type=control,
  }

  call caveCnPrep as tumor_caveCnPrep {
    inputs: type=tumor, 
  }

  scatter(i in caveManSplitIndices) {
    call caveman as caveman_split {
      inputs: process=split, 
              index=i,
    }
  }

  call caveman as caveman_split_concat {
    inputs: process=split_concat, 
  }

  call caveman as caveman_mstep {
    inputs: process=mstep, 
  }

  call caveman as caveman_merge {
    inputs: process=merge, 
  }

  call caveman as caveman_estep {
    inputs: process=estep, 
  }

  call caveman as caveman_merge_results {
    inputs: process=merge_results, 
  }

  call caveman as caveman_add_ids {
    inputs: process=add_ids, 
  }

  call caveman as caveman_flag {
    inputs: process=flag, 
  }

}
