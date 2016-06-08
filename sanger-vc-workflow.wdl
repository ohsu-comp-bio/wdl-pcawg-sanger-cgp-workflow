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
  File bbChrRefFile
  String outputDir

  command {
    execute_with_sample ${bamFile} alleleCounter \
    -l ${bbChrRefFile} \ 
    -o ${outputDir + "/" + bamName + "." + chr ".tsv"} \
    -b ${bamFile}
  }

  output {
    File alleleCounts = ${outputDir + "/" + bamName + "." + chr ".tsv"}
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

}
