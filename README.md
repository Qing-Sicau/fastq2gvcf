# fastq2gvcf
from fastq file to g.vcf file using GATK pipeline
make sure make sure bwa samblaster gatk java samtools R Rtools ggplot2 gplots picard seqkit are all properly installed.
change the dataDir, config file, refDir, ref, and step from which you begin, eg. from step 4 when you already have bwa which duplications were tagged.
need to paralell
## BEFORE RUN THE PIPELINE, RUN THIS:
```
samtools faidx $ref
picard CreateSequenceDictionary \
R="$ref" \
O="$refDir/Ppse.genome.dict"
```
once done, the pipeline can run properly.

The pipeline was modified refferring to:
https://github.com/josieparis/gatk-snp-calling
https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
