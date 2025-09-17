# The QC Navigator: Bacterial WGS Edition

By Boas van der Putten, Jennifer L. Guthrie, Nora Berg Blomseth, Ewelina Kamińska

## Scope

This document focuses on bioinformatics quality control (QC) for **bacterial isolate** whole-genome sequencing (WGS) using Illumina technology. The guidance applies to three main applications commonly used in clinical and public health microbiology, including SNP-based outbreak detection, genotyping, and antimicrobial-resistance identification/characterization.

Other sequencing contexts (e.g., metagenomics, viral or fungal genomics) are outside the scope of this QC framework.

## Contents


[Quality Control Metrics Overview](#quality-control-metrics-overview)

[Software to Perform Quality Control](#Software-to-Perform-Quality-Control)

[Selecting the Best Assembly](#Selecting-the-Best-Assembly)

[Organism-Specific Considerations](#organism-specific-considerations)

[Table 1. Table 1. Key Quality Control Metrics](#metrics-tbl)

[Table 2. Bioinformatics QC Tools, Metrics, and Purposes.](#tools-tbl)

[Table 3. Example-Sequences with Quality Control Metrics](#example-sequences)


## Quality Control Metrics Overview

Quality control metrics provide a structured way to assess the reliability, completeness, and accuracy of bacterial whole-genome sequencing data. Evaluating these metrics ensures that downstream analyses are based on high-quality data and reduces the risk of incorrect or misleading results.

Several QC metrics are correlated with each other, e.g. if a sample has a larger genome size than expected for that species, that is often due to contamination. Often samples fail multiple metrics at the same time, although it is considered good practice to check different metrics as in some circumstances a sample might fail only a single metric.

### <a name="metrics-tbl"></a>**Table 1. Key Quality Control Metrics for Bacterial Whole-Genome Sequencing.**

| **Metric name** | **Metric description** | **Recommended Thresholds** | **What is wrong if the threshold is not met?** | **Metric Type** | **Tools (\*see Table 2)** |
| --- | --- | --- | --- | --- | --- |
| N50 | Length of the shortest contig in the sequenced genome that together with other (longer) contigs make up at least 50% of the entire genome | Qualibact<sup>a</sup>/Aquamis<sup>b</sup> | Assembly may be too fragmented for reliable analysis | Assembly quality | 4 |
| Total assembly length | Sum of the lengths of all contigs or scaffolds in an assembled genome. Can be compared to the expected genome size to assess completeness | Qualibact<sup>a</sup>/Aquamis<sup>b</sup> | Too short may indicate missing data; too long may suggest contamination or misassembly | Assembly quality | 4 |
| GC Content | Percentage of the bases guanine and cytosine in the genome. The GC content differs between species | Qualibact<sup>a</sup>/Aquamis<sup>b</sup> | Suggests contamination with another species | Sequence quality / Assembly | 1, 3, 4 |
| Number of contigs | Count of contiguous assembled sequences in a genome assembly | Qualibact<sup>a</sup>/Aquamis<sup>b</sup>b | Fewer contigs generally indicate a more contiguous and higher-quality assembly, while many contigs may suggest fragmentation or assembly difficulties | Sequence quality / Assembly | 4 |
| Depth (/depth of coverage) | Average number of times each base in the genome is sequenced | Average 50–100× and >95% of the genome should meet a minimum depth (10× for Illumina, preference for 30×) to ensure reliable assembly and variant calling | Base calls may be unreliable and variants may be missed | Sequence quality | 5 |
| Coverage (/breadth of coverage) | Proportion of the reference genome (or expected genome size) that is represented by sequencing reads | Aim for ≥95% of the reference genome covered by high-quality reads | Parts of the genome may be missing from the assembly | Sequence quality | 5 |
| Completeness | Estimated completeness of genome as determined from the presence/absence of marker genes | Qualibact<sup>a</sup> | Parts of the genome may be missing from the assembly |  |  |
| Contamination | Estimated contamination of genome e.g. determined by the presence of multi-copy marker genes or by sequences assigned to species | Qualibact<sup>a</sup> (CheckM)/ Aquamis<sup>b</sup> (ConFindr) | Foreign or unexpected DNA may be present, which can compromise the accuracy of the assembly and downstream analyses | Sequence quality | 1, 3, 4 |
| Total coding sequences | Number of open reading frames identified in a genome | Qualibact<sup>a</sup> | Too low suggests part of the genome is missing, too high might suggest contamination | Assembly / Annotation metric | 4 |
| Strand bias | An imbalance in sequencing reads supporting a variant between the forward and reverse strands | Variants should be supported by reads from both forward and reverse strands, with a minimal number of reads per strand (typically ≥2–5); strong imbalance (e.g., >90% of reads on a single strand) should be filtered | Variant may be an artifact rather than a true biological signal | Sequence quality / Variant-calling | 6 |
| Read quality score (Phred score) | A measure of the accuracy of base calling in sequencing reads, where higher scores indicate lower error probabilities | Phred score (Q)>30 equals an accurency of more than 99.9% (Q20=99% accuracy) / minimum percent of bases with quality of Q30 (MiSeq): >65%, (NextSeq): ≥ 75 % [Guidance document for WGS-Benchmarking] | Sequencing errors may compromise downstream analyses | Sequence quality | 1 |
| Read length | Length of raw or trimmed FASTQ reads | Dependent on experimental setup. Mean length of Illumina reads should stay above 90% of intended read length. | Too low suggests sequencing experiment was not optimal (e.g. degraded DNA; library prep issues), or excessive trimming | Sequence quality | 4 |
| Number of reads mapping back to assembly | What proportion of sequencing reads can be mapped back to the de novo assembly made from these reads | Typically ≥95% of reads should map back to the assembly for a reliable genome | A large proportion of sequencing reads may not align to the assembly, suggesting incomplete or poor-quality assembly, contamination, or sequencing errors | Assembly metric | 5 |
| Coding sequence length | Whether coding sequences have the expected length | CDS lengths should generally fall within the expected range for the species (e.g., bacterial CDSs often 100–4,000 bp); extreme deviations should be investigated | Aberrant coding sequence length suggests presence of insertions/deletions (indels) compared to true sequence or missasembly or partial genes | Assembly / Annotation metric | 4 |


## Common Issues Addressed by Quality Control

There are a number of underlying issues that can be detected through quality control. The list below describes common issues and how these can be detected.

* **Interspecies Contamination** – DNA from other bacterial species can lead to misassemblies, inflated genome size, false gene predictions, and incorrect variant calls.
  + Detect by:
    - Taxonomic classification of reads or contigs using Kraken, Centrifuge, or Kaiju
    - Duplicated single-copy marker genes or unexpected contigs (CheckM, QUAST)
    - Deviations in GC content or abnormal contig characteristics
* **Intraspecies Contamination** – Presence of multiple strains of the same species in one sample can cause chimeric assemblies, inflated SNP counts, ambiguous variant calls, and difficulty resolving population structure.
  + Detect by:
    - Mixed allele frequencies at variant sites (e.g., >1 base consistently present at the same position)
    - Unusually high heterozygosity in bacterial WGS (normally clonal)
    - Multiple copies of “single-copy” genes (CheckM, assembly annotation)
    - Read mapping showing inconsistent coverage or mismatches
* **Homopolymer(s) error** – Sequencing errors occur in stretches of identical bases (e.g., AAAAA), leading to incorrect base counts, frameshifts, or false indels in assemblies and variant calls. More common in long-read platforms
  + Detect by:
    - Elevated indel rates in or near homopolymer regions
    - Discrepancies between sequencing technologies (e.g., Illumina vs. ONT)
    - Error patterns flagged by variant callers (low confidence indels)
* **Host contamination (e.g., human; cell line)** – Presence of host DNA reduces sequencing efficiency for the target bacteria, lowers coverage, and can confound assembly or downstream analyses.
  + Detect by:
    - taxonomic classification of reads using Kraken
    - Check for unexpected GC content or contigs in the assembly (QUAST, CheckM)
    - Look for a low proportion of reads mapping to the bacterial reference genome
* **Sequencer carryover** – Residual DNA from a previous sequencing run contaminates the current run, leading to false detection of organisms or spurious reads in the dataset.
  + Detect by:
    - Low-level presence of reads from unrelated organisms (detected via Kraken)
    - Detection of sequences from samples run in previous batches
    - Unexpected coverage of control/reference sequences not included in the run
* **“Kit-ome”** – Sequencing reagents contaminated
  + Detect by
    - Taxonomic classification of reads (Kraken) shows consistent low-level hits to “odd” organisms (e.g., Ralstonia, Acinetobacter, Pseudomonas) across negative controls or unrelated projects.
  + Resolve by contacting company, also sequencing to only 100X to avoid “oversampling” contaminants
* **Low yield** – Insufficient DNA quantity can lead to poor library preparation, low sequencing depth, uneven coverage, and incomplete assemblies.
  + Detect by looking for low read counts, low coverage/depth, or a low proportion of reads mapping back to assembly
* **Low base calling quality** – Inaccurate base calls can lead to sequencing errors, misassemblies, and false-positive or false-negative variant calls.
  + Detect by:
    - Examine read quality scores (Phred scores) using FastQC or fastp
    - Check per-base quality plots for drops at 3′ or 5′ ends
    - Assess proportion of bases above quality thresholds (e.g., Q20 or Q30)
* **Non-uniform/inconsistent coverage** – For various reasons, including that certain sequencing platforms or library preparation methods preferentially sequence regions with specific GC content. Regions with very high or very low GC may have low coverage. Overall, some genomic regions may have much higher or lower read depth than others, which can cause assembly gaps, missed variants, or misrepresentation of repetitive elements.
  + Detect by:
    - Coverage plots across the genome (samtools depth, bedtools)
    - Plotting coverage versus GC content across the genome (e.g., using Qualimap, samtools depth, or bedtools genomecov)
    - Regions with unusually high or low depth compared to the genome average
* **Mislabelling / Sample Mix-up** – The wrong sample identity can lead to completely incorrect conclusions, e.g., misassigned outbreak strains, incorrect AMR profiles, or wrong phylogenetic placement. Standard QC metrics (coverage, N50, read quality, contamination) **cannot reliably detect this**, because the data itself may be technically flawless.
  + Detect by:
    - Compare **genomic fingerprints** to previously sequenced isolates from the same source (e.g., MLST profiles, core-genome SNP distances).
    - Check for **unexpected relatedness** or incongruence with metadata (e.g., sample origin, collection date).
    - Include **negative and positive controls** in sequencing runs to help flag cross-run inconsistencies.

## Software to Perform Quality Control

<a name="tools-tbl"></a>**Table 2. Bioinformatics QC Tools, Metrics, and Purposes.**

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| **Tool Number** | **Purpose** | **Tool examples** | **Metrics to check** | **Purpose** |
| **1** | Raw read quality | FastQC, fastp | - Read quality scores (Phred scores)  - Per-base quality  - GC content  - (Adapter contamination, sequence duplication levels) | Identify and remove low-quality or contaminated reads before assembly |
| **2** | Trimming/filtering of reads | Trimmomatic, fastp |  | Trim adapters, low-quality bases, or reads below quality thresholds |
| **3** | Contamination check | Kraken, CheckM, QUAST (secondary), FastQC (adapter/foreign reads) | - Taxonomic composition of reads or contigs  - Unexpected GC content or contig anomalies  - Duplicate marker genes (CheckM) | Detect foreign DNA in raw reads of assemblies that could compromise analysis |
| **4** | Assembly evaluation | QUAST, CheckM | - Total assembly length  - Number of contigs  - N50  - GC content  - Total contig length / CDS length  - Contamination estimates | Assess completeness, contiguity, and accuracy of the assembly. |
| **5** | Read mapping statistics | BWA, Bowtie2 + samtools flagstat | - Coverage (breadth) and depth  - Number / proportion of reads mapping back | Verify that most reads are represented in the assembly and identify coverage gaps. |
| **6** | Variant-level quality (if performing variant analyses) | bcftools, GATK, FreeBayes | - Base quality at variant positions  - Strand bias  - Minimum read support per strand | Ensure variant calls are reliable and not artifacts. |

## Selecting the Best Assembly

For the same sample, when multiple sequencing runs or assemblies pass established QC metrics, the choice of which sequence to use should consider overall assembly quality, completeness, and coverage of regions of interest. Prefer assemblies with high chromosomal contiguity (e.g., higher N50 for the main chromosome) and a total genome length consistent with the expected reference. The presence of additional contigs corresponding to plasmids or other accessory elements should not be penalized, as these are biologically relevant. Sequences with more uniform coverage, higher read quality, and lower contamination are preferred. Additionally, prioritize assemblies that provide the most complete representation of genomic regions relevant to the application, such as MLST loci or antimicrobial resistance genes. When differences between assemblies are minor, select the sequence with stronger support for variant calling and other downstream analyses, and document the rationale for transparency and reproducibility. **Consultation with the laboratory team is recommended to interpret QC results, discuss potential sources of error, and consider iterative wet-lab or dry-lab improvements for future sequencing runs.**

![A diagram of a diagram

AI-generated content may be incorrect.](data:image/png;base64...)

## Organism-Specific Considerations

Bioinformatics analysis can be applied to a wide variety of microbes. While many QC and analysis methods are broadly applicable, some require organism-specific considerations.

The following are examples of well-studied organisms with genomic characteristics that may influence QC and analysis:

* ***Staphylococcus aureus*** many insertion sequences, plasmids, and phages leading to fragmented assemblies, misassemblies
* ***Mycobacterium tuberculosis*** (~65% GC) – Prone to GC bias leading to uneven coverage, missing regions with Illumina. Also, assemblies may miss or incorrectly place repetitive PE/PPE gene families.
* **Enterobacteriaceae (*E. coli*, Klebsiella, Salmonella)** – Large accessory genome and plasmids, leading to ambiguous AMR gene placement and generally more fragmented assemblies.
* ***Neisseria gonorrhoeae*, *Campylobacter jejuni*** – Naturally high recombination which complicates SNP-based phylogenies. May appear as “mixed” or contaminated if not interpreted carefully.
* ***Treponema pallidum***, ***Helicobacter pylori*** – Low biomass/difficult to extract DNA can lead to uneven coverage.
* ***Bordetella pertussis*** – Generally more fragmented assemblies due to the presence of very high copy number insertion sequence (IS) elements, when sequenced using short read technology.

<a name="example-sequences"></a> **Table 3. Example Quality Control Metrics Flagged or Failing for Bacterial Whole-Genome Assemblies.** This table provides illustrative examples of QC metrics that did not meet expected thresholds, highlighting potential issues in assembly contiguity, completeness, coverage, base composition, and contamination for bacterial isolates.

|  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **SRA ID** | **Organism** | **N50** | **Total assembly length** | **GC Content** | **Number of contigs** | **Depth** | **Coverage** | **Complete-ness** | **Contamina-tion** | **# coding sequences** |
| SRR1105999 | *Streptococcus pyogenes* | x |  |  | x |  |  |  |  |  |
| SRR30077323 | *Streptococcus pyogenes* |  |  |  |  |  |  | x |  |  |
| SRR2081599 | *Streptococcus pyogenes* |  |  |  |  |  |  |  | x |  |
| ERR13170183 | *Streptococcus pyogenes* | x | x |  |  |  |  | x |  | x |
| ERR2041950 | *Streptococcus pyogenes* | x | x | x | x |  |  |  | x | x |
| ASM4039145v1 | *Escherichia coli* |  | x | x | x | x | x |  | x |  |
| AUSP0129 | *Escherichia coli* |  | x | x |  |  |  |  |  |  |
| SAMN18191398 | *Escherichia coli* | x | x |  | x |  |  |  | x | x |
| SAMD00111247 | *Escherichia coli* | x | x | x |  |  |  | x |  | x |
| ASM77649v1 | *Escherichia coli* |  | x | x | x |  | x |  |  |  |
| SAL\_XA5511 | *Salmonella enterica* |  |  |  | x |  |  |  | x |  |
| SAMN12571641 | *Salmonella enterica* | x |  |  | x |  |  |  |  |  |
| SAMN24849049 | *Salmonella enterica* | x |  |  | x |  |  | x |  | x |
| SAMN33565140 | *Salmonella enterica* | x | x |  | x |  |  | x | x | x |
| SAMN11086324 | *Salmonella enterica* |  | x |  | x |  |  |  |  | x |
| GCA\_001058635.1\_ASM105863v1\_genomic | *Enterococcus faecium* | x | x | x | x |  |  |  | x |  |
| GCA\_964234405.1\_24-559 | *Enterococcus faecium* |  | x | x |  |  |  |  | x |  |