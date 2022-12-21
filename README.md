# HLA collections

Run HLA typing by different tools.

Tools includes

* arcasHLA=0.5.0 (RNA, https://github.com/RabadanLab/arcasHLA)
* ATHLATES=2014_04_26 (RNA, https://www.broadinstitute.org/viral-genomics/athlates, https://github.com/cliu32/athlates)
* bwakit=0.7.17 (https://github.com/lh3/bwa/tree/master/bwakit)
* graphtyper=2.7.5 (https://github.com/DecodeGenetics/graphtyper)
* HISAT-genotype=1.3.3 (https://daehwankimlab.github.io/hisat-genotype/)
* HLAforest=75edb46 (RNA, https://github.com/FNaveed786/hlaforest)
* HLA-HD=1.5.0 (download request, https://www.genome.med.kyoto-u.ac.jp/HLA-HD)
* HLA-LA=1.0.3 (https://github.com/DiltheyLab/HLA-LA)
* HLAminer=1.4 (https://github.com/bcgsc/HLAminer)
* HLAProfiler=1.0.5 (RNA, https://github.com/ExpressionAnalysis/HLAProfiler)
* HLAscan=2.1.4 (https://github.com/SyntekabioTools/HLAscan)
* HLAssign=2015 (https://www.ikmb.uni-kiel.de/resources/download-tools/software/hlassign)
* Kourami=0.9.6 (https://github.com/Kingsford-Group/kourami)
* OptiType=1.3.5 (RNA, https://github.com/FRED-2/OptiType)
* PHLAT=1.1 (https://sites.google.com/site/phlatfortype)
* POLYSOLVER=v4 (https://software.broadinstitute.org/cancer/cga/polysolver)
* seq2HLA=2.2--2 (RNA, https://github.com/TRON-Bioinformatics/seq2HLA)
* SOAP-HLA=1.0.0 (https://github.com/adefelicibus/soap-hla)
* STC-Seq=1.0 (https://ngdc.cncb.ac.cn/biocode/tools/BT007068)
* HLA-VBSeq=11/22/2018 (https://nagasakilab.csml.org/hla/)
* xHLA=34221ea (https://github.com/humanlongevity/HLA)

Some tools are able to update the index to the latest IMGT-HLA version
(https://www.ebi.ac.uk/ipd/imgt/hla/)

* arcasHLA (Default 3.24.0 -> Current 3.49.0)
* ATHLATES (Provided 3.9.0 -> Max 3.31.0)
* bwakit (Cannot change)
* graphtyper (Provided 3.23.0, cannot change)
* HISAT-genotype (Provided 3.26.0 -> Max 3.43.0)
* HLAforest=75edb46 (Provided 3.10.0 -> Current 3.49.0)
* HLA-HD (Default 3.32.0 -> Current 3.49.0)
* HLA-LA (Unknwon, cannot change)
* HLAminer (Provided 3.33.0 -> Current 3.49.0)
* HLAProfiler (Provided 3.24.0, maybe can change but require more time to try)
* HLAscan (Unknwon, cannot change, not open-source)
* HLAssign (Provided 3.12.0 -> Current 3.49.0)
* Kourami (Provided 3.24.0 -> Current 3.49.0)
* OptiType (Provided 3.14.0, cannot change)
* PHLAT (Provided 3.9.0, cannot change)
* POLYSOLVER (Provided 3.10.0, cannot change)
* seq2HLA (Unknwon, cannot change)
* SOAP-HLA (Provided 3.9.0, cannot change)
* STC-Seq (Provided 3.26.0, cannot change)
* HLA-VBSeq (Provided 3.31.0 -> Current 3.49.0)
* xHLA=34221ea (Unknwon -> Current 3.49.0)


## Requirements
* python>=3.10
* podman (You can use docker if replace 'podman' in `pipeline.py` to 'docker')
* Python package (by `pip install`)
    * pandas


## Usage

Our pipeline starts from fastq (pair-end reads),

and run the HLA typing process.

The tool installation is included in dokcer image from Dockerfile or online registry.

The tool index will be downloaded and created by our pipeline when the index is not existed.

All the individual tool's result will be transformed and merged into `{sample_name}.hla_result_merge.tsv`.

``` bash
mkdir -p data
# copy your fastq to here
cp xxxx data/cohort_name.xxxxx.read.1.fq.gz
cp xxxx data/cohort_name.xxxxx.read.2.fq.gz
python pipeline.py data/cohort_name.xxxxx --tools hisat                   # select hisat-genotyping to typing your sample
python pipeline.py data/cohort_name.xxxxx --tools hisat --version 3.49.0  # select and using IMGT/HLA 3.49.0 index if available.

# see more options via
python pipeline.py --help

# test all pipeline via
python pipeline.py example
```


## Example Data folder
```
data/NA12878.bwa_bwakit_hs37_fa.bam
data/NA12878.bwa_bwakit_hs37_fa.bam.bai
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.bam
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.read.1.fq.gz
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.read.2.fq.gz
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.sortn.bam
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.est.call.hla_result.tsv
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.est.call.txt
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.est.txt
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.sam
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3490.est.call.hla_result.tsv
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3490.est.call.txt
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3490.est.txt
data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3490.sam
data/NA12878.bwa_bwakit_hs38_fa.bam
data/NA12878.bwa_bwakit_hs38_fa.bam.bai
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.fq.gz
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.hla
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.hla.full
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.hla_result.tsv
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.json
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.tsv
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_hla_3490.tsv.dna
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.fq.gz
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.hla
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.hla.full
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.hla_result.tsv
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.json
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.tsv
data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin.tsv.dna
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.addunmap.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.addunmap.bam.bai
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.addunmap.hlala_hlala_PRG_MHC_GRCh38_withIMGT
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.addunmap.hlala_hlala_PRG_MHC_GRCh38_withIMGT.hla_result.tsv
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.bam.bai
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.hlascan_hlascan.HLA-A.txt
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.hlascan_hlascan.HLA-B.txt
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.hlascan_hlascan.HLA-C.txt
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.hlascan_hlascan_merge.hla_result.tsv
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess_extract_1.fq.gz
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess_extract_2.fq.gz
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess_on_KouramiPanel.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.bam.bai
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_A.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_B.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_C.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_DQA1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_DQB1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call_DRB1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call.log
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call.result
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3240.call.result.hla_result.tsv
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.bam
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.bam.bai
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_A.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_B.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_C.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_DQA1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_DQB1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call_DRB1.typed.fa
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call.log
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call.result
data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.kourami_preprocess.panel_kourami_hla_3490.call.result.hla_result.tsv
data/NA12878.bwakit_bwakit_hs38DH_fa.hla.all
data/NA12878.bwakit_bwakit_hs38DH_fa.hla.top
data/NA12878.bwakit_bwakit_hs38DH_fa.hla.top.hla_result.tsv
data/NA12878.bwakit_bwakit_hs38DH_fa.log.bwamem
data/NA12878.bwakit_bwakit_hs38DH_fa.log.hla
data/NA12878.bwakit_bwakit_hs38DH_fa.sh
data/NA12878.download.cram
data/NA12878.hisat_hisat_hla_3260
data/NA12878.hisat_hisat_hla_3260.hla_result.tsv
data/NA12878.hisat_hisat_hla_3430
data/NA12878.hisat_hisat_hla_3430.hla_result.tsv
data/NA12878.hla_result.bwa_bwakit_hs37_fa_vbseq_extract_vbseq_vbseq_hla_3310_est_call.tsv
data/NA12878.hla_result.bwa_bwakit_hs37_fa_vbseq_extract_vbseq_vbseq_hla_3490_est_call.tsv
data/NA12878.hla_result.bwa_bwakit_hs38_fa_xhla_xhla_hla_3490.tsv
data/NA12878.hla_result.bwa_bwakit_hs38_fa_xhla_xhla_origin.tsv
data/NA12878.hla_result.bwakit_bwakit_hs38DH_fa_aln_sort_addunmap_hlala_hlala_PRG_MHC_GRCh38_withIMGT.tsv
data/NA12878.hla_result.bwakit_bwakit_hs38DH_fa_aln_sort_hlascan_hlascan_merge.tsv
data/NA12878.hla_result.bwakit_bwakit_hs38DH_fa_aln_sort_kourami_preprocess_panel_kourami_hla_3240_call_result.tsv
data/NA12878.hla_result.bwakit_bwakit_hs38DH_fa_aln_sort_kourami_preprocess_panel_kourami_hla_3490_call_result.tsv
data/NA12878.hla_result.bwakit_bwakit_hs38DH_fa_hla_top.tsv
data/NA12878.hla_result.hisat_hisat_hla_3260.tsv
data/NA12878.hla_result.hisat_hisat_hla_3430.tsv
data/NA12878.hla_result.hlascan_hlascan_merge.tsv
data/NA12878.hla_result_merge.csv
data/NA12878.hla_result.vbseq_vbseq_hla_3310_est_call.tsv
data/NA12878.hla_result.vbseq_vbseq_hla_3490_est_call.tsv
data/NA12878.hlascan_hlascan.HLA-A.txt
data/NA12878.hlascan_hlascan.HLA-B.txt
data/NA12878.hlascan_hlascan.HLA-C.txt
data/NA12878.hlascan_hlascan_merge.hla_result.tsv
data/NA12878.read.1.fq.gz
data/NA12878.read.2.fq.gz
data/NA12878.sortn.bam
data/NA12878.vbseq_vbseq_hla_3310.est.call.hla_result.tsv
data/NA12878.vbseq_vbseq_hla_3310.est.call.txt
data/NA12878.vbseq_vbseq_hla_3310.est.txt
data/NA12878.vbseq_vbseq_hla_3310.sam
data/NA12878.vbseq_vbseq_hla_3490.est.call.hla_result.tsv
data/NA12878.vbseq_vbseq_hla_3490.est.call.txt
data/NA12878.vbseq_vbseq_hla_3490.est.txt
data/NA12878.vbseq_vbseq_hla_3490.sam
```

## Example Result

e.g.  `data/NA12878.hla_result_merge.tsv`
```
    gene                 1                 2                                                                         name
0      A     A*03:01:01:01     A*03:01:01:01            data/NA12878.bwa_bwakit_hs38_fa.graphtyper_graphtyper_hla_3230.{}
0      A     A*01:01:01:01        A*11:01:01                                         data/NA12878.bwakit_bwakit_hs38DH_fa
0      A     A*01:01:01:01        A*11:01:01                          data/NA12878.bwa_bwakit_hs37_fa.polysolver_hla_3100
0      A  A*11:01:01:01:x1  A*01:01:01:01:x1  data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.est.call
0      A     A*01:01:01:01     A*11:01:01:01                             data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin
0      A        A*01:01:01        A*11:01:01           data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.phlat_phlat_hla_3080
0      A               NaN               NaN                    data/NA12878.athlates_athlates_hla_3430.novoalign.call.{}
1      B          B*08:19N        B*56:01:01                                         data/NA12878.bwakit_bwakit_hs38DH_fa
1      B        B*08:01:01        B*55:02:01                          data/NA12878.bwa_bwakit_hs37_fa.polysolver_hla_3100
1      B     B*08:01:01:01  B*56:01:01:04:x1  data/NA12878.bwa_bwakit_hs37_fa.vbseq_extract.vbseq_vbseq_hla_3310.est.call
1      B        B*08:01:01     B*56:01:01:01                             data/NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin
1      B        B*08:01:01        B*56:01:01           data/NA12878.bwakit_bwakit_hs38DH_fa.aln.sort.phlat_phlat_hla_3080
15     B        B*42:01:01     B*56:01:01:01            data/NA12878.bwa_bwakit_hs38_fa.graphtyper_graphtyper_hla_3230.{}
```


## Notice
I not sure all the index building processes are correct.
Some of them are undocumented, I roughly review the code and find out the way to do it.

Also, the result reading part is written in the easily way.
In some cases, you should customize the code, add some features about
filtering, thresholding on tool provided parameters (e.g. abundance, quality).
Maybe require fine-tunning the parameters in running code.

Not every tools has the ability to limit the threads and memory consumption,
you may test your data which has size as real data first.


Some HLA tools are not considered as TODO:
* HLA-PRG (The precessor of HLA-LA)
* SNP2HLA (SNP array data as input)
* HIBAG (SNP array)
* HLAreporter (Hard to run it, I don't have enough enthusiasm to fix them. https://github.com/jiaolongsun/hlareporter)
* HLAssign2 (Just visualization update)


These tools require manully fills the download requests:
* HLA-HD (After accepted, copy to `hlahd/hlahd.1.5.0.tar.gz`)
* PHLAT (I use docker, https://sites.google.com/site/projectphlat/Downloads)
* ATHLATES (After accepted, copy to `athlates/athlates.zip`, https://sites.google.com/site/projectphlat/Downloads)
* PolyPheMe


## Reference
* **arcasHLA** Rose Orenbuch, Ioan Filip, Devon Comito, Jeffrey Shaman, Itsik Pe’er, Raul Rabadan, arcasHLA: high-resolution HLA typing from RNAseq, Bioinformatics, Volume 36, Issue 1, 1 January 2020, Pages 33–40, https://doi.org/10.1093/bioinformatics/btz474
* **ATHLATES** Chang Liu, Xiao Yang, Brian Duffy, Thalachallour Mohanakumar, Robi D. Mitra, Michael C. Zody, John D. Pfeifer, ATHLATES: accurate typing of human leukocyte antigen through exome sequencing, Nucleic Acids Research, Volume 41, Issue 14, 1 August 2013, Page e142, https://doi.org/10.1093/nar/gkt481
* **bwa** **bwakit** Li, Heng. "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM." arXiv preprint arXiv:1303.3997 (2013).
* **Graphtyper** Hannes P. Eggertsson, Hakon Jonsson, Snaedis Kristmundsdottir, Eirikur Hjartarson, Birte Kehr, Gisli Masson, Florian Zink, Kristjan E. Hjorleifsson, Aslaug Jonasdottir, Adalbjorg Jonasdottir, Ingileif Jonsdottir, Daniel F. Gudbjartsson, Pall Melsted, Kari Stefansson, Bjarni V. Halldorsson. Graphtyper enables population-scale genotyping using pangenome graphs. Nature Genetics 49, 1654–1660 (2017). doi:10.1038/ng.3964
* **hisat** **hisat-genotype** Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
* **HLAforest** Kim, H. J., & Pourmand, N. (2013). HLA typing from RNA-seq data using hierarchical read weighting [corrected]. PloS one, 8(6), e67885. https://doi.org/10.1371/journal.pone.0067885
* **HLA-HD** Kawaguchi, S., Higasa, K., Shimizu, M., Yamada, R., & Matsuda, F. (2017). HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. Human mutation, 38(7), 788–797. https://doi.org/10.1002/humu.23230
* **HLA-LA** Alexander T Dilthey, Alexander J Mentzer, Raphael Carapito, Clare Cutland, Nezih Cereb, Shabir A Madhi, Arang Rhie, Sergey Koren, Seiamak Bahram, Gil McVean, Adam M Phillippy, HLA\*LA—HLA typing from linearly projected graph alignments, Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4394–4396, https://doi.org/10.1093/bioinformatics/btz235
* **HLAminer** Warren, R.L., Choe, G., Freeman, D.J. et al. Derivation of HLA types from shotgun sequence datasets. Genome Med 4, 95 (2012). https://doi.org/10.1186/gm396
* **HLAProfiler** Buchkovich, M.L., Brown, C.C., Robasky, K. et al. HLAProfiler utilizes k-mer profiles to improve HLA calling accuracy for rare and common alleles in RNA-seq data. Genome Med 9, 86 (2017). https://doi.org/10.1186/s13073-017-0473-6
* **HLAscan** Ka, S., Lee, S., Hong, J. et al. HLAscan: genotyping of the HLA region using next-generation sequencing data. BMC Bioinformatics 18, 258 (2017). https://doi.org/10.1186/s12859-017-1671-3
* **HLAssign** Michael Wittig, Jarl A. Anmarkrud, Jan C. Kässens, Simon Koch, Michael Forster, Eva Ellinghaus, Johannes R. Hov, Sascha Sauer, Manfred Schimmler, Malte Ziemann, Siegfried Görg, Frank Jacob, Tom H. Karlsen, Andre Franke, Development of a high-resolution NGS-based HLA-typing and analysis pipeline, Nucleic Acids Research, Volume 43, Issue 11, 23 June 2015, Page e70, https://doi.org/10.1093/nar/gkv184
* **HLA-VBSeq** Nariai, N., Kojima, K., Saito, S. et al. HLA-VBSeq: accurate HLA typing at full resolution from whole-genome sequencing data. BMC Genomics 16 (Suppl 2), S7 (2015). https://doi.org/10.1186/1471-2164-16-S2-S7
* **Kourami** Lee, H., Kingsford, C. Kourami: graph-guided assembly for novel human leukocyte antigen allele discovery. Genome Biol 19, 16 (2018). https://doi.org/10.1186/s13059-018-1388-2
* **OptiType** András Szolek, Benjamin Schubert, Christopher Mohr, Marc Sturm, Magdalena Feldhahn, Oliver Kohlbacher, OptiType: precision HLA typing from next-generation sequencing data, Bioinformatics, Volume 30, Issue 23, 1 December 2014, Pages 3310–3316, https://doi.org/10.1093/bioinformatics/btu548
* **PHLAT** Bai, Y., Wang, D., Fury, W. (2018). PHLAT: Inference of High-Resolution HLA Types from RNA and Whole Exome Sequencing. In: Boegel, S. (eds) HLA Typing. Methods in Molecular Biology, vol 1802. Humana Press, New York, NY. https://doi.org/10.1007/978-1-4939-8546-3_13
* **POLYSOLVER** Shukla, S., Rooney, M., Rajasagi, M. et al. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol 33, 1152–1158 (2015). https://doi.org/10.1038/nbt.3344
* **seq2HLA** Boegel, S., Löwer, M., Schäfer, M. et al. HLA typing from RNA-Seq sequence reads. Genome Med 4, 102 (2013). https://doi.org/10.1186/gm403
* **SOAP-HLA** Cao H, Wu J, Wang Y, Jiang H, Zhang T, et al. (2013) An Integrated Tool to Study MHC Region: Accurate SNV Detection and HLA Genes Typing in Human MHC Region Using Targeted High-Throughput Sequencing. PLOS ONE 8(7): e69388. https://doi.org/10.1371/journal.pone.0069388`
* **STC-Seq** Jiao, Y., Li, R., Wu, C. et al. High-sensitivity HLA typing by Saturated Tiling Capture Sequencing (STC-Seq). BMC Genomics 19, 50 (2018). https://doi.org/10.1186/s12864-018-4431-5
* **xHLA** Xie, Chao, et al. "Fast and accurate HLA typing from short-read next-generation sequence data with xHLA." Proceedings of the National Academy of Sciences 114.30 (2017): 8059-8064.

* **IMGT** James Robinson, Dominic J Barker, Xenia Georgiou, Michael A Cooper, Paul Flicek, Steven G E Marsh, IPD-IMGT/HLA Database, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D948–D955, https://doi.org/10.1093/nar/gkz950


## LICENSE

MIT for the script in the repo.

But you have to check the licenses in each tools.

And don't forget to cite this repo.
