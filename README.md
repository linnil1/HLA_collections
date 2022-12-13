# HLA collections

Run HLA typing by different tools.

Tools includes

* arcasHLA=0.5.0 (RNA, https://github.com/RabadanLab/arcasHLA)
* bwakit=0.7.17 (https://github.com/lh3/bwa/tree/master/bwakit)
* graphtyper=2.7.5 (https://github.com/DecodeGenetics/graphtyper)
* OptiType=1.3.5 (RNA, https://github.com/FRED-2/OptiType)
* HISAT-genotype=1.3.3 (https://daehwankimlab.github.io/hisat-genotype/)
* HLA-HD=1.5.0 (manually download, https://www.genome.med.kyoto-u.ac.jp/HLA-HD)
* HLA-LA=1.0.3 (https://github.com/DiltheyLab/HLA-LA)
* HLAminer=1.4 (https://github.com/bcgsc/HLAminer)
* HLAscan=2.1.4 (https://github.com/SyntekabioTools/HLAscan)
* Kourami=0.9.6 (https://github.com/Kingsford-Group/kourami)
* POLYSOLVER=v4 (https://software.broadinstitute.org/cancer/cga/polysolver)
* seq2HLA=2.2--2 (RNA, https://github.com/TRON-Bioinformatics/seq2HLA)
* soap-HLA=1.0.0 (https://github.com/adefelicibus/soap-hla)
* HLA-VBSeq=11/22/2018 (https://nagasakilab.csml.org/hla/)
* xHLA=34221ea (https://github.com/humanlongevity/HLA)

Some tools are able to update the index to the latest IMGT-HLA version

* arcasHLA (Default 3.24.0 -> Current 3.49.0)
* bwakit (Cannot change)
* graphtyper (Provided 3.23.0, cannot change)
* OptiType (Provided 3.14.0, cannot change)
* HISAT-genotype (Provided 3.26.0 -> Max 3.43.0)
* HLA-HD (Default 3.32.0 -> Current 3.49.0)
* HLA-LA (Unknwon, cannot change)
* HLAminer (Provided 3.33.0 -> Current 3.49.0)
* HLAscan (Unknwon, cannot change, not open-source)
* Kourami (Provided 3.24.0 -> Current 3.49.0)
* POLYSOLVER (Provided 3.10.0, cannot change)
* seq2HLA (Unknwon, cannot change)
* soap-HLA (Provided 3.9.0, cannot change)
* HLA-VBSeq (Provided 3.31.0 -> Current 3.49.0)
* xHLA=34221ea (Unknwon -> Current 3.49.0)


## Requirements
* python>=3.10
* podman
* Python package (by pip install)
    * pandas


## Usage

Our pipeline started from fastq (pair-end reads).

``` bash
mkdir -p data
# copy your fastq to here
cp xxxx data/cohort_name.xxxxx.read.1.fq.gz
cp xxxx data/cohort_name.xxxxx.read.2.fq.gz
python pipeline.py data/cohort_name.xxxxx --version 3.49.0  # run one sample

# see more options via
python pipeline.py --help

# test all pipeline via
python pipeline.py example
```


## Example Output
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


## Notice
I not sure all the index building processes are correct.
Some of them are undocumented, I roughly review the code and find out the way to do it.

Also, the result reading part is written in the easily way.
In some cases, you should customize the code, add some features about
filtering, thresholding on tool provided parameters (e.g. abundance).


Some HLA tools are not considered as TODO:
* HLA-PRG (The precessor of HLA-LA)
* SNP2HLA (SNP array data as input)
* HLAreporter (The resources is GONE)


These tools require manully download requests:
* HLA-HD
* PHLAT
* PolyPheMe


TODO:
* HLAforest (https://github.com/FNaveed786/hlaforest)
* HLAssign (https://www.ikmb.uni-kiel.de/resources/download-tools/software/hlassign)
* HLAProfiler (https://github.com/ExpressionAnalysis/HLAProfiler)
* ATHLATES (https://github.com/cliu32/athlates)


## Reference
* Rose Orenbuch, Ioan Filip, Devon Comito, Jeffrey Shaman, Itsik Pe’er, Raul Rabadan, arcasHLA: high-resolution HLA typing from RNAseq, Bioinformatics, Volume 36, Issue 1, 1 January 2020, Pages 33–40, https://doi.org/10.1093/bioinformatics/btz474
* Li, Heng. "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM." arXiv preprint arXiv:1303.3997 (2013).
* Hannes P. Eggertsson, Hakon Jonsson, Snaedis Kristmundsdottir, Eirikur Hjartarson, Birte Kehr, Gisli Masson, Florian Zink, Kristjan E. Hjorleifsson, Aslaug Jonasdottir, Adalbjorg Jonasdottir, Ingileif Jonsdottir, Daniel F. Gudbjartsson, Pall Melsted, Kari Stefansson, Bjarni V. Halldorsson. Graphtyper enables population-scale genotyping using pangenome graphs. Nature Genetics 49, 1654–1660 (2017). doi:10.1038/ng.3964
* Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
* András Szolek, Benjamin Schubert, Christopher Mohr, Marc Sturm, Magdalena Feldhahn, Oliver Kohlbacher, OptiType: precision HLA typing from next-generation sequencing data, Bioinformatics, Volume 30, Issue 23, 1 December 2014, Pages 3310–3316, https://doi.org/10.1093/bioinformatics/btu548
* Kawaguchi, S., Higasa, K., Shimizu, M., Yamada, R., & Matsuda, F. (2017). HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. Human mutation, 38(7), 788–797. https://doi.org/10.1002/humu.23230
* Alexander T Dilthey, Alexander J Mentzer, Raphael Carapito, Clare Cutland, Nezih Cereb, Shabir A Madhi, Arang Rhie, Sergey Koren, Seiamak Bahram, Gil McVean, Adam M Phillippy, HLA\*LA—HLA typing from linearly projected graph alignments, Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4394–4396, https://doi.org/10.1093/bioinformatics/btz235
* Warren, R.L., Choe, G., Freeman, D.J. et al. Derivation of HLA types from shotgun sequence datasets. Genome Med 4, 95 (2012). https://doi.org/10.1186/gm396
* Ka, S., Lee, S., Hong, J. et al. HLAscan: genotyping of the HLA region using next-generation sequencing data. BMC Bioinformatics 18, 258 (2017). https://doi.org/10.1186/s12859-017-1671-3
* Lee, H., Kingsford, C. Kourami: graph-guided assembly for novel human leukocyte antigen allele discovery. Genome Biol 19, 16 (2018). https://doi.org/10.1186/s13059-018-1388-2
* Shukla, S., Rooney, M., Rajasagi, M. et al. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol 33, 1152–1158 (2015). https://doi.org/10.1038/nbt.3344
* Boegel, S., Löwer, M., Schäfer, M. et al. HLA typing from RNA-Seq sequence reads. Genome Med 4, 102 (2013). https://doi.org/10.1186/gm403
* Cao H, Wu J, Wang Y, Jiang H, Zhang T, et al. (2013) An Integrated Tool to Study MHC Region: Accurate SNV Detection and HLA Genes Typing in Human MHC Region Using Targeted High-Throughput Sequencing. PLOS ONE 8(7): e69388. https://doi.org/10.1371/journal.pone.0069388`
* Nariai, N., Kojima, K., Saito, S. et al. HLA-VBSeq: accurate HLA typing at full resolution from whole-genome sequencing data. BMC Genomics 16 (Suppl 2), S7 (2015). https://doi.org/10.1186/1471-2164-16-S2-S7
* Xie, Chao, et al. "Fast and accurate HLA typing from short-read next-generation sequence data with xHLA." Proceedings of the National Academy of Sciences 114.30 (2017): 8059-8064.


## LICENSE

MIT for the script in the repo.

But you have to check the licenses in each tools.
