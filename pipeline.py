"""
File: pipelilne.py of HLA_collections
Author: linnil1
Description: Write the build, running and read parsing sciprts
in the Python functions of each HLA tools.

Design Pattern:
* `folder  = createFolder("", "xxx")          `  Folder creation
* `folder  = xxxDownload(folder)              `  Download data that required to be used no matter the index is
* `folder  = checkAndBuildImage(folder)       `  Build docker image
* `index   = xxxBuild(folder, db_hla)         `  Build HLA index from db_hla (The path of IMGT'DB)
* `samples = xxxRun(samples_fq, index)        `  Run the tool by giving index
* `samples = xxxReadResult(samples)           `  Read the result from the tool's output and save as tsv
* `samples = renameResult(samples, samples_fq)`  Rename the hla result file for more organized

If you want to update the version of tool but not have differnet version tool at the same time,
just change the image name, folder name and xxx.dockerfile,
Mostly it will work.
"""

from glob import glob
from typing import Iterable, Callable
from pathlib import Path
from itertools import chain
from collections import defaultdict
import re
import gzip
import json
import uuid
import pandas as pd
import argparse
import subprocess


resources: dict[str, int] = {  # per sample
    "threads": 4,
    "memory": 7,  # unit: G
}

images = {
    # gerneral tools
    "bwa": "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
    "java": "docker.io/library/openjdk:11-jdk",
    "samtools": "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    # tools
    "arcas": "quay.io/biocontainers/arcas-hla:0.5.0--hdfd78af_1",
    "athlates": "localhost/linnil1/athlates:2014_04_26",
    "bwakit": "quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1",
    "graphtyper": "localhost/linnil1/graphtyper:2.7.5",
    "hisat": "localhost/linnil1/hisat2:1.3.3",
    "hisat2": "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "hlaforest": "localhost/linnil1/hlaforest",
    "hlahd": "localhost/linnil1/hlahd:1.5.0",
    "hlala": "quay.io/biocontainers/hla-la:1.0.3--hd03093a_0",
    "hlaminer": "localhost/linnil1/hlaminer:1.4",
    "hlaprofiler": "quay.io/biocontainers/hlaprofiler:1.0.5--hdfd78af_3",
    "hlascan": "localhost/linnil1/hlascan:2.1.4",
    "hlassign": "localhost/linnil1/hlassign",
    "kourami_preprocess": "localhost/linnil1/kourami_preprocess",
    "optitype": "quay.io/biocontainers/optitype:1.3.5--hdfd78af_2",
    "polysolver": "docker.io/sachet/polysolver:v4",
    "seq2hla": "localhost/linnil1/seq2hla:2.2--2",
    "soaphla": "localhost/linnil1/soaphla:1.0.0-pl526_3",
    "stcseq": "localhost/linnil1/stcseq:v1.0",
    "vbseq": "localhost/linnil1/vbseq:20181122",
    "xhla": "localhost/linnil1/xhla",
    "phlat_download": "docker.io/mgibio/phlat:1.1_withindex",
    "phlat": "docker.io/mgibio/phlat:1.1",
}

folders = {
    "arcas": "arcas",
    "athlates": "athlates",
    "bwakit": "bwakit",
    "graphtyper": "graphtyper",
    "hisat": "hisat",
    "hlaforest": "hlaforest",
    "hlahd": "hlahd",
    "hlala": "hlala",
    "hlaminer": "hlaminer",
    "hlaprofiler": "hlaprofiler",
    "hlascan": "hlascan",
    "hlassign": "hlassign",
    "kourami": "kourami",
    "optitype": "optitype",
    "phlat": "phlat",
    "seq2hla": "seq2hla",
    "soaphla": "soaphla",
    "stcseq": "stcseq",
    "vbseq": "vbseq",
    "xhla": "xhla",
}


def getThreads() -> int:
    """Get number of available threads"""
    return resources["threads"]


def setThreads(threads: int) -> None:
    """Set number of available threads"""
    global resources
    resources["threads"] = threads


def createFolder(base_folder: str = "", name: str = "") -> str:
    """Create folder at parent/base_folder"""
    output_name = base_folder + folders[name]
    Path(output_name).mkdir(exist_ok=True)
    return output_name


def name2Single(name: str) -> str:
    """Trun everything special character ('./') into '_'"""
    return name.replace("/", "_").replace(".", "_")


def runDocker(
    image: str, cmd: str, opts: str = "", chdir: str = ""
) -> subprocess.CompletedProcess[str]:
    """run docker container"""
    image = images.get(image, image)
    random_name = str(uuid.uuid4()).split("-", 1)[0]
    cmd_all = (
        f"podman run -it --rm -u root -w /app/{chdir} -v $PWD:/app "
        f"--name {random_name} "
        f"{opts} {image} {cmd}"
    )
    proc = runShell(cmd_all)
    return proc


def runShell(cmd: str) -> subprocess.CompletedProcess[str]:
    """Run command line"""
    print(cmd)
    proc = subprocess.run(
        cmd,
        shell=True,
        check=True,
        universal_newlines=True,
    )
    return proc


def buildImage(dockerfile: str, image: str, args: dict[str, str] = {}) -> None:
    """build docker image"""
    txt_arg = ""
    for key, value in args.items():
        txt_arg += f" --build-arg {key}={value} "
    runShell(f"podman build . -f {dockerfile} -t {image} {txt_arg}")


def checkImage(image: str) -> bool:
    """check image exists"""
    try:
        runShell(f"sh -c 'if [ ! $(podman image ls {image} -q) ]; then exit 1; fi'")
        return True
    except subprocess.CalledProcessError:
        return False


def checkAndBuildImage(folder: str, image: str = "") -> str:
    """Check image exists. If not existed, built it"""
    if not image:
        image = Path(folder).name
    if not checkImage(images[image]):
        buildImage(f"{image}.dockerfile", images[image], args={"folder": folder})
    return folder


def globAndRun(func: Callable[[str], str], input_name: str) -> str:
    """This is a temporary solution if not using my pipelilne"""
    for i in glob(input_name.replace(".{}", ".*.sh")):
        func(i[:-3])
    return input_name


def addSuffix(input_name: str, suffix: str) -> str:
    return input_name + suffix


def allelesToTable(
    alleles: Iterable[str], default_gene: list[str] = []
) -> pd.DataFrame:
    """
    Turn list of alleles into our hla_result format.

    Like so:
    ```
    gene	1	2	name
    A	A*11:126:01:01	A*01:01:01:01	data/NA12878.hisat_hisat_hla_3430
    B	B*08:01:01:45	B*56:01:01:03	data/NA12878.hisat_hisat_hla_3430
    C	C*01:02:01:13	C*01:02:01:06	data/NA12878.hisat_hisat_hla_3430
    ```
    """
    gene_allele: dict[str, list[str]] = defaultdict(list)
    for gene in default_gene:
        gene_allele[gene] = []

    for allele in alleles:
        gene = allele.split("*")[0]
        gene_allele[gene].append(allele)

    hla_list = []
    for gene, alleles in gene_allele.items():
        hla_list.append(
            {
                "gene": gene,
                "1": "",
                "2": "",
            }
        )
        if len(alleles) > 0:
            hla_list[-1]["1"] = alleles[0]
        if len(alleles) > 1:
            hla_list[-1]["2"] = alleles[1]
    df = pd.DataFrame(hla_list)
    print(df)
    return df


def bamSort(input_name: str, sam: bool = False) -> str:
    """sort the bamfile. If sam=True, remove samfile"""
    if sam:
        output_name = input_name
    else:
        output_name = input_name + ".sort"
    if Path(output_name + ".bam").exists():
        return output_name

    if sam:
        input_bam = input_name + ".sam"
    else:
        input_bam = input_name + ".bam"
    if not Path(input_bam).exists():
        raise ValueError(f"Not found {input_bam}")

    runDocker(
        "samtools", f"samtools sort -@{getThreads()} {input_bam} -o {output_name}.bam"
    )
    runDocker("samtools", f"samtools index {output_name}.bam")
    if sam:
        runShell(f"rm {input_bam}")
    return output_name


def bam2Fastq(input_name: str) -> str:
    """Extract paired fastq from bamfile"""
    output_name = input_name
    if Path(f"{output_name}.read.2.fq.gz").exists():
        return output_name
    runDocker(
        "samtools",
        f"samtools sort  -@{getThreads()} -n {input_name}.bam -o {output_name}.sortn.bam",
    )
    runDocker(
        "samtools",
        f"samtools fastq -@{getThreads()} -1 {output_name}.read.1.fq.gz -2 {output_name}.read.2.fq.gz "
        f"                                -0 /dev/null -s /dev/null -n {output_name}.sortn.bam",
    )
    return output_name


def downloadSample(folder: str = "data") -> str:
    """Download NA12878 example from HLA-LA"""
    name = folder + "/NA12878"
    if Path(f"{name}.read.2.fq.gz").exists():
        return name
    Path("data").mkdir(exist_ok=True)
    runShell(
        f"wget 'https://www.dropbox.com/s/xr99u3vqaimk4vo/NA12878.mini.cram?dl=0' -O {name}.download.cram"
    )
    runDocker(
        "samtools",
        f"samtools view -@{getThreads()} {name}.download.cram -o {name}.bam",
    )
    output_name = bam2Fastq(name)
    runShell(f"rm {name}.bam")
    return output_name


def bwaIndex(input_name: str) -> str:
    """bwa index"""
    if Path(f"{input_name}.bwt").exists():
        return input_name
    runDocker("bwa", f"bwa index {input_name}")
    return input_name


def bwaRun(input_name: str, index: str) -> str:
    """bwa mem"""
    output_name = input_name + ".bwa_" + name2Single(index)
    if Path(f"{output_name}.bam").exists():
        return output_name
    runDocker(
        "bwa",
        f"bwa mem -t {getThreads()} {index} "
        f"{input_name}.read.1.fq.gz {input_name}.read.2.fq.gz -o {output_name}.sam ",
    )
    output_name = bamSort(output_name, sam=True)
    return output_name


def downloadRef(folder: str = "bwakit", name: str = "hs38DH") -> str:
    """https://github.com/lh3/bwa/tree/master/bwakit"""
    if Path(f"{folder}/{name}.fa").exists():
        return f"{folder}/{name}.fa"

    runDocker("bwakit", f"run-gen-ref {name}")
    runShell(f"mv {name}.* {folder}")
    return f"{folder}/{name}.fa"


def bwakitRun(input_name: str, index: str) -> str:
    """https://github.com/lh3/bwa/tree/master/bwakit"""
    output_name = input_name + ".bwakit_" + name2Single(index)
    if Path(output_name + ".hla.top").exists():
        return output_name
    runDocker(
        "bwakit",
        f"run-bwamem -t {getThreads()} -H -o {output_name} -R '@RG\\tID:{input_name}\\tSM:{input_name}' "
        f"{index} {input_name}.read.1.fq.gz {input_name}.read.2.fq.gz > {output_name}.sh",
    )
    # BUG: dos format to linux format
    runShell(f"sed -i 's/\\r$//g' {output_name}.sh")
    runDocker("bwakit", f"bash {output_name}.sh")
    return output_name


def bwakitReadResult(input_name: str) -> str:
    """
    Turn bwakit HLA result format into our hla_result format

    bwakit HLA format:
    data/MMI001.bwakit_bwakit_hs38DH_fa.hla	HLA-A*02:07:01	HLA-A*11:01:01	0	0	6
    """
    output_name = input_name + ".hla.top.hla_result"
    if Path(output_name + ".tsv").exists():
        return output_name
    df = pd.read_csv(
        input_name + ".hla.top",
        sep="\t",
        names=["name", "1", "2", "unknown1", "unknown2", "unknown3"],
    )
    df = df[["name", "1", "2"]]
    df["1"] = df["1"].str.replace("HLA-", "")
    df["2"] = df["2"].str.replace("HLA-", "")

    df1 = allelesToTable(list(df["1"]) + list(df["2"]))
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def downloadHLA(folder: str = "", version: str = "3.49.0") -> str:
    """download database from https://github.com/ANHIG/IMGTHLA"""
    version_name = version.replace(".", "")
    if folder:
        folder = folder + "/hla_" + version_name
    else:
        folder = "hla_" + version_name
    if Path(folder).exists():
        return folder
    runShell(
        f"wget https://github.com/ANHIG/IMGTHLA/archive/refs/tags/v{version}-alpha.zip"
    )
    runShell(f"unzip v{version}-alpha.zip")
    runShell(f"rm v{version}-alpha.zip")
    runShell(f"mv IMGTHLA-{version}-alpha {folder}")
    return folder


def kouramiDownload(folder: str = "kourami") -> str:
    """https://github.com/Kingsford-Group/kourami"""
    if Path(folder + "/resources/hs38NoAltDH.fa").exists():
        return folder
    runShell(
        "wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9.6/kourami-0.9.6_bin.zip "
        f" -P {folder}"
    )
    runShell(f"unzip {folder}/kourami-0.9.6_bin.zip -d {folder}")
    runShell(f"mv {folder}/kourami-0.9.6/* {folder}")
    runDocker("bwa", f"bash {folder}/scripts/download_grch38.sh hs38NoAltDH")
    bwaIndex(folder + "/resources/hs38NoAltDH.fa")
    return folder


def kouramiBuild(folder: str, db_hla: str = "origin") -> str:
    """
    see https://github.com/Kingsford-Group/kourami/blob/master/preprocessing.md
    and formatIMGT.sh
    """
    # this line download v3.24 and save the db in {folder}/db
    if db_hla == "origin":
        output_name = f"{folder}/hla_3240"
    else:
        output_name = f"{folder}/{Path(db_hla).name}"

    if Path(f"{output_name}/All_FINAL_with_Decoy.fa.gz").exists():
        return output_name

    if db_hla == "origin":
        # scripts/download_panel.sh
        # runDocker("bwa", f"bash {folder}/scripts/download_panel.sh")
        runShell(
            f"wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9/kouramiDB_3.24.0.tar.gz -P {folder}"
        )
        runShell(f"tar -vxf {folder}/kouramiDB_3.24.0.tar.gz -C {folder}")
        runShell(f"mv {folder}/db {output_name}")
    else:
        # Kourami's bug?
        # runDocker("kourami", f"java -cp {folder}/build/Kourami.jar FormatIMGT {db_hla}/alignments/ . {output_name}")
        runShell(f"cp -r {db_hla}/alignments {folder}/tmp_alignments")
        runShell(
            f"wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/3240/alignments/Y_nuc.txt -O {folder}/tmp_alignments/Y_nuc.txt"
        )
        runShell(
            f"wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/3240/alignments/Y_gen.txt -O {folder}/tmp_alignments/Y_gen.txt"
        )
        runDocker(
            "java",
            f"java -cp {folder}/build/Kourami.jar FormatIMGT {folder}/tmp_alignments/ . {output_name}",
        )
        # runShell(f"rm -r {folder}/tmp_alignments")
        runShell(
            f"cat {output_name}/*.merged.fa {folder}/resources/HLA_decoys.fa | gzip > {output_name}/All_FINAL_with_Decoy.fa.gz"
        )
        runShell(f"cp {db_hla}/wmda/hla_nom_g.txt {output_name}")
    bwaIndex(output_name + "/All_FINAL_with_Decoy.fa.gz")
    return output_name


def kouramiPreprocess(input_name: str, index: str, kourami_folder: str = "") -> str:
    """https://github.com/Kingsford-Group/kourami"""
    output_name = input_name + ".kourami_preprocess"
    if Path(f"{output_name}._extract_2.fq.gz").exists():
        return output_name
    if not kourami_folder:
        kourami_folder = index + "/.."
    runDocker(
        "kourami_preprocess",
        f"bash {kourami_folder}/scripts/alignAndExtract_hs38DH.sh "
        f"-d {index} -r {kourami_folder}/resources/hs38NoAltDH.fa "
        f"{output_name}. "
        f"{input_name}.bam ",
    )
    return output_name


def kouramiRun(input_name: str, index: str, kourami_folder: str = "") -> str:
    """https://github.com/Kingsford-Group/kourami"""
    panel_name = input_name + ".panel_" + name2Single(index)
    output_name = panel_name + ".call"
    if Path(f"{output_name}.result").exists():
        return output_name
    # Because we use same preprocessed file
    # the last step of alignAndExtract_hs38DH generate extracted fastq
    if not kourami_folder:
        kourami_folder = index + "/.."
    runDocker(
        "kourami_preprocess",
        f"bwa mem -t {getThreads()} {index}/All_FINAL_with_Decoy.fa.gz "
        f"{input_name}._extract_1.fq.gz {input_name}._extract_2.fq.gz -o {panel_name}.sam ",
    )
    bamSort(panel_name, sam=True)
    runDocker(
        "java",
        f"java -jar {kourami_folder}/build/Kourami.jar -d {index} "
        f"{panel_name}.bam -o {output_name} ",
    )
    return output_name


def kouramiReadResult(input_name: str) -> str:
    """
    Turn Kourami HLA result into our hla_result.

    It's Format:
    A*11:01:01G	546	1.0	546	546	29.0	8.0	14.0
    """
    output_name = input_name + ".result.hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    df = pd.read_csv(
        input_name + ".result",
        sep="\t",
        names=[
            "allele",
            "matched_bases",
            "identity",
            "assembled_length",
            "matched_length",
            "bottleneck_sum",
            "bottleneck_1",
            "bottleneck_2",
        ],
    )
    df = allelesToTable(df["allele"], default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hisatDownload(folder: str = "hisat") -> str:
    """Download genome databse. see hisatgenotype_modules/hisatgenotype_typing_common.py"""
    if Path(f"{folder}/genome.fa.fai").exists():
        return folder
    # download_genotype_genome in hisatgenotype_modules/hisatgenotype_typing_common.py
    runShell(
        f"wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz -P {folder}"
    )
    runShell(f"tar -vxf {folder}/genotype_genome_20180128.tar.gz -C {folder}")
    # download_genome_and_index in hisatgenotype_modules/hisatgenotype_typing_common.py
    runShell(
        f"wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz -P {folder}"
    )
    runShell(f"tar -xvzf {folder}/grch38.tar.gz -C {folder}")
    runDocker(
        "hisat2", f"sh -c 'hisat2-inspect {folder}/grch38/genome > {folder}/genome.fa'"
    )
    runDocker("samtools", f"samtools faidx {folder}/genome.fa")
    return folder


def hisatBuild(folder: str, db_hla: str = "origin") -> str:
    """
    HISAT2 cannot read the IMGT-HLA higher than version 3.43.0.
    Which the A's exon length is not the same as A's nuc length
    (They have different number of gap insertions)
    """
    if db_hla == "origin":
        output_name = folder + "/hla_3260"
    else:
        output_name = folder + "/" + Path(db_hla).name
    Path(output_name).mkdir(exist_ok=True)

    if Path(f"{output_name}/hla.graph.8.ht2").exists():
        return output_name

    # link basic things
    files_genome = glob(folder + "/genotype_genome*")
    for file_genome in files_genome:
        runShell(f"ln -s ../{Path(file_genome).name} {output_name}/")
    runShell(f"    ln -s ../grch38                   {output_name}/")
    runShell(f"    ln -s ../genome.fa                {output_name}/")
    runShell(f"    ln -s ../genome.fa.fai            {output_name}/")

    # download HLA related things
    if db_hla == "origin":
        runShell(
            f"git clone https://github.com/DaehwanKimLab/hisatgenotype_db {output_name}/hisatgenotype_db"
        )
        runShell(f"cd {output_name}/hisatgenotype_db && git checkout d5d9b80")
    else:
        runShell(f"mkdir -p {output_name}/hisatgenotype_db/HLA")
        # Don't use all genes, risky to fail
        # runShell(f"cp -r {db_hla}/fasta {output_name}/hisatgenotype_db/HLA/")
        # runShell(f"cp -r {db_hla}/msf   {output_name}/hisatgenotype_db/HLA/")
        runShell(f"mkdir -p {output_name}/hisatgenotype_db/HLA/fasta")
        runShell(f"mkdir -p {output_name}/hisatgenotype_db/HLA/msf")
        genes = [
            "A",
            "B",
            "C",
            "DMA",
            "DMB",
            "DOA",
            "DOB",
            "DPA1",
            "DPB1",
            "DPB2",
            "DQA1",
            "DQB1",
            "DRA",
            "DRB",
            "E",
            "F",
            "G",
            "H",
            "HFE",
            "J",
            "K",
            "L",
            "MICA",
            "MICB",
            "P",
            "TAP1",
            "TAP2",
            "V",
        ]
        runShell(f"cp {db_hla}/hla.dat {output_name}/hisatgenotype_db/HLA")
        for gene in genes:
            runShell(
                f"cp -r {db_hla}/fasta/{gene}*  {output_name}/hisatgenotype_db/HLA/fasta/"
            )
            runShell(
                f"cp -r {db_hla}/msf/{gene}*  {output_name}/hisatgenotype_db/HLA/msf/"
            )

    # Build the index by running one example
    json.dump({"sanity_check": False}, open(f"{folder}/settings.json", "w"))
    with open(f"{folder}/read.fq", "w") as f:
        f.write(
            """
@B*82:02:01:02-746
ACCCACCCGGACTCAGAGTCTCCTCAGACGCCGAGATGCGGGTCACGGCACCCCGAACCCTCCTCCTGCTGCTCTGGGGGGCCCTGGCCCTGACCGAGACCTGGGCTGGTGAGTGCGGGGTCGGGAGGGAAATGGCCTCTGTGGGGAGGA
+
CCCGGGGGGGGGGJJJJJGCGJJ1G=GJJJJJGJCJGGJGGJ=JGJJJGJJGJJGJJGGGGJGGJGG=GCGGGGGC=G8CCGGGGCGGGGGGGGCGCGG1JCGGGGGGCGCCGG1GCGGCGGGGGGCGGCGCCGG=GGGGG1GCGGGCGC
@B*82:02:01:02-744
CCCTCACCCTGAGATGGGGTAAGGAGGGGGATGAGGGGTCATATCTCTTCTCAGGGAAAGCAGGAGCCCTTCTGGAGCCCTTCAGCAGGGTCAGGGCCCCTCATCTTCCCCTCCTTTCCCAGAGCCATCTTCCCAGTCCACCATCCCCAT
+
CC1CGGGGGGGGGGJGJJJJJJJGGJJGJGGJGGJGCJJJJGJJJJJJGGGGJ=GJJCGJGGGGJ=G=GGG==GGGCGGGC=GGCCGGGCGGG8GG=CGGJCGGCGGGGG1GGGGGGGGCGC1GCG=GGGCGGGCGGG8=C=GGGCGGGG
        """.strip()
        )
    runDocker(
        "hisat",
        f"hisatgenotype -z {output_name} --threads {getThreads()} "
        f"--base hla --out-dir /tmp -v -U {folder}/read.fq",
        opts=f"-v {Path(folder).absolute()}/settings.json:/opt/hisat-genotype/devel/settings.json",
    )
    return output_name


def hisatRun(input_name: str, index: str) -> str:
    """https://daehwankimlab.github.io/hisat-genotype/manual/"""
    output_name = input_name + ".hisat_" + name2Single(index)
    if len(glob(f"{output_name}/*.report")) > 0:
        return output_name
    runDocker(
        "hisat",
        f"hisatgenotype -z {index} --threads {getThreads()} "
        f"--base hla --out-dir {output_name} "
        f"--keep-alignment -v --keep-extract "
        f"-1 {input_name}.read.1.fq.gz "
        f"-2 {input_name}.read.2.fq.gz ",
    )

    # The folder has these three file
    # assembly_graph-hla.MMI001_read_1_fq_gz-hla-extracted-1_fq.report
    # MMI001.read.1.fq.gz-hla-extracted-1.fq.gz
    # MMI001.read.1.fq.gz-hla-extracted-2.fq.gz
    # and the current dir has these files
    # MMI001_read_1_fq_gz-hla-extracted-1_fq.bam*
    name = Path(f"{input_name}.read.1.fq.gz").name
    runShell(f"mv {name.replace('.', '_')}* {output_name}/")
    return output_name


def hisatReadResult(input_name: str) -> str:
    """
    Turn hisat HLA report into our hla_result format:

    It's format:
    ```
    1 ranked A*02:07:01 (abundance: 50.61%)
    2 ranked A*11:01:01:01 (abundance: 49.39%)
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    report_file = glob(f"{input_name}/*.report")[0]
    f = open(report_file)
    line_ranked = filter(lambda i: "ranked" in i, f)
    allele_abundance = map(
        lambda i: re.findall(r"ranked (.*) \(abundance: (.*)%\)", i)[0],  # type:ignore
        line_ranked,
    )
    df = pd.DataFrame(allele_abundance, columns=["allele", "abundance"])
    df["gene"] = df["allele"].str.split("*", expand=True)[0]
    df = (
        df.sort_values(["gene", "abundance"], ascending=[True, False])
        .groupby("gene")
        .head(2)
    )
    # print(df)
    df1 = allelesToTable(df["allele"], default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlascanBuild(folder: str = "hlascan", db_hla: str = "origin") -> str:
    """https://github.com/SyntekabioTools/HLAscan"""
    output_name = f"{folder}/origin"
    if Path(output_name + ".IMGT").exists():
        return output_name
    runShell(
        f"wget https://github.com/SyntekabioTools/HLAscan/releases/download/v2.0.0/dataset.zip -P {folder}"
    )
    runShell(f"unzip {folder}/dataset.zip -d {folder}")
    runShell(f"mv {folder}/db/HLA-ALL.IMGT {output_name}.IMGT")
    return output_name


def hlascanPreRun(input_name: str, index: str) -> str:
    """Split the gene before hlascanRun"""
    output_name = input_name + ".hlascan_" + name2Single(index)
    if Path(f"{output_name}.HLA-C").exists():
        return output_name + ".{}"
    genes = [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "HLA-F",
        "HLA-G",
        "MICA",
        "MICB",
        "HLA-DMA",
        "HLA-DMB",
        "HLA-DOA",
        "HLA-DOB",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DQA1",
        "HLA-DQB1",
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DRB5",
        "TAP1",
        "TAP2",
    ]
    for gene in genes:
        with open(f"{output_name}.{gene}.sh", "w") as f:
            if Path(input_name + ".bam").exists():
                if "hs37" in input_name:
                    version = "37"
                elif "hs38" in input_name:
                    version = "38"
                else:
                    raise ValueError("reference version cannot determine")
                f.write(
                    f"hlascan -g {gene} -t {getThreads()} -d {index}.IMGT "
                    f"-b {input_name}.bam -v {version} > {output_name}.{gene}.txt"
                )
            else:
                f.write(
                    f"hlascan -g {gene} -t {getThreads()} -d {index}.IMGT "
                    f"-l {input_name}.read.1.fq.gz -r {input_name}.read.2.fq.gz > {output_name}.{gene}.txt"
                )
    return output_name + ".{}"


def hlascanRun(input_name: str) -> str:
    """
    https://github.com/SyntekabioTools/HLAscan
    I'm not sure what input should given for hlascan:
    * selected Fastq
    * hs38 with alt (hs38DH)?
    * hs37 or hs37d5
    """
    output_name = input_name
    if Path(f"{output_name}.txt").exists():
        return output_name
    try:
        runDocker("hlascan", f"bash {input_name}.sh")
    except subprocess.CalledProcessError as e:
        pass
    return output_name


def hlascanReadResult(input_name: str) -> str:
    """
    Turn hlascan's HLA txt into our hla_result format.

    Note that the txt files is splitted by hla gene

    It's format:
    ```
    HLA gene : HLA-DPA1
    # of considered types : 40

    ----------- HLA-Types ----------
    [Type 1]	02:02:06	EX2_2.85366_53.3333	EX3_0_0	EX4_0_0
    [Type 2]	02:02:06	EX2_2.85366_53.3333	EX3_0_0	EX4_0_0
    ```
    """
    output_name = input_name.replace(".{}", "_merge") + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    alleles = []
    for gene in ["HLA-A", "HLA-B", "HLA-C"]:
        gene_name = gene.split("-")[1]
        for i in open(input_name.replace("{}", f"{gene}.txt")):
            if "[Type" in i:
                alleles.append(f"{gene_name}*{i.split()[2]}")

    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def vbseqBuild(folder: str, db_hla: str = "origin") -> str:
    """The document said it's able to use newest IMGT-HLA databse"""
    if db_hla == "origin":
        output_name = folder + "/hla_3310"
    else:
        output_name = folder + "/" + Path(db_hla).name
    if Path(f"{output_name}.all.fa").exists():
        return output_name
    # default version
    if db_hla == "origin":
        runShell(
            f"wget https://nagasakilab.csml.org/hla/Allelelist_v2.txt -O {output_name}.allelelist.txt"
        )
        runShell(
            f"wget https://nagasakilab.csml.org/hla/hla_all_v2.fasta -O {output_name}.all.fa"
        )
    else:
        # csv -> special txt
        df = pd.read_csv(open(f"{db_hla}/Allelelist.txt"), comment="#")
        df.to_csv(f"{output_name}.allelelist.txt", index=False, header=False, sep=" ")
        runShell(f"cp {db_hla}/hla_gen.fasta {output_name}.all.fa")
    bwaIndex(output_name + ".all.fa")
    return output_name


def vbseqPreprocess(input_name: str) -> str:
    output_name = input_name + ".vbseq_extract"
    if Path(f"{output_name}.bam").exists():
        return output_name
    runDocker(
        "samtools",
        f"""\
        samtools view {input_name}.bam \
            6:29907037-29915661 6:31319649-31326989 6:31234526-31241863 \
            6:32914391-32922899 6:32900406-32910847 6:32969960-32979389 6:32778540-32786825 \
            6:33030346-33050555 6:33041703-33059473 6:32603183-32613429 6:32707163-32716664 \
            6:32625241-32636466 6:32721875-32733330 6:32405619-32414826 6:32544547-32559613 \
            6:32518778-32554154 6:32483154-32559613 6:30455183-30463982 6:29689117-29699106 \
            6:29792756-29800899 6:29793613-29978954 6:29855105-29979733 6:29892236-29899009 \
            6:30225339-30236728 6:31369356-31385092 6:31460658-31480901 6:29766192-29772202 \
            6:32810986-32823755 6:32779544-32808599 6:29756731-29767588 \
            -o {output_name}.bam
        """,
    )
    return bam2Fastq(output_name)


def vbseqRun(input_name: str, index: str) -> str:
    """https://nagasakilab.csml.org/hla/"""
    output_name = input_name + ".vbseq_" + name2Single(index)
    output_name1 = output_name + ".est"
    output_name2 = output_name1 + ".call"
    if Path(f"{output_name2}.txt").exists():
        return output_name2
    runDocker(
        "bwa",
        f"bwa mem -t {getThreads()} -P -L 10000 -a {index}.all.fa "
        f"{input_name}.read.1.fq.gz {input_name}.read.2.fq.gz -o {output_name}.sam ",
    )
    runDocker(
        "vbseq",
        f"java -jar /opt/HLAVBSeq.jar {index}.all.fa "
        f"{output_name}.sam {output_name1}.txt --alpha_zero 0.01 --is_paired ",
    )
    runDocker(
        "vbseq",
        "perl /opt/parse_result.pl "
        f"{index}.allelelist.txt {output_name1}.txt > {output_name2}.txt",
    )
    return output_name2


def vbseqReadResult(input_name: str) -> str:
    """
    Read vbseq hla and coverage result into our hla_result format.

    It's format:
    ```
    A*02:03:01	0
    A*02:03:03	0
    A*02:03:04	0.996701649175412
    A*02:04	0
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    df = pd.read_csv(input_name + ".txt", sep="\t", names=["allele", "coverage"])

    if not len(df):
        df1 = allelesToTable([], default_gene=["A", "B", "C"])
        df1["name"] = input_name
        df1.to_csv(output_name + ".tsv", index=False, sep="\t")
        return output_name

    df["gene"] = df["allele"].str.split("*", expand=True)[0]
    df = df[df["coverage"] > 5]
    df = (
        df.sort_values(["gene", "coverage"], ascending=[True, False])
        .groupby("gene")
        .head(2)
    )

    df1 = allelesToTable(df["allele"], default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlalaBuild(folder: str = "hlala", db_hla: str = "origin") -> str:
    """https://github.com/DiltheyLab/HLA-LA"""
    output_name = f"{folder}/origin"
    if Path(f"{output_name}/serializedGRAPH").exists():
        return output_name
    runShell(
        f"wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz -P {folder}"
    )
    runShell(f"tar -vxf {folder}/PRG_MHC_GRCh38_withIMGT.tar.gz -C {folder}")
    runShell(f"mv {folder}/PRG_MHC_GRCh38_withIMGT {output_name}")
    runDocker(
        "hlala",
        "sh -c 'export PATH=/usr/local/opt/hla-la/bin:$PATH && "
        f"HLA-LA --workingDir ./ --action prepareGraph --PRG_graph_dir {output_name}'",
    )
    return output_name


def hlalaRun(input_name: str, index: str) -> str:
    """https://github.com/DiltheyLab/HLA-LA"""
    output_name = input_name + ".hlala_" + name2Single(index)
    if Path(f"{output_name}/data/hla/R1_bestguess_G.txt").exists():
        return output_name
    Path(output_name).mkdir(exist_ok=True)
    runDocker(
        "hlala",
        f"HLA-LA.pl --workingDir {output_name} --graph . --maxThreads {getThreads()} "
        f"--BAM {input_name}.bam --sampleID data ",
        opts=f" -v {Path(index).absolute()}:/usr/local/opt/hla-la/graphs/",
    )
    return output_name


def hlalaReadResult(input_name: str) -> str:
    """
    Read HLA-LA guess into our hla_result format

    It's format:
    ```
    Locus	Chromosome	Allele	Q1	...
    A	1	A*11:01:01G	1	-70	...
    A	2	A*01:01:01G	1	-70	...
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    df = pd.read_csv(f"{input_name}/data/hla/R1_bestguess_G.txt", sep="\t")
    df = df[df["AverageCoverage"] > 5]
    df1 = allelesToTable(df["Allele"], default_gene=["A", "B", "C"])
    # print(df1)
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def addUnmap(input_name: str) -> str:
    """
    HLA-LA has the step for extracting unmapped read
    and it require at least one unmapped reads.
    So I add this step to ensure HLA-LA doesn't fail.
    """
    output_name = input_name + ".addunmap"
    if Path(f"{output_name}.bam").exists():
        return output_name
    view = runDocker(
        "samtools", f"samtools view -h {input_name}.bam -o {output_name}.sam"
    )
    with open(output_name + ".sam", "a") as f:
        fields = ["ggggg", 77, "*", 0, 0, "*", "*", 0, 0, "A", "A", "AS:i:0 XS:i:0"]
        f.write("\t".join(map(str, fields)) + "\n")
        fields[1] = 141
        f.write("\t".join(map(str, fields)) + "\n")
    bamSort(output_name, sam=True)
    return output_name


def xhlaBuild(folder: str, db_hla: str = "origin") -> str:
    """
    The HLA data is saved in image.
    If you want to update the HLA index,
    you should build it by script proived in github,
    and wait for 7 hours to built.
    """
    if db_hla == "origin":
        output_name = folder + "/origin"
    else:
        output_name = f"{folder}/{Path(db_hla).name}"
    if Path(output_name).exists():
        return output_name

    if db_hla == "origin":
        runDocker("xhla", f"cp -r /opt/data {output_name}")
    else:
        output_name1 = output_name + "_tmp"
        runShell(f"git clone https://github.com/humanlongevity/HLA {output_name1}")
        runShell(f"cd {output_name1} && git checkout 34221ea")
        runShell(f"mkdir -p {output_name1}/raw")
        runShell(f"cp -r {db_hla}/alignments {output_name1}/raw")
        for i in ["nuc", "exon", "dna", "temp", "align"]:
            runShell(f"mkdir -p {output_name1}/data/{i}")
        runDocker("xhla", f"bash script/batch.sh", chdir=f"{output_name1}/data")
        runDocker(
            "xhla",
            "diamond makedb --in hla.faa --db hla.dmnd",
            chdir=f"{output_name1}/data/",
        )
        runShell(f"mv {output_name1}/data {output_name}")
        # runShell(f"rm -rf {output_name1}")
        return output_name
    return output_name


def xhlaRun(input_name: str, index: str) -> str:
    """
    https://github.com/humanlongevity/HLA
    Don't input hg19 or hs38DH bam file
    """
    output_name = input_name + ".xhla_" + name2Single(index)
    if Path(f"{output_name}.json").exists():
        return output_name
    # Full typing result is not in output_path. It's in hla-{id}.
    # see bin/run.py and bin/typer.sh
    id = Path(output_name).name
    runDocker(
        "xhla",
        f"typer.sh {input_name}.bam {id} full",
        opts=f"-v {Path(index).absolute()}:/opt/data ",
    )
    runShell(f"mv hla-{id}/{id}* {Path(output_name).parent}")
    runShell(f"rm hla-{id} -r")
    return output_name


def xhlaReadResult(input_name: str) -> str:
    """
    Read xhla json into our hla_result format

    Its format:
    ```
    {
     "sample_id": "NA12878.bwa_bwakit_hs38_fa.xhla_xhla_origin",
     "hla": {
      "alleles": [
       "A*01:01",
       "A*11:01",
       "B*08:01",
      ]
     }
    }
    or full
        type	full	total	good	rank
        A*01:01	A*01:01:01:01	3	3	1
        A*11:01	A*11:01:01:01	3	3	1
        B*08:01	B*08:01:01	2	2	1
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    # 2 fields
    # result = json.load(open(f"{input_name}.json"))
    # df = allelesToTable(result["hla"]["alleles"], default_gene=["A", "B", "C"])
    # 4 fields
    df_full = pd.read_csv(f"{input_name}.hla.full", sep="\t")
    # print(df_full)
    df = allelesToTable(df_full["full"], default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def graphtyperBuild(
    folder: str, db_hla: str = "origin", index_hs38: str = "bwakit/hs38.fa"
) -> str:
    """
    https://github.com/DecodeGenetics/graphtyper

    IMGT version 3.23.0 written in their paper.
    I'm not sure if database is updated.

    Ignore db_hla, It doesn't provided scripts for updating index
    """
    output_name = folder + "/hla_3230"
    if Path(f"{output_name}/hs38.fa.fai").exists():
        return output_name
    runShell(f"mkdir -p {output_name}")
    runShell(
        f"wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.5/HLA_data.tar.gz -P {output_name}"
    )
    runShell(f"tar -vxf {output_name}/HLA_data.tar.gz -C {output_name}")
    runDocker("samtools", f"samtools faidx {index_hs38}")
    runShell(f"ln -s ../../{index_hs38} {output_name}/hs38.fa")
    runShell(f"ln -s ../../{index_hs38}.fai {output_name}/hs38.fa.fai")
    return output_name


def graphtyperPreRun(input_name: str, index: str) -> str:
    """Split genes before graphtyperRun"""
    output_name = input_name + ".graphtyper_" + name2Single(index)
    if Path(output_name + ".HLA_A.sh").exists():
        return output_name + ".{}"
    id = Path(output_name).name
    df = pd.read_csv(
        f"{index}/regions.tsv", sep="\t", names=["chrom", "start", "end", "gene"]
    )
    for i in df.itertuples():
        with open(f"{output_name}.{i.gene}.sh", "w") as f:
            f.write(
                f"graphtyper genotype_hla {index}/hs38.fa --verbose --threads={getThreads()} "
                f"{index}/{i.gene}.vcf.gz --region={i.chrom}:{i.start}-{i.end} "
                f"--sam={input_name}.bam --output={output_name}"
            )
            f.write("\n")
            f.write(
                f"ln -s ../{output_name}/{i.chrom}/{i.start:09d}-{i.end:09d}.vcf.gz {output_name}.{i.gene}.vcf.gz"
            )
    return output_name + ".{}"


def graphtyperRun(input_name: str) -> str:
    """https://github.com/DecodeGenetics/graphtyper/wiki/HLA-genotyping"""
    output_name = input_name
    if Path(output_name + ".vcf.gz").exists():
        return output_name
    runDocker("graphtyper", f"sh {input_name}.sh")
    return output_name


def graphtyperReadResult(input_name: str) -> str:
    """
    Read graphtyper result in each vcfs into our hla_result format

    Its format:
    ```
    ...
    ##FILTER=<ID=LowPratio,Description="Ratio of PASSed calls was too low.">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878
    chr6    29724785        chr6:29724785:H.all     <HLA-F*01:03:01:01>     <HLA-F*01:01:01:01>,<HLA-F*01:01:02:01> 255     .       AC=1,1;AF=0.5,0.5;AN=2;MQ=0;NHet=0;NHomAlt=1;NHomRef=0;PASS_AC=1,1;PASS_AN=2;PASS_ratio=1;RefLen=19;VarType=H  GT:GQ:PL        1/2:99:255,255,200,150,0,200
    ```
    """
    output_name = input_name.replace(".{}", "_merge") + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name

    alleles: list[str] = []
    for vcf_file in glob(input_name.replace(".{}", ".*.vcf.gz")):
        for vcf_line in gzip.open(vcf_file, "rt"):
            if vcf_line.startswith("#"):
                continue
            # "H.2digit"   # 1 fields
            # "H.4digit"   # 2 fields
            if "H.all" not in vcf_line:  # full resolution
                continue
            vcf_field = vcf_line.split("\t")
            gts = [vcf_field[3], *vcf_field[4].split(",")]
            gts = [i.replace("<HLA-", "").replace(">", "") for i in gts]
            sample_gt_txt = vcf_field[9].split(":")[0].split("/")
            sample_gt = [int(i) for i in sample_gt_txt if i != "."]
            alleles.extend(gts[gt] for gt in sample_gt)
    df = allelesToTable(alleles)
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def optitypeBuild(folder: str = "optitype", db_hla: str = "origin") -> str:
    """Get index from dockerfile"""
    output_name = folder + "/hla_3140"
    if Path(output_name).exists():
        return output_name
    runDocker("optitype", f"cp -r /usr/local/bin/data/ {output_name}")
    with open(f"{output_name}/config.ini", "w") as f:
        f.write(
            f"""
[mapping]
razers3=/usr/local/bin/razers3
threads={getThreads()}

[ilp]
solver=glpk
threads={getThreads()}

[behavior]
deletebam=false
unpaired_weight=0
use_discordant=false
    """
        )
    return output_name


def optitypeRun(input_name: str, index: str) -> str:
    """https://github.com/FRED-2/OptiType"""
    output_name = input_name + ".optitype_" + name2Single(index)
    if Path(f"{output_name}._result.tsv").exists():
        return output_name
    parent = Path(output_name).parent
    name = Path(output_name).name
    runDocker(
        "optitype",
        f"OptiTypePipeline.py  --dna --verbose"
        f" --input {input_name}.read.1.fq.gz {input_name}.read.2.fq.gz"
        f" --outdir {parent} --prefix {name}.",
        opts=(
            f" -v {Path(index).absolute()}:/usr/local/bin/data "
            f" -v {Path(index).absolute()}/config.ini:/usr/local/bin/config.ini "
        ),
    )
    return output_name


def optitypeReadResult(input_name: str) -> str:
    """
    Read optitype json into our hla_result format

    Its format:
    ```
        A1      A2      B1      B2      C1      C2      Reads   Objective
    0       A*01:01 A*11:01 B*08:01 B*56:01 C*01:02 C*07:01 1529.0  1446.4
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    df = pd.read_csv(f"{input_name}._result.tsv", sep="\t")
    alleles = df[["A1", "A2", "B1", "B2", "C1", "C2"]].dropna(axis=1).values.tolist()[0]
    df1 = allelesToTable(filter(None, alleles), default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def arcasBuild(folder: str, db_hla: str) -> str:
    """ArcasHLA doesn't provide downloadable HLA"""
    output_name = f"{folder}/{Path(db_hla).name}"
    if Path(f"{output_name}/hla_partial.idx").exists():
        return output_name
    Path(output_name).mkdir(exist_ok=True)
    runShell(f"cp {db_hla}/hla.dat {output_name}/")
    runDocker(
        "arcas",
        "arcasHLA reference --rebuild",
        opts=(
            f" -v {Path(db_hla).absolute()}:/usr/local/share/arcas-hla-0.5.0-1/dat/IMGTHLA/:ro "
            f" -v {Path(output_name).absolute()}:/usr/local/share/arcas-hla-0.5.0-1/dat/ref/ "
        ),
    )
    return output_name


def arcasPreprocess(input_name: str) -> str:
    """https://github.com/RabadanLab/arcasHLA"""
    output_name = input_name + ".arcas_extract"
    # output_name = input_name + ".arcas_" + name2Single(index)
    if Path(f"{output_name}.read.2.fq.gz").exists():
        return output_name
    runDocker(
        "arcas",
        f"arcasHLA extract -t {getThreads()} -v " f"{input_name}.bam -o {output_name}",
    )
    name = Path(input_name).name
    runShell(
        f"ln -s ../{output_name}/{name}.extracted.1.fq.gz {output_name}.read.1.fq.gz"
    )
    runShell(
        f"ln -s ../{output_name}/{name}.extracted.2.fq.gz {output_name}.read.2.fq.gz"
    )
    return output_name


def arcasRun(input_name: str, index: str) -> str:
    """https://github.com/RabadanLab/arcasHLA"""
    output_name = input_name + ".arcas_" + name2Single(index)
    if len(glob(f"{output_name}/*.json")) >= 1:
        return output_name
    runDocker(
        "arcas",
        f"arcasHLA genotype -o {output_name} -t {getThreads()} -v "
        f"{input_name}.read.1.fq.gz {input_name}.read.2.fq.gz",
        opts=(
            f" -v {Path(index).absolute()}:/usr/local/share/arcas-hla-0.5.0-1/dat/ref:ro "
            f" -v {Path(index).absolute()}/hla.dat:/usr/local/share/arcas-hla-0.5.0-1/dat/IMGTHLA/hla.dat:ro "
        ),
    )
    return output_name


def arcasReadResult(input_name: str) -> str:
    """
    Read arcas json into our hla_result format

    Its format:
    ```
    {"A": ["A*01:01:71", "A*03:01:01"], "B": ["B*07
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    f = glob(f"{input_name}/*.genotype.json")[0]
    data = json.load(open(f))
    df = allelesToTable(
        chain.from_iterable(data.values()), default_gene=["A", "B", "C"]
    )
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def polysolverRun(input_name: str, index: str) -> str:
    """https://software.broadinstitute.org/cancer/cga/polysolver_run"""
    output_name = input_name + ".polysolver_hla_3100"
    if Path(f"{output_name}/winners.hla.nofreq.txt").exists():
        return output_name
    runDocker(
        "polysolver",
        f"bash /home/polysolver/scripts/shell_call_hla_type"
        f" {input_name}.bam Unknown 0 hg19 STDFQ 0 {output_name}",
    )
    return output_name


def polysolverReadResult(input_name: str) -> str:
    """
    Read polysolver txt into our hla_result format

    Its format:
    ```
    HLA-A   hla_a_01_01_01_01       hla_a_01_01_01_01
    HLA-B   hla_b_07_02_01  hla_b_82_02
    HLA-C   hla_c_01_02_01  hla_c_01_02_01
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    txt = pd.read_csv(f"{input_name}/winners.hla.nofreq.txt", sep="\t", header=None)
    alleles = txt[[1, 2]].values.flatten().tolist()
    alleles = [i.replace("hla_", "").upper().replace("_", ":") for i in alleles]
    alleles = [i.replace(":", "*", 1) for i in alleles]
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def soaphlaBuild(folder: str, db_hla: str = "origin") -> str:
    """
    https://github.com/adefelicibus/soap-hla
    Get index from git
    """
    output_name = f"{folder}/hla_3090"
    if Path(output_name).exists():
        return output_name
    runShell(f"git clone https://github.com/adefelicibus/soap-hla.git {folder}/soaphla")
    runShell(f"cd {folder}/soaphla && git checkout 407812d")
    runShell(f"cp -r {folder}/soaphla/data {output_name}")
    return output_name


def soaphlaRun(input_name: str, index: str) -> str:
    """https://github.com/adefelicibus/soap-hla"""
    output_name = input_name + ".soaphla_" + name2Single(index)
    name = Path(input_name).name.split(".")[0]
    if Path(f"{output_name}/{name}/{name}.type").exists():
        return output_name
    runDocker(
        "soaphla",
        f"MHC_autopipeline -i {input_name}.bam -od {output_name} -v hg19",
        opts=f" -v {Path(index).absolute()}:/opt/conda/share/database:ro ",
    )
    return output_name


def soaphlaReadResult(input_name: str) -> str:
    """
    Read soaphla typing result into our hla_result format

    Its format:
    ```
    ---     ---     50.00
    B*82:02 ---     99.03
    B*07:02:01      B*07:04 68.58
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    name = Path(input_name).name.split(".")[0]
    txt = pd.read_csv(
        f"{input_name}/{name}/{name}.type",
        sep="\t",
        names=["final_type", "second_type", "total_type_score", "tmp"],
    )
    alleles = [i for i in txt["final_type"] if i != "---"]
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlaminerBuild(folder: str, db_hla: str = "origin") -> str:
    """https://github.com/bcgsc/HLAminer"""
    if db_hla == "origin":
        output_name = f"{folder}/hla_3330"
        if Path(output_name + "/HLA-I_II_GEN.fasta").exists():
            return output_name
        runShell(
            f"wget https://github.com/bcgsc/HLAminer/releases/download/v1.4/HLAminer_1-4.tar.gz -P {folder}"
        )
        runShell(f"tar -vxf {folder}/HLAminer_1-4.tar.gz -C {folder}")
        runShell(f"cp -r {folder}/HLAminer-1.4/HLAminer_v1.4/database {output_name}")
    else:
        output_name = f"{folder}/{Path(db_hla).name}"
        if Path(output_name + "/HLA-I_II_GEN.fasta").exists():
            return output_name
        runShell(f"mkdir -p {output_name}")
        # see updateHLA-I_II_genomic.sh
        genes = ["A", "B", "C", "F", "G", "H", "DP*", "DQ*", "DR*"]
        with open(f"{output_name}/updateI_II_gen.sh", "w") as f:
            f.write("cat ")
            for gene in genes:
                f.write(f"{db_hla}/fasta/{gene}_gen.fasta ")
            f.write(
                ' | perl -ne \'chomp;if(/\>\S+\s+(\S+)/){print ">$1\\n";}else{print "$_\\n";}\' '
            )
            f.write(f"> {output_name}/HLA-I_II_GEN.fasta")
        runDocker("hlaminer", f"sh {output_name}/updateI_II_gen.sh")
        runShell(f"cp {db_hla}/wmda/hla_nom_p.txt {output_name}")
        runDocker("hlaminer", f"formatdb -p F -i {output_name}/HLA-I_II_GEN.fasta")

    # Index it (The index from downloaded files is very old bwa format)
    runDocker("bwa", f"bwa index {output_name}/HLA-I_II_GEN.fasta")
    # file path hacking
    runShell(f"mkdir -p HLAminer_HPRA_")
    return output_name


def hlaminerRun(input_name: str, index: str) -> str:
    """https://github.com/bcgsc/HLAminer"""
    output_name = input_name + ".hlaminer_" + name2Single(index)
    if Path(output_name + ".csv").exists():
        return output_name
    index_ref = f"{index}/HLA-I_II_GEN.fasta"
    runDocker(
        "bwa",
        f"bwa mem -t {getThreads()} {index_ref} "
        f"{input_name}.read.1.fq.gz {input_name}.read.2.fq.gz -o {output_name}.sam ",
    )
    # file path hacking
    runDocker(
        "hlaminer",
        f"HLAminer.pl -h {index_ref} -p {index}/hla_nom_p.txt -i 90 -q 10 -s 100 "
        f" -a {output_name}.sam -l /../{output_name} ",
    )
    return output_name


def hlaminerReadResult(input_name: str) -> str:
    """
    Read hlaminer csv result into our hla_result format

    Its format:
    ```
    HLA-C
        Prediction #1 - C*06
            C*06:103,1676.94,1.63e-28,277.9
            C*06:02P,202.00,8.18e-04,30.9
        Prediction #2 - C*01
            C*01:99,1010.00,3.65e-16,154.4
            C*01:02P,994.00,3.65e-16,154.4

    HLA-E
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    alleles = []
    for i in open(f"{input_name}.csv"):
        if "Prediction" in i:
            alleles.append(i.split("-")[1].strip())
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def seq2hlaBuild(folder: str, db_hla: str = "origin") -> str:
    """https://github.com/TRON-Bioinformatics/seq2HLA"""
    output_name = f"{folder}/origin"
    if Path(output_name).exists():
        return output_name
    runDocker(
        "seq2hla", f"cp -r /usr/local/share/seq2hla-2.2-2/references {output_name}"
    )
    return output_name


def seq2hlaRun(input_name: str, index: str) -> str:
    """https://github.com/TRON-Bioinformatics/seq2HLA"""
    output_name = input_name + ".seq2hla_" + name2Single(index)
    if Path(f"{output_name}.-ClassI.HLAgenotype4digits").exists():
        return output_name
    runDocker(
        "seq2hla",
        f"seq2HLA -r {output_name}. -p {getThreads()}"
        f" -1 {input_name}.read.1.fq.gz -2 {input_name}.read.2.fq.gz",
        opts=f" -v {Path(index).absolute()}:/usr/local/share/seq2hla-2.2-2/references ",
    )
    return output_name


def seq2hlaReadResult(input_name: str) -> str:
    """
    Read seq2HLA 4digits txt into our hla_result format

    Its format:
    ```
    A   A*36:01 0.520369    A*11:02'    0.0
    B   B*08:01 0.0007307464    B*56:05'    0.0
    C   C*07:01'    0.0 C*07:01 NA
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    txt1 = pd.read_csv(f"{input_name}.-ClassI.HLAgenotype4digits", sep="\t")
    txt2 = pd.read_csv(f"{input_name}.-ClassII.HLAgenotype4digits", sep="\t")
    txt = pd.concat([txt1, txt2])
    allele1 = txt[["Allele 1", "Confidence"]]
    allele2 = txt[["Allele 2", "Confidence.1"]]
    allele1 = allele1.rename(columns={"Allele 1": "Allele"})
    allele2 = allele2.rename(
        columns={"Allele 2": "Allele", "Confidence.1": "Confidence"}
    )
    df = pd.concat([allele1, allele2])
    df = df[df["Confidence"] > 0.1]
    df1 = allelesToTable(df["Allele"], default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlahdDownload(folder: str = "hlahd", version: str = "1.5.0") -> str:
    """https://www.genome.med.kyoto-u.ac.jp/HLA-HD/download-request/"""
    if not Path(f"{folder}/hlahd.{version}.tar.gz").exists():
        raise ValueError(
            f"Require manually download requests and place the file here: {folder}/hlahd.{version}.tar.gz"
        )
    if Path(f"{folder}/hlahd").exists():
        return folder
    runShell(f"tar -vxf {folder}/hlahd.{version}.tar.gz -C {folder}")
    runShell(f"mv {folder}/hlahd.{version} {folder}/hlahd")
    # install.sh include index building So i separate it
    runShell(
        f"grep -h g++ {folder}/hlahd/install.sh {folder}/hlahd/update.dictionary.sh > {folder}/hlahd/compile.sh"
    )
    return folder


def hlahdBuild(folder: str, db_hla: str = "origin") -> str:
    """https://www.genome.med.kyoto-u.ac.jp/HLA-HD/"""
    if db_hla == "origin":
        output_name = f"{folder}/hla_3320"
        if Path(output_name).exists():
            return output_name
        runShell(f"cp -r {folder}/hlahd/dictionary {output_name}")
        runShell(
            f"cp -r {folder}/hlahd/HLA_gene.split.3.32.0.txt {output_name}/HLA_gene.split.txt"
        )
        runDocker("hlahd", "bash bw_build.sh", chdir=output_name)
    else:
        output_name = f"{folder}/{Path(db_hla).name}"
        if Path(output_name).exists():
            return output_name
        runShell(f"mkdir -p {output_name}")
        runShell(f"cp {db_hla}/hla.dat {output_name}/")
        # see update.dictionary.sh
        runDocker(
            "hlahd", f"create_fasta_from_dat {output_name}/hla.dat {output_name}/ 150"
        )
        # why this script not exists
        runDocker("hlahd", f"bash create_dir.sh", chdir=output_name)
        runDocker("hlahd", f"bash move_file.sh", chdir=output_name)
        runDocker("hlahd", f"bash bw_build.sh", chdir=output_name)
        runShell(
            f"cp -r {folder}/hlahd/HLA_gene.split.3.32.0.txt {output_name}/HLA_gene.split.txt"
        )
    return output_name


def hlahdRun(input_name: str, index: str) -> str:
    """https://www.genome.med.kyoto-u.ac.jp/HLA-HD/"""
    output_name = input_name + ".hlahd_" + name2Single(index)
    if Path(output_name + "/data/result").exists():
        return output_name
    runShell(f"mkdir -p {output_name}/data")
    runDocker(
        "hlahd",
        f"hlahd.sh -t {getThreads()} {input_name}.read.1.fq.gz {input_name}.read.2.fq.gz "
        f"{index}/HLA_gene.split.txt {index}/ data {output_name}",
    )
    return output_name


def hlahdReadResult(input_name: str) -> str:
    """
    Read hlahd final result into our hla_result format

    Its format:
    ```
    A       Not typed       Not typed
    B       HLA-B*07:02:01  HLA-B*82:02
    C       Not typed       Not typed
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    txt = pd.read_csv(
        f"{input_name}/data/result/data_final.result.txt",
        sep="\t",
        names=["gene", "allele1", "allele2"],
        usecols=[0, 1, 2],
    )
    txt = txt.dropna()
    # print(txt)
    alleles = list([*txt["allele1"], *txt["allele2"]])
    alleles = [i.replace("HLA-", "") for i in alleles if i != "Not typed" and i != "-"]
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def phlatBuild(folder: str, db_hla: str = "origin") -> str:
    """Their website is gone. So download from docker"""
    output_name = f"{folder}/hla_3080"  # described in paper (2014)
    if Path(output_name).exists():
        return output_name
    # resources/preload1.pk is also a index file, but whatever.
    runDocker("phlat_download", f"cp -r /usr/bin/phlat-release/b2folder/ {output_name}")
    return output_name


def phlatRun(input_name: str, index: str, index_bam: str) -> str:
    """Input hs38DH and using run.b38.sh instead of directly call PHLAT.py"""
    output_name = input_name + ".phlat_" + name2Single(index)
    if Path(output_name + "/data_HLA.sum").exists():
        return output_name
    runShell(f"mkdir -p {output_name}/data")
    runDocker(
        "phlat",
        f"bash /usr/bin/run.b38.sh --data-dir {output_name} --rs-dir {output_name} "
        f" --ref-fasta {index_bam} --bam {input_name}.bam --tag data --index-dir {index} ",
    )
    return output_name


def phlatReadResult(input_name: str) -> str:
    """
    Read phlat result into our hla_result format

    Its format:
    ```
    Locus   Allele1 Allele2 LLtot   pval1   pval2
    HLA_B   B*07:02:01      B*82:02 -2356.54        1.6e-04 3.9e-04
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    df = pd.read_csv(f"{input_name}/data_HLA.sum", sep="\t")
    # print(txt)
    alleles = list([*df["Allele1"], *df["Allele2"]])
    alleles = [i.replace("HLA-", "") for i in alleles if i != "Not typed" and i != "-"]
    df1 = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlaforestBuild(folder: str, db_hla: str = "origin") -> str:
    """https://github.com/FNaveed786/hlaforest"""
    if db_hla == "origin":
        output_name = f"{folder}/hla_3100"  # described in paper
    else:
        output_name = f"{folder}/{Path(db_hla).name}"
    if Path(output_name + "/config.sh").exists():
        return output_name

    if not Path(f"{folder}/hlaforest").exists():
        runShell(
            f"git clone https://github.com/FNaveed786/hlaforest.git {folder}/hlaforest"
        )
        runShell(f"cd {folder}/hlaforest && git checkout 75edb46")

    if db_hla == "origin":
        runShell(f"mv {folder}/hlaforest/reference {output_name}")
    else:
        # Find the building step by comparing it's nospace and filtered fasta
        runShell(f"mkdir -p {output_name}")
        with (
            open(f"{db_hla}/hla_nuc.fasta") as f_in,
            open(f"{output_name}/hla_nuc_nospace_filtered.fasta", "w") as f_out,
        ):
            skip = False
            for line in f_in:
                if line.startswith(">"):
                    if line.split(" ")[1].endswith("N"):
                        skip = True
                    else:
                        skip = False
                        f_out.write("_".join(line.split(" ")))
                elif not skip:
                    f_out.write(line)
        runDocker(
            "hlaforest",
            f"bowtie-build {output_name}/hla_nuc_nospace_filtered.fasta "
            f"             {output_name}/hla_nuc_nospace_filtered",
        )
    with (
        open(f"{folder}/hlaforest/scripts/config.sh") as f_in,
        open(f"{output_name}/config.sh", "w") as f_out,
    ):
        for i in f_in.readlines():
            if i.startswith("SCRIPT_PATH"):
                f_out.write("SCRIPT_PATH=/usr/local/bin\n")
            elif i.startswith("BOWTIE_INDEX"):
                f_out.write(
                    f"BOWTIE_INDEX=/app/{output_name}/hla_nuc_nospace_filtered\n"
                )
            elif i.startswith("NUM_THREADS"):
                f_out.write("NUM_THREADS=$NUM_THREADS\n")
            else:
                f_out.write(i)
    return output_name


def hlaforestRun(input_name: str, index: str) -> str:
    """https://github.com/FNaveed786/hlaforest"""
    output_name = input_name + ".hlaforest_" + name2Single(index)
    if Path(output_name + "/haplotypes.txt").exists():
        return output_name
    runDocker(
        "hlaforest",
        f"CallHaplotypesPE.sh {output_name} {input_name}.read.1.fq.gz {input_name}.read.2.fq.gz",
        opts=(
            f" -v {Path(index).absolute()}/config.sh:/root/hla/hlaforest/scripts/config.sh:ro "
            f" -e NUM_THREADS={getThreads()} "
        ),
    )
    return output_name


def hlaforestReadResult(input_name: str) -> str:
    """
    Read hlaforest haplotype result into our hla_result format

    Its format:
    ```
    ROOT:A  ROOT:A:01       ROOT:A:01:03
    ROOT:A  ROOT:A:11       ROOT:A:11:01    ROOT:A:11:01:18
    ROOT:B  ROOT:B:08       ROOT:B:08:20
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    alleles = []
    with open(f"{input_name}/haplotypes.txt") as f:
        for i in f:
            alleles.append(i.split("\t")[-1])
    alleles = [
        i.replace("ROOT:", "").replace(":", "*", 1).strip()
        for i in alleles
        if i.strip()
    ]
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlaprofilerBuild(folder: str, db_hla: str = "origin") -> str:
    """https://expressionanalysis.github.io/HLAProfiler/"""
    output_name = f"{folder}/hla_3240"  # described in paper
    if Path(output_name).exists():
        return output_name

    runShell(
        f"wget https://github.com/ExpressionAnalysis/HLAProfiler/archive/refs/tags/v1.0.0-db_only.tar.gz -P {folder}"
    )
    runShell(f"tar -vxf {folder}/v1.0.0-db_only.tar.gz -C {folder}")
    runShell(f"mv {folder}/HLAProfiler-1.0.0-db_only/hla_database {output_name}")
    runShell(
        f"wget https://github.com/ExpressionAnalysis/HLAProfiler/releases/download/v1.0.0-db_only/database.idx -P {output_name}"
    )
    runShell(
        f"wget https://github.com/ExpressionAnalysis/HLAProfiler/releases/download/v1.0.0-db_only/database.kdb -P {output_name}"
    )
    return output_name


def hlaprofilerRun(input_name: str, index: str) -> str:
    """https://github.com/ExpressionAnalysis/HLAProfiler"""
    output_name = input_name + ".hlaprofiler_" + name2Single(index)
    name = Path(input_name).name + ".read.1.fq.gz"
    if Path(output_name + f"/{name}/{name}.HLATypes.txt").exists():
        return output_name
    runDocker(
        "hlaprofiler",
        f"HLAProfiler.pl predict -kp /usr/local/share/kraken-ea-0.10.5ea.3-3/ -intermediate_files"
        f" -intermediate_files -threads {getThreads()}"
        f" -database_name . -database_dir {index} -r {index}/data/reference/hla.ref.merged.fa"
        f" -output_dir {output_name} -l {output_name}.log"
        f" -fastq1 {input_name}.read.1.fq.gz "
        f" -fastq2 {input_name}.read.2.fq.gz ",
    )
    return output_name


def hlaprofilerReadResult(input_name: str) -> str:
    """
    Read hlaprofilerhaplotype result into our hla_result format

    Its format:
    ```
    Allele1_Accession       Allele2_Accession       Allele1 Allele2 Proportion_reads...
    HLA     HLA     A       A       -       -       -       -       -       -       Not enough reads to make call
    HLA11935        HLA05253        B*82:02:02      B*07:02:20      0.292442325075655 ...
    HLA02519        HLA11935        B*07:49N        B*82:02:02      0.267367979772837 ...
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    name = Path(list(Path(input_name).iterdir())[0]).name
    txt = pd.read_csv(f"{input_name}/{name}/{name}.HLATypes.txt", sep="\t")
    txt = txt[txt["Final_score"] != "-"]
    txt["Final_score"] = txt["Final_score"].astype(float)
    txt["gene"] = txt["Allele1"].str.split("*", expand=True)[0]
    df = (
        txt.sort_values(["gene", "Final_score"], ascending=[True, False])
        .groupby("gene")
        .head(1)
    )
    df1 = allelesToTable([*df["Allele1"], *df["Allele2"]], default_gene=["A", "B", "C"])
    df1["name"] = input_name
    df1.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def athlatesDownload(folder: str) -> str:
    """
    Website (https://www.broadinstitute.org/viral-genomics/athlates)
    has binary execute file, but require download submittion.
    And script building and whole pipeline are written in
    https://github.com/cliu32/athlates
    """
    if Path(f"{folder}/athlates").exists():
        return folder

    if not Path(f"{folder}/athlates.zip").exists():
        raise ValueError(f"Please download to {folder}/athlates.zip")

    runShell(f"unzip {folder}/athlates.zip -d {folder}")
    runShell(f"mv {folder}/Athlates_2014_04_26 {folder}/athlates")
    return folder


def athlatesBuild(folder: str, db_hla: str = "origin") -> str:
    """https://www.broadinstitute.org/viral-genomics/athlates"""
    if db_hla == "origin":
        output_name = f"{folder}/hla_3090"  # written in A_nuc.txt
    else:
        output_name = f"{folder}/{Path(db_hla).name}"
    if Path(output_name + "/ref").exists():
        return output_name

    if db_hla == "origin":
        runShell(f"cp -r {folder}/athlates/db {output_name}")
        runShell(f"ln -s DRB_nuc.txt {output_name}/msa/DRB1_nuc.txt")
        runShell(f"ln -s DQB_nuc.txt {output_name}/msa/DQB1_nuc.txt")
        runDocker(
            "athlates",
            f"novoindex {output_name}/ref/hla.clean.nix {output_name}/ref/hla.clean.fasta",
        )
    else:
        runShell(f"mkdir -p {output_name}")
        runShell(f"cp -r {db_hla}/alignments {output_name}/msa")
        runShell(f"cp {db_hla}/hla_gen.fasta {output_name}/")
        runDocker(
            "athlates",
            f"perl /usr/local/bin/hla_ref_clean.pl "
            f"-ref {output_name}/hla_gen.fasta  -oprefix {output_name}/hla",
        )
        runShell(f"mkdir -p {output_name}/bed {output_name}/ref")
        runShell(f"mv {output_name}/*.bed {output_name}/bed")
        runShell(f"mv {output_name}/*.fasta {output_name}/ref")
        runDocker(
            "athlates",
            f"novoindex {output_name}/ref/hla.clean.nix {output_name}/ref/hla.clean.fasta",
        )
    return output_name


def athlatesPreRun(input_name: str, index: str) -> str:
    """Run read mapping, and separate each athlates gene task"""
    output_name = input_name + ".athlates_" + name2Single(index) + ".novoalign"
    output_call = output_name + ".call.{}"
    if Path(output_name + ".bam").exists():
        return output_call
    runShell(f"gunzip -f --keep {input_name}.read.1.fq.gz")
    runShell(f"gunzip -f --keep {input_name}.read.2.fq.gz")
    runDocker(
        "athlates",
        "sh -c 'novoalign -o SAM -r Random -i PE 100-1400 -H -t 30 -n 100 "
        f"  -d {index}/ref/hla.clean.nix"
        f"  -f {input_name}.read.1.fq {input_name}.read.2.fq"
        f"  >  {output_name}.sam'",
    )
    bamSort(output_name, sam=True)

    # same as "perl /usr/local/bin/runAthlates_usr_update.pl "
    # but using `samtools sort -n `
    # f"  -ibam {output_name2}.bam "
    # f"  -bed {index}/bed/hla.{gene}.bed -nbed {index}/bed/hla.non-{gene}.bed "
    # f"  -msa {index}/msa/{gene}_nuc.txt "
    # f"  -odir {parent} -oprf {name}.{gene}"
    for gene in ["A", "B", "C", "DRB1", "DQB1"]:
        output_name1 = output_call.format(gene) + ".target"
        output_name2 = output_call.format(gene) + ".nontarget"

        with open(f"{output_call.format(gene)}.sh", "w") as f:
            f.write("set -x\n")
            f.write("set -o pipefail\n")
            #  f.write("set -e\n")  # enable this will break most of samples
            f.write(
                f"samtools view -h -F 4 -L {index}/bed/hla.{gene}.bed {output_name}.bam -o {output_name1}.bam"
            )
            f.write("\n")
            f.write(
                f"samtools sort -n {output_name1}.bam -o {output_name1}.sortname.bam"
            )
            f.write("\n")
            f.write(
                f"samtools view -h -F 4 -L {index}/bed/hla.non-{gene}.bed {output_name}.bam -o {output_name2}.bam"
            )
            f.write("\n")
            f.write(
                f"samtools sort -n {output_name2}.bam -o {output_name2}.sortname.bam"
            )
            f.write("\n")
            f.write(
                f"typing -bam {output_name1}.sortname.bam -exlbam {output_name2}.sortname.bam "
                f"       -msa {index}/msa/{gene}_nuc.txt -hd 2 "
                f"       -o {output_call.format(gene)}"
            )
    return output_call


def athlatesRun(input_name: str) -> str:
    """https://github.com/cliu32/athlates"""
    output_name = input_name
    if Path(input_name + ".typing.txt").exists():
        return output_name
    try:
        runDocker("athlates", f"bash {input_name}.sh")
    except subprocess.CalledProcessError:
        pass
    return output_name


def athlatesReadResult(input_name: str) -> str:
    """
    Read athlates result into our hla_result format

    Its format:
    ```
     ------------ Inferred Allelic Pairs -------------

     A*01:01:01:02N      A*11:126        0
     A*01:01:01:01       A*11:126        0
    ```
    """
    output_name = input_name.replace(".{}", "_merge") + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    alleles = []
    for f in glob(f"{input_name.format('*')}.typing.txt"):
        skip = True
        for i in open(f):
            if "Inferred" in i:
                skip = False
                continue
            if skip or not i.strip():
                continue
            # first infer candiated
            alleles.extend(i.strip().split()[:2])
            break

    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def hlassignBuild(folder: str, db_hla: str = "origin") -> str:
    """
    HLAssign v2 doesn't provided database
    So I use first version: which is 3.12.0

    Also, in alignments/B_nuc.txt after 3.34.0,
    B*13:123Q introduced and break the format of txt,
    I temporary remove the tailing sequences.
    """
    if db_hla == "origin":
        output_name = f"{folder}/hla_3120"  # folder name
    else:
        output_name = f"{folder}/{Path(db_hla).name}"

    if Path(f"{output_name}/DPB1_nuc.hlp.txt").exists():
        return output_name

    if db_hla == "origin":
        runShell(
            f"wget https://www.ikmb.uni-kiel.de/sites/default/files/downloads/hlassign/pilot_source_code.zip -P {folder}"
        )
        runShell(f"unzip {folder}/pilot_source_code.zip -d {folder}")
        runShell(f"cp -r {folder}/pilot_source_code/IMGTv312 {output_name}")
        runShell(f"cp {folder}/pilot_source_code/ambiguity_table.txt {output_name}")
        runShell(f"cp {folder}/pilot_source_code/incomplete_alleles.txt {output_name}")
        runShell(f"ln -s DRB_gen.txt {output_name}/DRB1_gen.txt")
        runShell(f"ln -s DPA_nuc.txt {output_name}/DPA1_nuc.txt")
        runShell(f"ln -s DPB_nuc.txt {output_name}/DPB1_nuc.txt")
        runShell(f"ln -s DQA_nuc.txt {output_name}/DQA1_nuc.txt")
        runShell(f"ln -s DQB_nuc.txt {output_name}/DQB1_nuc.txt")
    else:
        runShell(f"cp -r {db_hla}/alignments {output_name}")
        runShell(f"ln -s DRB_nuc.txt {output_name}/DRB1_nuc.txt")

        # modifiy nuc
        for gene in ["A", "B", "C", "DRB1", "DPA1", "DPB1", "DQA1", "DQB1"]:
            txt = open(f"{output_name}/{gene}_nuc.txt").read()
            txt_arr = txt.split("cDNA")
            txt_arr_new = [txt_arr[0]]
            num_line = txt_arr[0].count("\n")
            for i in txt_arr[1:]:
                if i.count("\n") < num_line - 10:
                    continue
                num_line = i.count("\n")
                txt_arr_new.append(i)
            open(f"{output_name}/{gene}_nuc.txt", "w").write("cDNA".join(txt_arr_new))

        # I don't know how the generated
        runShell(f"cp {folder}/pilot_source_code/incomplete_alleles.txt {output_name}")
        # I think this is G group
        # runShell(f"cp {folder}/pilot_source_code/ambiguity_table.txt {output_name}")
        arrs = []
        for line in open("/home/linnil1//HLA/hla_3490/wmda/hla_nom_g.txt"):
            if line.startswith("#"):
                continue
            gene, group, name = line.strip().split(";")
            if not name:
                name = group.split("/")[0]
            arrs.append([gene + name, *[gene + i for i in group.split("/")]])
        df = pd.DataFrame(arrs)  # it will become same length
        df.to_csv(
            f"{output_name}/ambiguity_table.txt", index=False, sep="\t", header=False
        )

    # ANY gene
    for gene in ["A", "B", "C", "DRB1", "DPA1", "DPB1", "DQA1", "DQB1"]:
        runDocker(
            "hlassign",
            f"gethlamultialignexonwise {output_name}/{gene}_gen.txt "
            f"  > {output_name}/{gene}_gen.hlp.txt",
        )
        runDocker(
            "hlassign",
            f"gethlamultialignexonwise {output_name}/{gene}_nuc.txt "
            f"  > {output_name}/{gene}_nuc.hlp.txt",
        )
    return output_name


def hlassignPreRun(input_name: str, index: str) -> str:
    """
    Split each gene for typing by hlassign.
    https://www.ikmb.uni-kiel.de/resources/download-tools/software/hlassign
    """
    output_name_base = input_name + ".hlassign_" + name2Single(index)
    output_name_pattern = output_name_base + ".gene.{}"

    if Path(output_name_pattern.format("C") + ".sh").exists():
        return output_name_pattern

    runShell(f"zcat {input_name}.read.*.fq.gz > {output_name_base}.fq")
    # get read length by read one line in readfile
    read_f = open(output_name_base + ".fq")
    next(read_f)  # at least one read exists
    read = next(read_f)
    read_length = len(read.strip())

    # for gene in ["A", "B", "C"]:
    for gene in ["B"]:
        output_name = output_name_pattern.format(gene)
        if Path(output_name + ".sh").exists():
            continue
        with open(output_name + ".sh", "w") as f:
            f.write("set -x\n")
            f.write("set -o pipefail\n")
            f.write(
                f"prefilterhlareads --Ref {index}/{gene}_gen.hlp.txt "
                f"  --read_length {read_length} --min_alignment 30 "
                f"  --fastq {output_name_base}.fq "
                f"  --on-target {output_name}.target.fq --off-target {output_name}.offtarget.fq"
            )
            output_name = output_name + ".target"
            runShell(f"mkdir -p {output_name}/{Path(output_name).parent}")
            # Usage said cd {output_name}
            # I don't like it
            parent_str = "../" * len(Path(output_name).parent.parents)
            runShell(f"ln -s ./{parent_str} {output_name}/{output_name}")
            f.write("\n")
            f.write(f"sed -n '1~4s/^@/>/p;2~4p' {output_name}.fq > {output_name}.fa")
            f.write("\n")
            f.write(
                f"nextcallhla --cDNA {index}/{gene}_nuc.txt --gDNA {index}/{gene}_gen.txt"
                f"  --bias-alleles {index}/incomplete_alleles.txt "
                f"  --amb-table {index}/ambiguity_table.txt "
                f"  --read-length {read_length}"
                f"  --reads {output_name}.fa "
                f"  --out {output_name}.result.txt "
                f"  --scratch {output_name} --outPic {output_name} "
            )
    return output_name_pattern


def hlassignRun(input_name: str) -> str:
    """Each task takes very long time"""
    if Path(input_name + ".result.txt").exists():
        return input_name
    runDocker("hlassign", f"bash {input_name}.sh")
    return input_name


def hlassignReadResult(input_name: str) -> str:
    """
    Read hlassign result into our hla_result format

    Its format:
    ```
    #[Result]
    B*82:02 B*82:02 -   1
    ```
    """
    output_name = input_name.replace(".{}", "_merge") + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    alleles = []
    for f in glob(f"{input_name.format('*')}.result.txt"):
        skip = True
        for i in open(f):
            if "Result" in i:
                skip = False
                continue
            if skip or not i.strip():
                continue
            # first infer candiated
            alleles.extend(i.strip().split()[:2])
            break

    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def stcseqBuild(folder: str, db_hla: str = "origin") -> str:
    """https://ngdc.cncb.ac.cn/biocode/tools/BT007068/releases/v1.0"""
    if db_hla == "origin":
        output_name = f"{folder}/hla_3260"  # written in paper

    if Path(output_name + "/exon-70bp.fastq").exists():
        return output_name

    runShell(
        f"wget https://ngdc.cncb.ac.cn/biocode/tools/7068/releases/v1.0/file/download?filename=STC-Seq.tar.gz -O {folder}/STC-Seq.tar.gz"
    )
    runShell(f"tar -vxf {folder}/STC-Seq.tar.gz -C {folder}")
    runShell(f"mv {folder}/STC-Seq/data {output_name}")
    runDocker(
        "stcseq",
        f"getArtifical_reads {output_name}/hla.exon.fasta {output_name}/exon-70bp.fastq",
    )
    return output_name


def stcseqRun(input_name: str, index: str) -> str:
    """ "https://ngdc.cncb.ac.cn/biocode/tools/BT007068/manual"""
    output_folder = input_name + ".stcseq_" + name2Single(index)
    output_name = output_folder + "/data"
    if Path(output_name + "_report.txt").exists():
        return output_name
    runShell(f"mkdir -p {output_folder}")
    runShell(f"zcat {input_name}.read.*.fq.gz > {output_name}.fq")
    runDocker("stcseq", f"step_1.sh {output_name} {index}")
    runDocker("stcseq", f"step_2.sh {output_name} {index}")
    runDocker(
        "stcseq",
        f"Rscript /usr/local/bin/Initial_screening.R "
        f"{output_name}_position-array.txt 7 15 {output_name}_second_retain_allele.txt",
    )
    runDocker("stcseq", f"step_3.sh {output_name} {index}")
    return output_name


def stcseqReadResult(input_name: str) -> str:
    """
    Read STC-seq into our hla_result format

    Its format:
    ```
    ##Valid allele pair
    Locus   Genotype    Alternative_genotype    Quality
    B   B*07:02:01,B*82:02:01   null    PASS
    ```
    """
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name
    txt = pd.read_csv(f"{input_name}_report.txt", sep="\t", comment="#")
    alleles_pair = txt["Genotype"]
    alleles = [j for i in alleles_pair for j in i.split(",")]
    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def renameResult(input_name: str, sample_name: str) -> str:
    """rename xx.*.hla_result.csv into xx.hla_result.{method}.csv"""
    suffix = input_name.split(sample_name)[1]
    suffix = suffix.split(".hla_result")[0]
    suffix = suffix[1:].replace(".", "_").replace("/", "_")
    output_name = f"{sample_name}.hla_result.{suffix}"
    print(output_name)
    # "../"  for data/
    runShell(f"ln -sf ../{input_name}.tsv {output_name}.tsv")
    return sample_name + ".hla_result.{}"


def mergeResult(input_name: str) -> str:
    """Merge hla_result data and print it"""
    files = glob(input_name.replace(".{}", ".*.tsv"))
    print(files)
    df = pd.concat(pd.read_csv(i, sep="\t") for i in files)
    print(input_name.replace(".{}", "_merge"))
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.max_colwidth",
        None,
        "display.width",
        200,
    ):
        print(df)
        print(df.sort_values("gene"))

    output_name = input_name.replace(".{}", "_merge")
    df.to_csv(f"{output_name}.csv", index=False)
    return output_name


def readArgument() -> argparse.Namespace:
    """Read command line arguments"""
    parser = argparse.ArgumentParser(
        prog="hla_collections",
        description="Run all HLA tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "sample_name",
        help="samples name. Use {} to indicate wildcard.",
    )
    parser.add_argument("--version", default="origin", help="IMGT-HLA version")
    parser.add_argument("--thread", default=4, help="Number of threads")
    parser.add_argument(
        "--tools",
        default=[
            # "arcashla",  # for RNA, easy to fail
            # "athlates",  # for RNA, easy to fail
            "bwakit",
            "graphtyper",
            "hisat",
            # 'hlaforest',  # for RNA, easy to fail
            # 'hlahd',  # require download requests, disable as default
            "hlala",
            "hlaminer",
            "hlaprofiler",
            "hlascan",
            "kourami",
            "optitype",
            # "phlat",  # require download requests
            "polysolver",
            "seq2hla",
            "soaphla",
            "vbseq",
            "xhla",
        ],
        nargs="+",
        help="HLA tools to execute",
    )
    args = parser.parse_args()
    # print(args)
    return args


if __name__ == "__main__":
    args = readArgument()
    setThreads(args.thread)
    if args.sample_name == "example":
        samples_fq = downloadSample()
    else:
        samples_fq = str(args.sample_name)
    version = args.version
    if version == "origin":
        db_hla = "origin"
    else:
        db_hla = downloadHLA("", version=version)

    # hs37
    folder_bwa = createFolder("", "bwakit")
    if "vbseq" in args.tools or "polysolver" in args.tools or "soaphla" in args.tools:
        index_hs37 = downloadRef(folder_bwa, name="hs37")
        index_hs37 = bwaIndex(index_hs37)
        samples_hs37 = bwaRun(samples_fq, index_hs37)

    if "polysolver" in args.tools:
        samples = polysolverRun(samples_hs37, "")
        samples = polysolverReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "soaphla" in args.tools:
        folder = createFolder("", "soaphla")
        folder = checkAndBuildImage(folder)
        index = soaphlaBuild(folder, db_hla)
        samples = soaphlaRun(samples_hs37, index)
        samples = soaphlaReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "vbseq" in args.tools:
        folder = createFolder("", "vbseq")
        folder = checkAndBuildImage(folder)
        index = vbseqBuild(folder, db_hla)
        samples = vbseqPreprocess(samples_hs37)
        # samples_new = samples  # directly mapped on HLA
        samples = vbseqRun(samples, index)
        samples = vbseqReadResult(samples)
        samples = renameResult(samples, samples_fq)

    # hs38
    if "graphtyper" in args.tools or "xhla" in args.tools:
        index_hs38 = downloadRef(folder_bwa, name="hs38")
        index_hs38 = bwaIndex(index_hs38)
        samples_hs38 = bwaRun(samples_fq, index_hs38)

    if "graphtyper" in args.tools:
        folder = createFolder("", "graphtyper")
        folder = checkAndBuildImage(folder)
        index = graphtyperBuild(folder, db_hla, index_hs38)
        samples = graphtyperPreRun(samples_hs38, index)
        samples = globAndRun(graphtyperRun, samples)
        samples = graphtyperReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "xhla" in args.tools:
        folder = createFolder("", "xhla")
        folder = checkAndBuildImage(folder)
        index = xhlaBuild("xhla", db_hla)
        samples = xhlaRun(samples_hs38, index)
        samples = xhlaReadResult(samples)
        samples = renameResult(samples, samples_fq)

    # hs38DH
    if (
        "bwakit" in args.tools
        or "arcashla" in args.tools
        or "phlat" in args.tools
        or "hlala" in args.tools
        or "hlascan" in args.tools
        or "kourami" in args.tools
    ):
        index_bwa = downloadRef(folder_bwa, name="hs38DH")
        index_bwa = bwaIndex(index_bwa)
        samples = bwakitRun(samples_fq, index_bwa)
        samples_new = bwakitReadResult(samples)
        samples_new = renameResult(samples_new, samples_fq)
        samples = addSuffix(samples, ".aln")
        samples_hs38dh = bamSort(samples)

    if "arcashla" in args.tools:
        db_hla_for_arcas = db_hla
        if version == "origin":
            db_hla_for_arcas = downloadHLA("", version="3.24.0")
        folder = createFolder("", "arcas")
        index = arcasBuild(folder, db_hla_for_arcas)
        # I found the decoy in hs38DH wll not be extracted
        samples = arcasPreprocess(samples_hs38dh)
        # samples = samples_fq  # directly mapped
        samples = arcasRun(samples, index)
        samples = arcasReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlala" in args.tools:
        folder = createFolder("", "hlala")
        index = hlalaBuild(folder)
        samples = addUnmap(samples_hs38dh)
        samples = hlalaRun(samples, index)
        samples = hlalaReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlascan" in args.tools:
        folder = createFolder("", "hlascan")
        folder = checkAndBuildImage(folder)
        index = hlascanBuild(folder)
        # directly mapped on HLA from fastq
        # samples = hlascanPreRun(samples_fq, index)
        # samples = globAndRun(hlascanRun, samples)
        # samples = hlascanReadResult(samples)
        # samples = renameResult(samples, samples_fq)
        samples = hlascanPreRun(samples_hs38dh, index)
        samples = globAndRun(hlascanRun, samples)
        samples = hlascanReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "kourami" in args.tools:
        folder = createFolder("", "kourami")
        folder = kouramiDownload(folder)
        folder = checkAndBuildImage(folder, "kourami_preprocess")
        index = kouramiBuild(folder, db_hla)
        samples = kouramiPreprocess(samples_hs38dh, index, folder)
        samples = kouramiRun(samples, index, folder)
        samples = kouramiReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "phlat" in args.tools:
        folder = createFolder("", "phlat")
        index = phlatBuild(folder, db_hla)
        samples = phlatRun(samples_hs38dh, index, index_bwa)
        samples = phlatReadResult(samples)
        samples = renameResult(samples, samples_fq)

    # fastq
    if "hisat" in args.tools:
        # maximum version 3.43.0
        db_hla_for_hisat = db_hla
        if version != "origin" and tuple(map(int, version.split("."))) > (3, 43, 0):
            db_hla_for_hisat = downloadHLA("", version="3.43.0")
        folder = createFolder("", "hisat")
        folder = hisatDownload(folder)
        folder = checkAndBuildImage(folder)
        index = hisatBuild(folder, db_hla_for_hisat)
        samples = hisatRun(samples_fq, index)
        samples = hisatReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlassign" in args.tools:
        folder = createFolder("", "hlassign")
        folder = checkAndBuildImage(folder)
        index = hlassignBuild(folder, db_hla)
        samples = hlassignPreRun(samples_fq, index)
        samples = globAndRun(hlassignRun, samples)
        samples = hlassignReadResult(samples)

    # Below tools mapped all the reads on HLA sequences
    # without any filtering (e.g. remove unrelated read)
    # it may cause some problems
    if "athlates" in args.tools:
        db_hla_for_athlates = db_hla
        if version != "origin" and tuple(map(int, version.split("."))) > (3, 31, 0):
            db_hla_for_athlates = downloadHLA("", version="3.31.0")
        folder = createFolder("", "athlates")
        folder = athlatesDownload(folder)
        folder = checkAndBuildImage(folder)
        index = athlatesBuild(folder, db_hla_for_athlates)
        samples = athlatesPreRun(samples_fq, index)
        samples = globAndRun(athlatesRun, samples)
        samples = athlatesReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlaforest" in args.tools:
        folder = createFolder("", "hlaforest")
        folder = checkAndBuildImage(folder)
        index = hlaforestBuild(folder, db_hla)
        samples = hlaforestRun(samples_fq, index)
        samples = hlaforestReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlahd" in args.tools:
        folder = createFolder("", "hlahd")
        folder = hlahdDownload(folder)
        folder = checkAndBuildImage(folder)
        # require about 40G
        index = hlahdBuild(folder, db_hla)
        # require about 40G
        samples = hlahdRun(samples_fq, index)
        samples = hlahdReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlaminer" in args.tools:
        folder = createFolder("", "hlaminer")
        folder = checkAndBuildImage(folder)
        index = hlaminerBuild(folder, db_hla)
        samples = hlaminerRun(samples_fq, index)
        samples = hlaminerReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "hlaprofiler" in args.tools:
        folder = createFolder("", "hlaprofiler")
        index = hlaprofilerBuild(folder, db_hla)
        samples = hlaprofilerRun(samples_fq, index)
        samples = hlaprofilerReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "optitype" in args.tools:
        folder = createFolder("", "optitype")
        index = optitypeBuild(folder)
        samples = optitypeRun(samples_fq, index)
        samples = optitypeReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "seq2hla" in args.tools:
        folder = createFolder("", "seq2hla")
        folder = checkAndBuildImage(folder)
        index = seq2hlaBuild(folder, db_hla)
        samples = seq2hlaRun(samples_fq, index)
        samples = seq2hlaReadResult(samples)
        samples = renameResult(samples, samples_fq)

    if "stcseq" in args.tools:
        folder = createFolder("", "stcseq")
        folder = checkAndBuildImage(folder)
        index = stcseqBuild(folder, db_hla)
        samples = stcseqRun(samples_fq, index)
        samples = stcseqReadResult(samples)
        samples = renameResult(samples, samples_fq)

    mergeResult(samples_fq + ".hla_result.{}")
