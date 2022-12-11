from glob import glob
from pathlib import Path
from itertools import chain
from collections import defaultdict
from typing import Iterable
import re
import uuid
import json
import gzip
import argparse
import subprocess
import pandas as pd


resources: dict[str, int] = {  # per sample
    "threads": 4,
    "memory": 7,  # unit: G
}

images = {
    "samtools": "quay.io/biocontainers/samtools:1.15.1--h1170115_0",
    "hisat": "localhost/linnil1/hisat2:1.3.3",
    # "hisat": "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4",
    "bwa": "quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
    "bwakit": "quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1",
    "java": "docker.io/library/openjdk:11-jdk",
    "ubuntu": "docker.io/library/ubuntu:22.04",
    "vbseq": "localhost/linnil1/vbseq:20181122",
    "hlala": "quay.io/biocontainers/hla-la:1.0.3--hd03093a_0",
    # "hlala": " localhost/linnil1/hlala:alpine",  # linnil1-built for fun
    "kourami_preprocess": "localhost/linnil1/kourami_preprocess",
    "xhla": "localhost/linnil1/xhla",
    "graphtyper": "localhost/linnil1/graphtyper:2.7.5",
}


def getThreads() -> int:
    return resources["threads"]


def setThreads(threads: int) -> None:
    global resources
    resources["threads"] = threads


def name2Single(name: str) -> str:
    """Trun everything special character ('./') into '_'"""
    return name.replace("/", "_").replace(".", "_")


def runDocker(image: str, cmd: str, opts: str = "") -> subprocess.CompletedProcess[str]:
    """run docker container"""
    image = images.get(image, image)
    random_name = str(uuid.uuid4()).split("-", 1)[0]
    cmd_all = (
        f"podman run -it --rm -u root -w /app -v $PWD:/app "
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


def buildImage(dockerfile: str, image: str) -> None:
    """build docker image"""
    runShell(f"podman build . -f {dockerfile} -t {image}")


def checkImage(image: str) -> bool:
    """check image exists"""
    try:
        runShell(
            f"sh -c 'if [ ! $(podman image ls {images[image]} -q) ]; then exit 1; fi'"
        )
        return True
    except subprocess.CalledProcessError:
        return False


def addSuffix(input_name: str, suffix: str) -> str:
    return input_name + suffix


def allelesToTable(
    alleles: Iterable[str], default_gene: list[str] = []
) -> pd.DataFrame:
    """Turn list of alleles into our hla_result format"""
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
        f"samtools view -@{getThreads()} -n {name}.download.cram -o {name}.bam",
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
    Path(folder).mkdir(exist_ok=True)
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


def kouramiBuildImage(folder: str) -> str:
    """
    see https://github.com/Kingsford-Group/kourami/blob/master/preprocessing.md
    and script/alignAndExtract_hs38DH
    docker build -f kourami_preprocess.dockerfile
    """
    if checkImage("kourami_preprocess"):
        return folder
    buildImage("kourami_preprocess.dockerfile", images["kourami_preprocess"])
    return folder


def kouramiDownload(folder: str = "kourami") -> str:
    """https://github.com/Kingsford-Group/kourami"""
    if Path(folder + "/resources/hs38NoAltDH.fa").exists():
        return folder
    runShell(
        "wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9.6/kourami-0.9.6_bin.zip"
    )
    runShell("unzip kourami-0.9.6_bin.zip")
    runShell(f"mv kourami-0.9.6 {folder}")
    runShell(f"mv kourami-0.9.6_bin.zip {folder}")
    # runDocker("bwakit", f"bash {folder}/scripts/download_grch38.sh hs38DH")
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
        kourami_db = f"{folder}/hla_3240"
    else:
        kourami_db = f"{folder}/{Path(db_hla).name}"

    if Path(f"{kourami_db}/All_FINAL_with_Decoy.fa.gz").exists():
        return kourami_db

    if db_hla == "origin":
        # scripts/download_panel.sh
        # runDocker("bwa", f"bash {folder}/scripts/download_panel.sh")
        runShell(
            f"wget https://github.com/Kingsford-Group/kourami/releases/download/v0.9/kouramiDB_3.24.0.tar.gz -P {folder}"
        )
        runShell(f"tar -vxf {folder}/kouramiDB_3.24.0.tar.gz -C {folder}")
        runShell(f"mv {folder}/db {kourami_db}")
    else:
        # Kourami's bug?
        # runDocker("kourami", f"java -cp {folder}/build/Kourami.jar FormatIMGT {db_hla}/alignments/ . {kourami_db}")
        runShell(f"cp -r {db_hla}/alignments {folder}/tmp_alignments")
        runShell(
            f"wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/3240/alignments/Y_nuc.txt -O {folder}/tmp_alignments/Y_nuc.txt"
        )
        runShell(
            f"wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/3240/alignments/Y_gen.txt -O {folder}/tmp_alignments/Y_gen.txt"
        )
        runDocker(
            "java",
            f"java -cp {folder}/build/Kourami.jar FormatIMGT {folder}/tmp_alignments/ . {kourami_db}",
        )
        # runShell(f"rm -r {folder}/tmp_alignments")
        runShell(
            f"cat {kourami_db}/*.merged.fa {folder}/resources/HLA_decoys.fa | gzip > {kourami_db}/All_FINAL_with_Decoy.fa.gz"
        )
        runShell(f"cp {db_hla}/wmda/hla_nom_g.txt {kourami_db}")
    bwaIndex(kourami_db + "/All_FINAL_with_Decoy.fa.gz")
    return kourami_db


def kouramiPreprocess(input_name: str, index: str, kourami_folder: str = "") -> str:
    """https://github.com/Kingsford-Group/kourami"""
    output_name = input_name + ".kourami_preprocess"
    if Path(f"{output_name}_extract_2.fq.gz").exists():
        return output_name
    if not kourami_folder:
        kourami_folder = index + "/.."
    runDocker(
        "kourami_preprocess",
        f"bash {kourami_folder}/scripts/alignAndExtract_hs38DH.sh "
        f"-d {index} -r {kourami_folder}/resources/hs38NoAltDH.fa "
        f"{output_name} "
        f"{input_name}.bam ",
    )
    return output_name


def kouramiRun(input_name: str, index: str, kourami_folder: str) -> str:
    """https://github.com/Kingsford-Group/kourami"""
    panel_name = input_name + ".panel_" + name2Single(index)
    output_name = panel_name + ".call"
    if Path(f"{output_name}.result").exists():
        return output_name
    # Because we use same preprocessed file
    # the last step of alignAndExtract_hs38DH generate extracted fastq
    runDocker(
        "kourami_preprocess",
        f"bwa mem -t {getThreads()} {index}/All_FINAL_with_Decoy.fa.gz "
        f"{input_name}_extract_1.fq.gz {input_name}_extract_2.fq.gz -o {panel_name}.sam ",
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
    Path(folder).mkdir(exist_ok=True)
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
        "hisat", f"sh -c 'hisat2-inspect {folder}/grch38/genome > {folder}/genome.fa'"
    )
    runDocker("samtools", f"samtools faidx {folder}/genome.fa")
    return folder


def hisatBuildImage(folder: str) -> str:
    """docker build . -f hisat.dockerfile"""
    if checkImage("hisat"):
        return folder
    runShell(
        f"git clone --recurse-submodules https://github.com/DaehwanKimLab/hisat-genotype {folder}/hisat-genotype"
    )
    buildImage("hisat.dockerfile", images["hisat"])
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

    if Path(f"{output_name}/hisatgenotype_db/HLA/msf").exists():
        return output_name

    # link basic things
    files_genome = glob(folder + "/genotype_genome*")
    for f in files_genome:
        runShell(f"ln -s ../{Path(f).name} {output_name}/")
    runShell(f"    ln -s ../grch38         {output_name}/")
    runShell(f"    ln -s ../genome.fa      {output_name}/")
    runShell(f"    ln -s ../genome.fa.fai  {output_name}/")

    # download HLA related things
    if db_hla == "origin":
        runShell(
            f"git clone https://github.com/DaehwanKimLab/hisatgenotype_db {output_name}/hisatgenotype_db"
        )
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
    return output_name


def hisatRun(input_name: str, index: str) -> str:
    """https://daehwankimlab.github.io/hisat-genotype/manual/"""
    output_name = input_name + ".hisat_" + name2Single(index)
    if len(glob(f"{output_name}/*.report")) > 0:
        return output_name
    json.dump({"sanity_check": False}, open(f"{index}/settings.json", "w"))
    runDocker(
        "hisat",
        f"hisatgenotype -z {index} --threads {getThreads()} "
        f"--base hla --out-dir {output_name} "
        f"--keep-alignment -v --keep-extract "
        f"-1 {input_name}.read.1.fq.gz "
        f"-2 {input_name}.read.2.fq.gz ",
        opts=f"-v $PWD/{index}/settings.json:/opt/hisat-genotype/devel/settings.json",
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


def hlascanDownload(folder: str = "hlascan") -> str:
    """https://github.com/SyntekabioTools/HLAscan"""
    if Path(f"{folder}/hlascan").exists():
        return folder
    Path(folder).mkdir(exist_ok=True)
    runShell(
        "wget https://github.com/SyntekabioTools/HLAscan/releases/download/v2.0.0/dataset.zip -P {folder}"
    )
    runShell(f"unzip {folder}/dataset.zip -d {folder}")
    runShell(
        f"wget https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4 -P {folder}"
    )
    runShell(f"ln -s hla_scan_r_v2.1.4 {folder}/hlascan")
    return folder


def hlascanRun(input_name: str, index: str) -> str:
    """
    https://github.com/SyntekabioTools/HLAscan
    I'm not sure what input should given for hlascan:
    * selected Fastq
    * hs38 with alt (hs38DH)?
    * hs37 or hs37d5
    """
    output_name = input_name + ".hlascan_" + name2Single(index)
    if Path(f"{output_name}.HLA-C.txt").exists():
        return output_name + ".{}"
    for gene in ["HLA-A", "HLA-B", "HLA-C"]:
        try:
            if Path(input_name + ".bam").exists():
                if "hs37" in input_name:
                    version = "37"
                elif "hs38" in input_name:
                    version = "38"
                else:
                    raise ValueError("reference version cannot determine")
                runDocker(
                    "ubuntu",
                    f"./{index}/hlascan -g {gene} -t {getThreads()} -d {index}/db/HLA-ALL.IMGT "
                    f"-b {input_name}.bam -v {version} > {output_name}.{gene}.txt",
                )
            else:
                runDocker(
                    "ubuntu",
                    f"./{index}/hlascan -g {gene} -t {getThreads()} -d {index}/db/HLA-ALL.IMGT "
                    f"-l {input_name}.read.1.fq.gz -r {input_name}.read.2.fq.gz > {output_name}.{gene}.txt",
                )
        except subprocess.CalledProcessError as e:
            continue
    return output_name + ".{}"


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
        for i in open(input_name.replace(".{}", f".{gene}.txt")):
            if "[Type" in i:
                alleles.append(f"{gene_name}*{i.split()[2]}")

    df = allelesToTable(alleles, default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def vbseqDownload(folder: str = "vbseq") -> str:
    """https://nagasakilab.csml.org/hla/"""
    Path(folder).mkdir(exist_ok=True)
    return folder


def vbseqBuildImage(folder: str) -> str:
    """docker build -f vbseq.dockerfile"""
    if checkImage("vbseq"):
        return folder
    runShell(f"wget https://nagasakilab.csml.org/hla/HLAVBSeq.jar -P {folder}")
    runShell(f"wget https://nagasakilab.csml.org/hla/parse_result.pl -P {folder}")
    buildImage("vbseq.dockerfile", images["vbseq"])
    return folder


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


def hlalaDownload(folder: str = "hlala") -> str:
    """https://github.com/DiltheyLab/HLA-LA"""
    db = f"{folder}/PRG_MHC_GRCh38_withIMGT"
    if Path(f"{db}/serializedGRAPH").exists():
        return db
    Path(folder).mkdir(exist_ok=True)
    # buildImage("hlala.dockerfile", images['hlala'])  # linnil1-built
    runShell(
        f"wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz -P {folder}"
    )
    runShell(f"tar -vxf {folder}/PRG_MHC_GRCh38_withIMGT.tar.gz -C {folder}")
    runDocker(
        "hlala",
        "sh -c 'export PATH=/usr/local/opt/hla-la/bin:$PATH && "
        f"HLA-LA --workingDir ./ --action prepareGraph --PRG_graph_dir {db}'",
    )
    return db


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
        opts=f" -v $PWD/{index}:/usr/local/opt/hla-la/graphs/",
        # opts=f" -v $PWD/{index}:/usr/local/bin/HLA-LA/graphs/",  # linnil1-built
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


def xhlaDownload(folder: str = "xhla") -> str:
    """https://github.com/humanlongevity/HLA"""
    if Path(f"{folder}/HLA").exists():
        return folder
    Path(folder).mkdir(exist_ok=True)
    runShell(f"git clone https://github.com/humanlongevity/HLA {folder}/HLA")
    return folder


def xhlaBuildImage(folder: str) -> str:
    """docker build -f xhla.dockerfile"""
    if checkImage("xhla"):
        return folder
    buildImage("xhla.dockerfile", images["xhla"])
    return folder


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
        runShell(f"mkdir -p {output_name1}/raw")
        runShell(f"cp -r {folder}/HLA/data {output_name1}/data")
        runShell(f"cp -r {db_hla}/alignments {output_name1}/raw")
        for i in ["nuc", "exon", "dna", "temp", "align"]:
            runShell(f"mkdir -p {output_name1}/data/{i}")
        runDocker("xhla", f"sh -c 'cd {output_name1}/data/ && bash script/batch.sh'")
        runDocker(
            "xhla",
            f"sh -c 'cd {output_name1}/data/ && diamond makedb --in hla.faa --db hla.dmnd'",
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
        opts=f"-v $PWD/{index}:/opt/data ",
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
    print(df_full)
    df = allelesToTable(df_full["full"], default_gene=["A", "B", "C"])
    df["name"] = input_name
    df.to_csv(output_name + ".tsv", index=False, sep="\t")
    return output_name


def graphtyperDownload(folder: str = "graphtyper") -> str:
    """https://github.com/DecodeGenetics/graphtyper"""
    if Path(f"{folder}").exists():
        return folder
    Path(folder).mkdir(exist_ok=True)
    return folder


def graphtyperBuildImage(folder: str) -> str:
    """docker build -f graphtyper.dockerfile"""
    if checkImage("graphtyper"):
        return folder
    buildImage("graphtyper.dockerfile", images["graphtyper"])
    return folder


def graphtyperBuild(
    folder: str, db_hla: str = "origin", index_hs38: str = "bwakit/hs38.fa"
) -> str:
    """
    IMGT version 3.23.0 written in their paper.
    I'm not sure if database is updated.

    Ignore db_hla, It doesn't provided scripts for updating index
    """
    output_name = folder + "/hla_3230"
    name = Path(index_hs38).name
    if Path(f"{output_name}/{name}.fai").exists():
        return output_name
    runShell(f"mkdir -p {output_name}")
    runShell(
        f"wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.5/HLA_data.tar.gz -P {output_name}"
    )
    runShell(f"tar -vxf {output_name}/HLA_data.tar.gz -C {output_name}")
    runDocker("samtools", f"samtools faidx {index_hs38}")
    runShell(f"ln -s ../../{index_hs38} {output_name}/{name}")
    runShell(f"ln -s ../../{index_hs38}.fai {output_name}/{name}.fai")
    return output_name


def graphtyperRun(input_name: str, index: str) -> str:
    """https://github.com/DecodeGenetics/graphtyper/wiki/HLA-genotyping"""
    output_name = input_name + ".graphtyper_" + name2Single(index)
    if Path(output_name).exists():
        return output_name
    id = Path(output_name).name
    df = pd.read_csv(
        f"{index}/regions.tsv", sep="\t", names=["chrom", "start", "end", "gene"]
    )
    for i in df.itertuples():
        open(f"{output_name}.gene.{i.gene}.txt", "w").write(f"{input_name}.bam")
        runDocker(
            "graphtyper",
            f"graphtyper genotype_hla {index}/hs38.fa --verbose --threads={getThreads()} "
            f"{index}/{i.gene}.vcf.gz --region={i.chrom}:{i.start}-{i.end} "
            f"--sam={input_name}.bam --output={output_name}",
        )
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
    output_name = input_name + ".hla_result"
    if Path(f"{output_name}.tsv").exists():
        return output_name

    alleles: list[str] = []
    for vcf_file in glob(input_name + "/chr6/*.vcf.gz"):
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


def renameResult(input_name: str, sample_name: str) -> str:
    """rename xx.*.hla_result.csv into xx.hla_result.{method}.csv"""
    suffix = input_name.split(sample_name)[1]
    suffix = suffix.split(".hla_result")[0]
    suffix = suffix[1:].replace(".", "_")
    output_name = f"{sample_name}.hla_result.{suffix}"
    print(output_name)
    # "../"  for data/
    runShell(f"ln -sf ../{input_name}.tsv {output_name}.tsv")
    return sample_name + ".hla_result.{}"


def mergeResult(input_name: str) -> str:
    """Merge hla_result data and print it"""
    files = glob(input_name.replace("{}", "*.tsv"))
    print(files)
    df = pd.concat(pd.read_csv(i, sep="\t") for i in files)
    print(input_name.replace("{}", "_merge"))
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
    parser.add_argument(
        "--tools",
        default=[
            "bwakit",
            "graphtyper",
            "hisat",
            "hlala",
            "hlascan",
            "kourami",
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
    if args.sample_name == "example":
        samples = downloadSample()
    else:
        samples = str(args.sample_name)
    version = args.version
    if version == "origin":
        db_hla = "origin"
    else:
        db_hla = downloadHLA("", version=version)

    # bwa
    if "xhla" in args.tools or "graphtyper" in args.tools:
        index_hs38 = downloadRef("bwakit", name="hs38")
        index_hs38 = bwaIndex(index_hs38)
        samples_hs38 = bwaRun(samples, index_hs38)
    if "vbseq" in args.tools:
        index_hs37 = downloadRef("bwakit", name="hs37")
        index_hs37 = bwaIndex(index_hs37)
        samples_hs37 = bwaRun(samples, index_hs37)

    # bwakit
    if (
        "bwakit" in args.tools
        or "kourami" in args.tools
        or "hlascan" in args.tools
        or "hlala" in args.tools
    ):
        index_bwakit = downloadRef("bwakit", name="hs38DH")
        index_bwakit = bwaIndex(index_bwakit)
        samples_new = bwakitRun(samples, index_bwakit)
        samples_new_1 = bwakitReadResult(samples_new)
        renameResult(samples_new_1, samples)
        samples_new = addSuffix(samples_new, ".aln")
        samples_hs38dh = bamSort(samples_new)

    # bwakit -> kourami
    if "kourami" in args.tools:
        folder_kourami = kouramiDownload("kourami")
        folder_kourami = kouramiBuildImage(folder_kourami)
        index_kourami = kouramiBuild(folder_kourami, db_hla)
        samples_new = kouramiPreprocess(samples_hs38dh, index_kourami, folder_kourami)
        samples_new = kouramiRun(samples_new, index_kourami, folder_kourami)
        samples_new = kouramiReadResult(samples_new)
        renameResult(samples_new, samples)

    # bwakit -> hlala
    if "hlala" in args.tools:
        index_hlala = hlalaDownload("hlala")
        samples_new = addUnmap(samples_hs38dh)
        samples_new = hlalaRun(samples_new, index_hlala)
        samples_new = hlalaReadResult(samples_new)
        renameResult(samples_new, samples)

    # bwakit -> hlascan
    if "hlascan" in args.tools:
        index_hlascan = hlascanDownload("hlascan")
        samples_new = hlascanRun(samples, index_hlascan)
        samples_new = hlascanReadResult(samples_new)
        renameResult(samples_new, samples)
        samples_new = hlascanRun(samples_hs38dh, index_hlascan)
        samples_new = hlascanReadResult(samples_new)
        renameResult(samples_new, samples)

    if "hisat" in args.tools:
        # maximum version 3.43.0
        db_hla_for_hisat = db_hla
        if version != "origin" and tuple(map(int, version.split("."))) > (3, 43, 0):
            db_hla_for_hisat = downloadHLA("", version="3.43.0")
        folder_hisat = hisatDownload("hisat")
        folder_hisat = hisatBuildImage(folder_hisat)
        index_hisat = hisatBuild(folder_hisat, db_hla_for_hisat)
        samples_new = hisatRun(samples, index_hisat)
        samples_new = hisatReadResult(samples_new)
        renameResult(samples_new, samples)

    if "vbseq" in args.tools:
        folder_vbseq = vbseqDownload("vbseq")
        folder_vbseq = vbseqBuildImage(folder_vbseq)
        index_vbseq = vbseqBuild(folder_vbseq, db_hla)
        samples_new = vbseqPreprocess(samples_hs37)
        # samples_new = samples  # directly mapped on HLA
        samples_new = vbseqRun(samples_new, index_vbseq)
        samples_new = vbseqReadResult(samples_new)
        renameResult(samples_new, samples)

    # hs38 -> xhla
    if "xhla" in args.tools:
        folder_xhla = xhlaDownload("xhla")
        folder_xhla = xhlaBuildImage(folder_xhla)
        index_xhla = xhlaBuild(folder_xhla, db_hla)
        samples_new = xhlaRun(samples_hs38, index_xhla)
        samples_new = xhlaReadResult(samples_new)
        renameResult(samples_new, samples)

    if "graphtyper" in args.tools:
        folder_graphtyper = graphtyperDownload("graphtyper")
        folder_graphtyper = graphtyperBuildImage(folder_graphtyper)
        index_graphtyper = graphtyperBuild(folder_graphtyper, db_hla, index_hs38)
        samples_new = graphtyperRun(samples_hs38, index_graphtyper)
        samples_new = graphtyperReadResult(samples_new)
        renameResult(samples_new, samples)

    mergeResult(samples + ".hla_result.{}")
