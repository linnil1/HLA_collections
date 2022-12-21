"""
File: Resource part of pipeline.py
Author: linnil1
Description: Modify this if you want to do more
The resouces include:
* Memory
* Threads
* Images name
* Folder location
* Code to running shell script
* Code to running container
* Code to building container
* The code to check image existance

"""
from pathlib import Path
import uuid
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
    "novoalign": "quay.io/biocontainers/novoalign:3.09.04--h82c745c_3",
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


def runDocker(
    image: str,
    cmd: str,
    mounts: list[tuple[str, str]] = [],
    envs: list[tuple[str, str]] = [],
    chdir: str = "",
) -> subprocess.CompletedProcess[str]:
    """run docker container"""
    image = images.get(image, image)
    random_name = str(uuid.uuid4()).split("-", 1)[0]
    opts = ""
    for src, dst in mounts:
        opts += f" -v {Path(src).absolute()}:{dst}"
    for src, val in envs:
        opts += f" -e {src}={val}"
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
