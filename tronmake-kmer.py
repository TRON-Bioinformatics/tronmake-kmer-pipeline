#!/usr/bin/env python3

import os
import sys
import argparse
import pathlib
import tempfile
import yaml
import subprocess
from logzero import logger


__version__ = "2.0.0"
__pipeline__ = pathlib.Path(__file__).parent / 'workflow' / 'Snakefile'

epilog = "Copyright (c) 2023 TRON gGmbH (See LICENSE for licensing details)"

def execute_cmd(cmd, working_dir = "."):
    """This function runs a command into a subprocess."""
    logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
    p = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, cwd = working_dir, shell=False)
    if p.returncode != 0:
        logger.error(p.stderr)
    return p.returncode


def indexing_pipeline(args):
    wf_config = {}
    wf_config["modus"] = {"query": False,
                          "indexing": True}
    wf_config["indexing"] = {
        "samples": args.samples,
        "kmer_size": args.kmer,
        "cutoff": args.cutoff,
        "method": args.method,
        "fpr": args.fpr,
        "quantitative_index": args.quantitative
    }
    with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=args.workdir) as temp_config:
        yaml.dump(wf_config, temp_config)
        temp_config.close()
        cmd = ['snakemake',
               '--snakefile', str(__pipeline__),
               '--local-cores', str(args.jobs),
               '--jobs', str(args.jobs),
               '--configfile', str(temp_config.name),
               '--use-conda',
               '--directory', str(args.workdir),
               '--rerun-triggers', 'mtime']
        if args.slurm:
            cmd.extend(['--executor', 'slurm'])
        returncode = execute_cmd(cmd)

        if returncode != 0:
            logger.error("-> Command \"{}\" returned non-zero exit status".format(cmd))
            sys.exit(1)
        else:
            logger.info("-> Pipeline finished")

def query_pipeline(args):
    wf_config = {}
    wf_config["modus"] = {"query": True,
                          "indexing": False}
    wf_config["query"] = {
        "index": args.index_manifest,
        "query_fasta": args.fasta ,
        "kmer_ratio": args.detection_ratio,
        "findere_z": args.findere
    }
    with tempfile.NamedTemporaryFile(mode="w", delete=False, dir=args.workdir) as temp_config:
        yaml.dump(wf_config, temp_config)
        temp_config.close()
        cmd = ['snakemake',
               '--snakefile', str(__pipeline__),
               '--local-cores', str(args.jobs),
               '--jobs', str(args.jobs),
               '--configfile', str(temp_config.name),
               '--use-conda',
               '--directory', str(args.workdir),
               '--rerun-triggers', 'mtime']
        if args.slurm:
            cmd.extend(['--executor', 'slurm'])
        returncode = execute_cmd(cmd)

        if returncode != 0:
            logger.error("-> Command \"{}\" returned non-zero exit status".format(cmd))
            sys.exit(1)
        else:
            logger.info("-> Pipeline finished")

def add_index_parser_args(parser):
    parser.add_argument(
        "--samples",
        dest="samples",
        help="Sample sheet (tsv)",
        required=True
    )
    parser.add_argument(
        "--kmer",
        dest="kmer",
        help="kmer size used for index construction",
        default=21,
        type=int
    )
    parser.add_argument(
        "--cutoff",
        dest="cutoff",
        help="Cutoff to define solid and weak k-mers. Only solid k-mers are included in index",
        default=2,
        type=int
    )
    parser.add_argument(
        "--method",
        dest="method",
        help="indexing method to use",
        default="raptor"
    )
    parser.add_argument(
        "--fpr",
        dest="fpr",
        help="FPR for kmindex/Raptor indexing",
        default=0.05, type=float
    )
    parser.add_argument(
        "--workdir",
        dest="workdir",
        help="Work directory for pipeline execution",
        default=pathlib.Path(__file__).parent
    )
    parser.add_argument(
        "--quantitative",
        dest="quantitative",
        help="Create quantitative kmindex index. EXPERIMENTAL",
        action='store_true'
    )
    parser.add_argument(
        "--jobs",
        dest="jobs",
        help="Number of local CPUs or number of jobs for slurm submission",
        default=16
    )
    parser.set_defaults(func=indexing_pipeline)


def add_query_parser_args(parser):
    parser.add_argument(
        "--index-manifest",
        dest="index_manifest",
        help="k-mer index to query",
        required=True,
    )
    parser.add_argument(
        "--fasta",
        dest="fasta",
        help="Query sequences in fasta format",
        required=True,
    )
    parser.add_argument(
        "--detection-ratio",
        dest="detection_ratio",
        help="k-mer detection ratio / Percentage of shared k-mers for presence/absence detection",
        default=0.7,
        type=float
    )
    parser.add_argument(
        "--findere",
        dest="findere",
        help="Z-value of findere algorithm. Only relevant when method = kmindex",
        default=2,
        type=int
    )
    parser.add_argument(
        "--workdir",
        dest="workdir",
        help="Work directory for pipeline execution",
        default=pathlib.Path(__file__).parent
    )
    parser.add_argument(
        "--jobs",
        dest="jobs",
        help="Number of local CPUs or number of jobs for slurm submission",
        default=16
    )
    parser.set_defaults(func=query_pipeline)



def tronmake_cli():
    parser = argparse.ArgumentParser(
        description="TronMake k-mer pipeline v{}".format(__version__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--slurm",
        dest="slurm",
        help="Execute snakemake pipeline with slurm support",
        action="store_true"
    )

    subparsers = parser.add_subparsers(description="Commands")

    indexing_parser = subparsers.add_parser(
        "index", 
        description="Runs the index construction pipeline",
        epilog=epilog, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_index_parser_args(indexing_parser)

    query_parser = subparsers.add_parser(
        "query",
        description="Runs the index query pipeline",
        epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_query_parser_args(query_parser)


    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError as e:
        logger.exception(e)
        parser.parse_args(["--help"])

if __name__ == "__main__":
    tronmake_cli()
