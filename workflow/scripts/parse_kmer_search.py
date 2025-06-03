import os
import sys
import csv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from loguru import logger
import pandas as pd

epilog = "Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)"


class IndexResultParser:
    """
    Class provides functions to parse table formats
    returned by kmer indexing tools and map to sample names
    if required.
    """

    def __init__(
        self,
        search_results: str,
        method: str,
        raptor_sample_mapping: str = None,
        kmer_ratio: float = 0.7,
    ) -> None:

        self.search_results = search_results
        self.method = method
        # Parameters specific for prediction tools
        self.raptor_sample_mapping = raptor_sample_mapping
        if self.method == "raptor":
            assert (
                self.raptor_sample_mapping is not None
                and self.raptor_sample_mapping != ""
            ), "Parsing Raptor results requires a sample/index mapping file"
        self.kmer_ratio = kmer_ratio

    def parse_results(self) -> pd.DataFrame:
        """
        Choose parsing method for specified method
        """
        result = {}
        match self.method:
            case "kmindex":
                logger.info("-> Parsing KMINDEX index query results...")
                result = self._parse_kmindex()
            case "raptor":
                logger.info("-> Parsing RAPTOR index query results...")
                result = self._parse_raptor()
            case _:
                logger.error(f"Tool {self.method} is unknown. Cannot parse results")
        return result

    @staticmethod
    def write_result(
        results: pd.DataFrame, out_file: str, out_type: str = "tsv"
    ) -> None:
        """
        Write parsed results into tabular format to be processed by user
        :param results: A dictionary qith cts as key and sample hits as value
        :param out_file: Output file
        :return:
        """
        assert out_type in [
            "tsv",
            "parquet",
        ], "Supported output types are parquet and tsv"
        with open(out_file, "wb") as file_handle:
            if out_type == "tsv":
                results.to_csv(file_handle, sep="\t")
            else:
                results.to_parquet(file_handle, compression="snappy", index=True)

    def _parse_kmindex(self) -> pd.DataFrame:
        """
        Parse tabular output format of kmindex
        """
        results = {}
        with open(self.search_results) as file_handle:
            reader = csv.DictReader(file_handle, delimiter="\t")
            for line in reader:
                cts_id = line["samples"].split(":")[1]
                results[cts_id] = {}
                for sample, prediction in line.items():
                    if sample == "samples":
                        continue
                    prediction = round(float(prediction), 2)
                    results[cts_id][sample] = prediction
        results = pd.DataFrame.from_dict(
            results, orient="index", dtype="Sparse[float64]"
        )
        logger.info(f"-> Parsed {len(results)} query sequences")
        return results

    def _parse_raptor(self) -> pd.DataFrame:
        """
        Parse raptor search results

        :return: Result mapping
        """
        dataset_mapping = {}
        sample_name_mapping = {}
        results = {}

        with open(self.raptor_sample_mapping, "r") as file_handle:
            logger.info(
                "-> Reading minimiser2sample mapping file to match raptor bin ids to sample names"
            )
            reader = csv.DictReader(file_handle, delimiter="\t")
            for row in reader:
                sample_name_mapping[row["minimiser_id"]] = row["sample_name"]

        with open(self.search_results) as file_handle:
            for line in file_handle:
                elements = line.rstrip().split("\t")
                # Skip config section returned in raptor output file -> Starts with '##'
                if line.startswith("##"):
                    continue
                # Get key/sample mapping from header
                elif line.startswith("#"):
                    # Stop collecting samples when header of results section start
                    if elements[0] != "#QUERY_NAME":
                        dataset_mapping[int(elements[0][1:])] = elements[1].rstrip()
                # Parse index hits -> Each query has a number of bins assigned Q1   3,4,5,6
                else:
                    elements = line.rstrip().split("\t")
                    cts_id = elements[0]
                    results[cts_id] = {}
                    # By default not detected in any sample
                    not_detected_samples = set(dataset_mapping.keys())
                    # If at least one bin is reported by raptor
                    if len(elements) > 1:
                        detected_samples = set()
                        # Iterate over samples bins
                        for this_sample in elements[1].split(","):
                            # Save detected bins and directly translate into sample identifier
                            detected_samples.add(int(this_sample))
                            results[cts_id][
                                sample_name_mapping[dataset_mapping[int(this_sample)]]
                            ] = self.kmer_ratio
                            # results[cts_id].append(sample_name_mapping[dataset_mapping[int(this_sample)]])
                        # Update not detected samples by removing bins with at least >= k-mer fraction
                        not_detected_samples = not_detected_samples - detected_samples
                    # Write annotation status for samples without a hit
                    for this_sample in not_detected_samples:
                        results[cts_id][
                            sample_name_mapping[dataset_mapping[int(this_sample)]]
                        ] = None
        results = pd.DataFrame.from_dict(
            results, orient="index", dtype="Sparse[float64, nan]"
        )
        logger.info(f"-> Parsed {len(results)} query sequences")
        return results


def main():
    parser = ArgumentParser(
        description="parseKmer - Parse search results of Raptor and Kmindex",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        "--search-results",
        dest="search_results",
        action="store",
        help="Output file of k-mer index search",
        required=True,
    )
    parser.add_argument(
        "--output",
        dest="output",
        action="store",
        help="Output table listing detection status in individual RNA-seq samples",
        required=True,
    )
    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        required=True,
        help="Input file is from Raptor/kmindex",
    )
    parser.add_argument(
        "--raptor-sample-mapping",
        dest="raptor_sample_mapping",
        action="store",
        required=False,
        help="Mapping of Raptor HIBF bins to real samples identifiers",
    )
    parser.add_argument(
        "--kmer-ratio",
        dest="kmer_ratio",
        action="store",
        required=False,
        default=0.7,
        type=float,
        help="K-mer ratio used for search",
    )

    args = parser.parse_args()
    logger.info(f"-> Parsing k-mer query file {args.search_results}")
    if args.kmer_ratio > 1.0 or args.kmer_ratio < 0.0:
        raise ValueError("k-mer ratio needs to be [0,1)")

    parser = IndexResultParser(
        search_results=args.search_results,
        method=args.method,
        raptor_sample_mapping=args.raptor_sample_mapping,
        kmer_ratio=args.kmer_ratio,
    )
    parsed_results = parser.parse_results()
    parser.write_result(parsed_results, args.output)


if __name__ == "__main__":
    main()
