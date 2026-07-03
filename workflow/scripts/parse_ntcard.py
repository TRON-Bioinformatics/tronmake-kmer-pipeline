#!/usr/bin/env python3


def parse_ntcard(ntcard_histo, parsed_histo, sample_name):
    with open(ntcard_histo, "r") as file_handle, open(
        parsed_histo, "w"
    ) as write_handle:
        f0, f1 = 0, 0
        for line in file_handle:
            elements = line.rstrip().split("\t")
            if elements[0] == "F0":
                f0 = int(elements[1])
            elif elements[0] == "1":
                f1 = int(elements[1])
            else:
                continue
        result = f"{sample_name}\t{f0}\t{f1}\t{f0- f1}\n"
        write_handle.write(result)


parse_ntcard(
    snakemake.input["histo"],
    snakemake.output["parsed_histo"],
    snakemake.wildcards["sample"],
)
