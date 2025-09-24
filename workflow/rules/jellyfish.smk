rule jellyfish_query:
    input:
        ,
    output:
        "search.fasta",
    log:
        "{prefix}.jf.log",
    threads: 2
    shell:
       "jellyfish query "
