def contig_with_most_reads(file: str):
    max_reads = 0
    num_reads = 0
    last_contig = False
    with open(file) as arq:
        line = arq.readline()
        while not line.startswith('******************* Contig'):
            line = arq.readline()
        actual_contig = line.strip()
        actual_contig = actual_contig.strip("*")
        actual_contig = actual_contig.strip(" ")
        while line:
            if (line.startswith('******************* Contig') and num_reads != 0) or last_contig:
                print(actual_contig)
                print(num_reads)
                if num_reads > max_reads:
                    max_reads = num_reads
                    max_contig = actual_contig
                num_reads = 0
                actual_contig = line.strip()
                actual_contig = actual_contig.strip("*")
                actual_contig = actual_contig.strip(" ")
                if last_contig:
                    break
            elif not line.startswith('*'):
                num_reads += 1
            line = arq.readline()
            if line.startswith('DETAILED DISPLAY OF CONTIGS'):
                last_contig = True
                num_reads -= 1
    return (max_contig, max_reads)

print(contig_with_most_reads("consenso"))