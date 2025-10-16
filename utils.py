from alinhamento import Contig, Alinhamento

def conversor_alinhamento(cap3_out_file: str, previous_contig: str) -> Alinhamento:
    max_reads = 0
    num_reads = 0
    max_contig = None
    rev = None
    # Checa se é o último contig
    last_contig = False
    with open(cap3_out_file) as file:
        line = file.readline()
        # Lê até achar o começo da descrição dos contigs
        while line and not line.startswith('******************* Contig'):
            line = file.readline()
        # actual_contig pega o nome do contig atual
        actual_contig = line.strip()
        actual_contig = actual_contig.strip("*")
        actual_contig = actual_contig.strip(" ")
        # Variável que diz se o contig possui o contig anterior
        contains = False
        while line:
            # Checa se já chegou no próximo contig
            if (line.startswith('******************* Contig') and num_reads != 0) or last_contig:
                if num_reads > max_reads and contains:
                    max_reads = num_reads
                    max_contig = actual_contig
                    rev = rev_temp
                num_reads = 0
                actual_contig = line.strip()
                actual_contig = actual_contig.strip("*")
                actual_contig = actual_contig.strip(" ")
                # Reseta variável dizendo que contém o contig anterior
                contains = False
                if last_contig:
                    break
            elif not line.startswith('*'):
                line = line.strip(" ")
                if (previous_contig is not None and line.startswith(previous_contig)) or (previous_contig is None):
                    contains = True
                    if line[len(previous_contig):len(previous_contig)+1] == '+':
                        print("A linha está no reverso com essa cara: ")
                        print(line)
                        rev_temp = False
                    else:
                        print("A linha NÃO está no reverso com essa cara: ")
                        print(line)
                        rev_temp = True
                num_reads += 1
            line = file.readline()
            if line.startswith('DETAILED DISPLAY OF CONTIGS'):
                last_contig = True
                num_reads -= 1
    if max_contig is not None:
        return (max_contig, rev)
    else:
        print("Nessa iteração não foi possível achar contig com o contig anterior")
        return self.contig_with_most_reads("temp/consenso", previous_contig=None)