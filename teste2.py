def extract_contig_with_coverage(cap3_output, contig_name, coverage_threshold=3):

    # Ler o arquivo de saída do CAP3
    with open(cap3_output, 'r') as file:
        line = file.readline()
        cont = 0
        while cont < 2:
            if line.startswith(f'******************* {contig_name} ********************'):
                cont += 1
            line = file.readline()
        line = file.readline()

        contig_final = ""
        
        # Variável para guardar sequências acima do consenso
        seqs = []

        achou_comeco = False

        while not (line.startswith('******************* Contig')):
            
            if not(line.startswith("*") or line.startswith(" ") or len(line.strip()) == 0 or line.startswith("consensus")):
                seqs.append(line.strip("\n"))

            if line.startswith("consensus"):
                line = line.strip()
                contig_seq = line.split(" ")
                contig_seq = contig_seq[len(contig_seq)-1]
                comeco_analise = len(line) - len(contig_seq)
                fim_analise = len(line)
                contig_final += line[comeco_analise:fim_analise]


                # Calcular a cobertura para cada posição
                coverage = [0] * len(contig_seq)
                for seq in seqs:
                    j = 0
                    for i in range(comeco_analise, fim_analise):
                        base = seq[i]
                        if base != ' ':  # Ignorar espaços (bases não alinhadas)
                            coverage[j] += 1
                        j = j + 1
                if not achou_comeco:
                    # Encontrar a posição inicial onde a cobertura é >= 3
                    start_position = 0
                    for i, cov in enumerate(coverage):
                        if cov >= coverage_threshold:
                            achou_comeco = True
                            start_position = i
                            break

                # Encontrar a posição final onde a cobertura é >= 3
                end_position = len(contig_seq)
                for i in range(len(coverage) - 1, -1, -1):
                    if coverage[i] >= coverage_threshold:
                        end_position = i + 1 + len(contig_final) - len(contig_seq)
                        break


                # Apagar sequências para vir próximo consensus.
                seqs = []

            line = file.readline()

    # Cortar o contig até as posições de início e fim
    trimmed_contig = contig_final[start_position:end_position]

    return trimmed_contig

# Exemplo de uso
cap3_output = 'consenso'  # Substitua pelo nome do seu arquivo de saída do CAP3
contig_name = 'Contig 1'  # Nome do contig que você deseja analisar
trimmed_contig = extract_contig_with_coverage(cap3_output, contig_name, coverage_threshold=3)
print("Contig cortado até a cobertura 3:", trimmed_contig)