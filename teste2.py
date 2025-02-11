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

        achou_comeco_anterior = False
        comeco_anterior = 0
        fim_anterior = 0

        while line and not (line.startswith('******************* Contig')):
            
            if not(line.startswith("*") or line.startswith(" ") or len(line.strip()) == 0 or line.startswith("consensus")):
                seqs.append(line.strip("\n"))


            if line.startswith(contig_name) and not achou_comeco_anterior:
                linha = line.strip("\n")
                for i in range(len(line)-len(contig_seq), len(line)):
                    if seq[i] != ' ':
                        comeco_anterior = i + len(contig_final)
                        achou_comeco_anterior = True
                        break
            
            if line.startswith(contig_name):
                linha = line.strip("\n")
                for i in range(len(line)-1, len(line)-len(contig_seq)-1, -1):
                    if seq[i] != ' ':
                        fim_anterior = i + len(contig_final)
                        break




            if line.startswith("consensus"):
                line = line.strip()
                contig_seq = line.split(" ")
                contig_seq = contig_seq[len(contig_seq)-1]
                comeco_analise = len(line) - len(contig_seq)
                fim_analise = len(line)


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
                            start_position = i + len(contig_final)
                            break

                # Encontrar a posição final onde a cobertura é >= 3
                for i in range(len(coverage) - 1, -1, -1):
                    if coverage[i] >= coverage_threshold:
                        end_position = i + 1 + len(contig_final)
                        break

                contig_final += line[comeco_analise:fim_analise]

                # Apagar sequências para vir próximo consensus.
                seqs = []

            line = file.readline()

    start_position = min(comeco_anterior, start_position)
    end_position = max(fim_anterior, end_position)
    # Cortar o contig até as posições de início e fim
    trimmed_contig = contig_final[start_position:end_position]

    return trimmed_contig

print(extract_contig_with_coverage("consenso", "Contig 1"))