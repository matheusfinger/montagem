def extract_contig_with_coverage(cap3_output, contig_name, coverage_threshold=3):
    # Variáveis para armazenar a sequência do contig e as sequências alinhadas
    contig_sequence = ''
    aligned_sequences = []

    # Ler o arquivo de saída do CAP3
    with open(cap3_output, 'r') as file:
        line = file.readline()
        cont = 0
        while cont < 2:
            if line.startswith(f'******************* {contig_name} ********************'):
                cont += 1
            line = file.readline()
        line = file.readline()
        
        consenso = False
        while not (line.startswith('******************* Contig')):
            # Coletar as sequências alinhadas de consenso em consenso

            """
            Tô precisando pensar aqui.
            """
            if not (line.startswith(' ')):
                aligned_sequences.append(line.split(" ")[1])  # Extrair a sequência
            line = file.readline()

    print(aligned_sequences)
    # Calcular a cobertura para cada posição
    coverage = [0] * len(contig_sequence)
    for seq in aligned_sequences:
        for i, base in enumerate(seq):
            if base != ' ':  # Ignorar espaços (bases não alinhadas)
                coverage[i] += 1
    print(coverage)
    # Encontrar a posição inicial onde a cobertura é >= 3
    start_position = 0
    for i, cov in enumerate(coverage):
        if cov >= coverage_threshold:
            start_position = i
            break

    # Encontrar a posição final onde a cobertura é >= 3
    end_position = len(contig_sequence)
    for i in range(len(coverage) - 1, -1, -1):
        if coverage[i] >= coverage_threshold:
            end_position = i + 1
            break

    # Cortar o contig até as posições de início e fim
    trimmed_contig = contig_sequence[start_position:end_position]

    return trimmed_contig

# Exemplo de uso
cap3_output = 'consenso'  # Substitua pelo nome do seu arquivo de saída do CAP3
contig_name = 'Contig 1'  # Nome do contig que você deseja analisar
trimmed_contig = extract_contig_with_coverage(cap3_output, contig_name, coverage_threshold=3)
print("Contig cortado até a cobertura 3:", trimmed_contig)