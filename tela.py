import os
import datetime
import PySimpleGUI as sg
import logging

diminuicoes = 1

# Configuração do logging
logging.basicConfig(filename='error_log.txt', level=logging.ERROR, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

class TelaPython:
    def __init__(self):
        sg.change_look_and_feel('LightGreen3')

        # Layout da interface gráfica
        layout = [
            [sg.Image('icone2.png', size=(500, 231))],
            [sg.Text('Contig', size=(13, 1)), sg.Multiline(default_text="Insira aqui o contig", size=(200, 7), key='seq')],
            [sg.Text('Localização da biblioteca', size=(20, 2)),
             sg.Input(size=(50, 0), key='diretorio'), sg.FileBrowse()],
            [sg.Frame('', [[
                sg.Checkbox('A biblioteca é paired end?', key='paired'),
                sg.Text('', size=(15, 2)),
                sg.Text('Localização da segunda biblioteca', size=(25, 2)),
                sg.Input(size=(50, 0), key='diretorio2'), sg.FileBrowse()]])],
            [
                sg.Text('Qual será o tamanho do k-mer utilizado?', size=(40, 2)),
                sg.Input(size=(15,0), key='kmer_size')
            ],
            [sg.Text("Deseja alongar a extremidade 3' ou 5'?"),
             sg.Radio("5'", "extremidade", key="5", default=True), sg.Radio("3'", "extremidade", key="3")],
            [sg.Frame('', [[
                sg.Checkbox('Deseja que o k-mer seja antes das extremidades?', key='voltar'),
                sg.Text('', size=(15, 2)),
                sg.Text('Se você marcou sim, quantos nucleotídeos deseja se afastar da extremidade?', size=(45, 3)),
                sg.Input(size=(15, 0), key='voltar_nuc')]])],
            [sg.Checkbox('Deseja usar o modo TURBO?', key='turbo')],
            [sg.Button('Enviar dados')],
            [sg.Output(size=(200, 15), key='saida')]
        ]

        self.janela = sg.Window("Montagem", layout)

    def iniciar(self):
        while True:
            event, values = self.janela.read()

            if event in (sg.WINDOW_CLOSED, 'Exit'):
                break

            try:
                self.processar_sequenciamento(values)
            except Exception as e:
                logging.error(f"Erro ao processar sequenciamento: {e}", exc_info=True)
                sg.popup_error(f"Ocorreu um erro: {e}\nVerifique o arquivo error_log.txt para mais detalhes.")

    def processar_sequenciamento(self, values):
        seq = values['seq']
        diretorio = values['diretorio']
        try:
            kmer_size = int(values['kmer_size'])
        except ValueError:
            sg.popup_error("Insira um número válido para o tamanho do k-mer.")
            return

        if not diretorio or not os.path.exists(diretorio):
            sg.popup_error("Arquivo da biblioteca não encontrado!")
            return

        paired = values['paired']
        diretorio2 = values['diretorio2']
        extremidade5 = values['5']
        extremidade3 = values['3']
        voltar = values['voltar']
        voltar_nuc = 0

        if values['voltar_nuc']:
            try:
                voltar_nuc = int(values['voltar_nuc'])
            except ValueError:
                sg.popup_error("Insira um número válido em 'voltar nucleotídeos'.")
                return

        turbo = values['turbo']
        self.janela['saida'].update('')  # Limpa o output

        contig = self.preparar_contig(seq)

        if not turbo and extremidade5:
            kmer = self.definir_kmer(contig, kmer_size, voltar, voltar_nuc, inicio=True)
        elif not turbo:
            kmer = self.definir_kmer(contig, kmer_size, voltar, voltar_nuc, inicio=False)
        else:
            kmer1 = self.definir_kmer(contig, kmer_size, voltar, voltar_nuc, inicio=True)
            kmer2 = self.definir_kmer(contig, kmer_size, voltar, voltar_nuc, inicio=False)
            kmer = (kmer1, kmer2)

        self.remover_arquivos_passados()

        if turbo:
            contig_final, tamanho = self.modo_turbo(diretorio, diretorio2, paired, kmer, kmer_size, contig, voltar, voltar_nuc)
        else:
            contig_final, tamanho = self.modo_normal(diretorio, diretorio2, paired, kmer, kmer_size, contig)

        print(f"O maior contig gerado foi:\n{contig_final}\nO contig possui tamanho {tamanho}.")

    def preparar_contig(self, seq):
        if seq[0] == '>':
            comeco = 1
        else:
            comeco = 0

        contig = ""
        quebras = seq.split("\n")

        for i in range(comeco, len(quebras)):
            contig += quebras[i]

        return contig

    def definir_kmer(self, contig, kmer_size, voltar, voltar_nuc, inicio=True):
        if inicio:
            if voltar:
                return contig[0 + voltar_nuc:kmer_size + voltar_nuc]
            else:
                return contig[0:kmer_size]
        else:
            if voltar:
                return contig[len(contig) - kmer_size - voltar_nuc:len(contig) - voltar_nuc]
            else:
                return contig[len(contig) - kmer_size:len(contig)]

    def remover_arquivos_passados(self):
        print("Apagando resultados anteriores...")
        if os.path.exists("matches.fa.cap.contigs"):
            os.remove("matches.fa.cap.contigs")
        if os.path.exists("resultado.fa"):
            os.remove("resultado.fa")
        if os.path.exists("resultado1.fa"):
            os.remove("resultado1.fa")
        if os.path.exists("resultado2.fa"):
            os.remove("resultado2.fa")

    def reverso_complemento(self, seq):
        res = ""
        for i in range(len(seq)-1, -1, -1):
            if seq[i].upper() == 'A':
                res += 'T'
            elif seq[i].upper() == 'T':
                res += 'A'
            elif seq[i].upper() == 'C':
                res += 'G'
            elif seq[i].upper() == 'G':
                res += 'C'
            else:
                res += seq[i]
        return res


    def modo_turbo(self, diretorio, diretorio2, paired, kmer, kmer_size, contig, voltar, voltar_nuc):
        continuar = True
        volta = 0
        pula = False
        num_pula = 0
        id_contig = None
        while continuar:
            print(f'{volta + 1}ª iteração iniciada.')
            print(datetime.datetime.now())
            
            if not pula:
                match = self.executar_bbduk(diretorio, diretorio2, paired, kmer, kmer_size)

            print(f'Análise das reads finalizada. {match} reads mapearam.')
            contigs = self.realizar_montagem(contig, id_contig = id_contig, turbo=True, iteracao = volta + 1)
            contig_fim = contigs[1]
            id_contig = contigs[2]
            contig = contig_fim
            tam_contig = len(contig_fim)
            print(f'Tamanho do contig gerado nessa iteração: {tam_contig}')

            if match >= 3 and volta < 50 and tam_contig < 3000 and num_pula < 4:
                volta += 1
                antigo = kmer
                kmer1 = self.definir_kmer(contig_fim, kmer_size, voltar, voltar_nuc, inicio=True)
                kmer2 = self.definir_kmer(contig_fim, kmer_size, voltar, voltar_nuc, inicio=False)
                kmer = (kmer1, kmer2)
                kmer0_rev = self.reverso_complemento(kmer[0])
                kmer1_rev = self.reverso_complemento(kmer[1])
                if antigo[0] == kmer[0]:
                    if antigo[1] == kmer[1] or antigo[1] == kmer1_rev:
                        if pula:
                            num_pula += 1
                        else:
                            num_pula = 1
                        pula = True
                    else:
                        pula = False
                elif antigo[0] == kmer0_rev:
                    if antigo[1] == kmer[1] or antigo[1] == kmer1_rev:
                        if pula:
                            num_pula += 1
                        else:
                            num_pula = 1
                        pula = True
                    else:
                        pula = False
                else:
                    pula = False
            else:
                contig_fim = contigs[0]
                continuar = False
                if match < 3:
                    print(f'Parou na {volta + 1}ª iteração por falta de matches.')
                elif tam_contig >= 3000:
                    print(f'Parou pq o tamanho do contig passou do corte')
                elif num_pula >= 4:
                    print(f'Ta repetindo os kmers')
                else:
                    print('Número máximo de iterações alcançado.')

        return contig_fim, len(contig_fim)

    def modo_normal(self, diretorio, diretorio2, paired, kmer, kmer_size, contig):
        print("Executando bbduk")
        match = self.executar_bbduk(diretorio, diretorio2, paired, kmer, kmer_size)
        print("Executando cap3")
        contig_fim = self.realizar_montagem(contig, id_contig=None, turbo=False, iteracao=1)

        return contig_fim, len(contig_fim)

    def executar_bbduk(self, diretorio, diretorio2, paired, kmer, kmer_size):
        if type(kmer) is tuple:
            kmer1, kmer2 = kmer
            if paired:
                print("Analisando kmer 5'")
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado1.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer1}')
                print("Analisando kmer 3'")
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado2.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer2}')
            else:
                print("Analisando kmer 5'")
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer1}')
                print("Analisando kmer 3'")
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer2}')
            return self.contar_reads("resultado1.fasta") + self.contar_reads("resultado2.fasta")
        else:
            if paired:
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer}')
            else:
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp=True mm=f k={kmer_size} threads=6 literal={kmer}')

            return self.contar_reads("resultado.fasta")

    def contar_reads(self, arquivo):
        match = 0

        with open(arquivo, 'r') as file:
            for linha in file:
                if linha.startswith('>'):
                    match += 1

        return match

    def tamanho_contigs(self, arquivo: str):
        tamanhos = []

    def contig_with_most_reads(self, arquivo, contig_anterior):
        max_reads = 0
        num_reads = 0
        max_contig = None
        # Checa se é o último contig
        last_contig = False
        with open(arquivo) as arq:
            line = arq.readline()
            # Lê até achar o começo da descrição dos contigs
            while line and not line.startswith('******************* Contig'):
                line = arq.readline()
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
                    num_reads = 0
                    actual_contig = line.strip()
                    actual_contig = actual_contig.strip("*")
                    actual_contig = actual_contig.strip(" ")
                    # Reseta variável dizendo que contém o contig anterior
                    contains = False
                    if last_contig:
                        break
                elif not line.startswith('*'):
                    if (contig_anterior is not None and line.startswith(contig_anterior)) or (contig_anterior is None):
                        contains = True
                    num_reads += 1
                line = arq.readline()
                if line.startswith('DETAILED DISPLAY OF CONTIGS'):
                    last_contig = True
                    num_reads -= 1
        if max_contig is not None:
            return max_contig
        else:
            print("Nessa iteração não foi possível achar contig com o contig anterior")
            return self.contig_with_most_reads("consenso", contig_anterior=None)
    
    def extract_contig_with_coverage(self, cap3_output, contig_name, contig_anterior, iteracao, coverage_threshold=3):
        # Ler o arquivo de saída do CAP3
        with open(cap3_output, 'r') as file:
            line = file.readline()
            cont = 0
            while line and cont < 2:
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
            if contig_anterior is not None:
                print(f'O contig anterior era: {contig_anterior}')

            while line and not (line.startswith('******************* Contig')):
                
                if not(line.startswith("*") or line.startswith(" ") or len(line.strip()) == 0 or line.startswith("consensus")):
                    seqs.append(line.strip("\n"))


                if contig_anterior is not None:
                    if line.startswith(contig_anterior) and not achou_comeco_anterior:
                        line = line.strip("\n")
                        seq = line
                        comeco_seq = line.split(" ")
                        comeco_seq = " ".join(comeco_seq[1:])
                        comeco_seq = len(comeco_seq)
                        for i in range(len(line)-comeco_seq, len(line)):
                            if seq[i] != ' ':
                                comeco_anterior = i + len(contig_final)
                                achou_comeco_anterior = True
                                break
                    
                    if line.startswith(contig_anterior):
                        line = line.strip("\n")
                        seq = line
                        comeco_seq = line.split(" ")
                        comeco_seq = " ".join(comeco_seq[1:])
                        comeco_seq = len(comeco_seq)
                        for i in range(len(line)-1, len(line)-comeco_seq-1, -1):
                            if seq[i] != ' ':
                                fim_anterior = i + len(contig_final)
                                break




                if line.startswith("consensus"):
                    line = line.strip()
                    contig_seq = line.split(" ")
                    contig_seq = contig_seq[len(contig_seq)-1]
                    comeco_analise = len(line) - len(contig_seq)
                    fim_analise = len(line)
                    end_position = None

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

        if achou_comeco_anterior:
            if end_position is None:
                end_position = fim_anterior
            if contig_anterior is not None:
                print(f'Contig anterior: {contig_anterior}')
            print(f'Começo anterior: {comeco_anterior}')
            print(f'Fim anterior: {comeco_anterior}')
            print(f'Start_position: {start_position}')
            if end_position is not None:
                print(f'End position: {end_position}')
            start_position = min(comeco_anterior, start_position)
            end_position = max(fim_anterior, end_position)
        elif iteracao != 1:
            print("O contig da iteração anterior não está presente nos contigs da iteração atual.")
        # Cortar o contig até as posições de início e fim
        trimmed_contig = contig_final[start_position:end_position]

        return trimmed_contig

    def realizar_montagem(self, contig, id_contig, turbo, iteracao):
        if turbo:
            arq1 = open("resultado1.fasta")
            filtradas1 = arq1.readlines()
            arq1.close()
            arq2 = open("resultado2.fasta")
            filtradas2 = arq2.readlines()
            arq2.close()

        else:
            arq = open("resultado.fasta")
            filtradas = arq.readlines()
            arq.close()
        
        if (not turbo) or (turbo and iteracao == 1):
            with open('matches.fa', 'w') as arquivo:
                arquivo.write(">contig1\n")
                arquivo.write(f"{contig}\n")
                arquivo.write(">contig2\n")
                arquivo.write(f"{contig}\n")
                if turbo:
                    for linha in filtradas1:
                        arquivo.write(linha)
                    for linha in filtradas2:
                        arquivo.write(linha)
                if not turbo:
                    for linha in filtradas:
                        arquivo.write(linha)
        else:
            with open('matches.fa', 'w') as arquivo_matches:
                with open('matches.fa.cap.contigs', 'r') as arquivo_contigs:
                    for linha in arquivo_contigs:
                        arquivo_matches.write(linha)
                for linha in filtradas1:
                    arquivo_matches.write(linha)
                for linha in filtradas2:
                    arquivo_matches.write(linha)

        print("Iniciando montagem.")
        os.system("cap3 matches.fa -p 98 > consenso")

        id_contig_formado = self.contig_with_most_reads("consenso", id_contig)


        if turbo:
            contig_trimmado = self.extract_contig_with_coverage("consenso", id_contig_formado, id_contig, iteracao)

            if len(contig_trimmado) < len(contig):
                print("DIMINUIU DIMINUIU DIMINUIU")
                raise ValueError
        
        id_contig_formado = id_contig_formado.replace(" ", "")
        contig_formado = ""

        with open('matches.fa.cap.contigs', 'r') as arquivo:
            linha = arquivo.readline()
            while linha and not linha.startswith(f'>{id_contig_formado}'):
                linha = arquivo.readline()
            linha = arquivo.readline()
            while linha and not linha.startswith('>'):
                contig_formado += linha.strip()
                linha = arquivo.readline()

        print("Montagem feita.")

        if turbo:
            return (contig_formado, contig_trimmado, id_contig_formado)
        else: return contig_formado

if __name__ == "__main__":
    tela = TelaPython()
    tela.iniciar()