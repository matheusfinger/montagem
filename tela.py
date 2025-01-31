import os
import time
import datetime
import PySimpleGUI as sg
import logging

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
            [sg.Checkbox('Considerar reverso complemento', key='rev', default=True)],
            [sg.Text("Deseja alongar a extremidade 3' ou 5'?"),
             sg.Radio("5'", "extremidade", key="5", default=True), sg.Radio("3'", "extremidade", key="3")],
            [sg.Frame('', [[
                sg.Checkbox('Deseja que o 25-mer seja antes das extremidades?', key='voltar'),
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

        if not diretorio or not os.path.exists(diretorio):
            sg.popup_error("Arquivo da biblioteca não encontrado!")
            return

        paired = values['paired']
        diretorio2 = values['diretorio2']
        rev = values['rev']
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
            kmer = self.definir_kmer(contig, voltar, voltar_nuc, inicio=True)
        elif not turbo:
            kmer = self.definir_kmer(contig, voltar, voltar_nuc, inicio=False)
        else:
            kmer1 = self.definir_kmer(contig, voltar, voltar_nuc, inicio=True)
            kmer2 = self.definir_kmer(contig, voltar, voltar_nuc, inicio=False)
            kmer = (kmer1, kmer2)

        self.remover_arquivos_passados()

        if turbo:
            contig_final, tamanho = self.modo_turbo(diretorio, diretorio2, paired, rev, kmer, contig)
        else:
            contig_final, tamanho = self.modo_normal(diretorio, diretorio2, paired, rev, kmer, contig)

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

    def definir_kmer(self, contig, voltar, voltar_nuc, inicio=True):
        if inicio:
            if voltar:
                return contig[0 + voltar_nuc:25 + voltar_nuc]
            else:
                return contig[0:25]
        else:
            if voltar:
                return contig[len(contig) - 25 - voltar_nuc:len(contig) - voltar_nuc]
            else:
                return contig[len(contig) - 25:len(contig)]

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


    def modo_turbo(self, diretorio, diretorio2, paired, rev, kmer, contig):
        continuar = True
        volta = 0
        pula = False
        num_pula = 0
        while continuar:
            print(f'{volta + 1}ª iteração iniciada.')
            print(datetime.datetime.now())
            
            if not pula:
                match = self.executar_bbduk(diretorio, diretorio2, paired, rev, kmer)

            print(f'Análise das reads finalizada. {match} reads mapearam.')
            contig_fim = self.realizar_montagem(contig, turbo=True)
            contig = contig_fim
            tam_contig = len(contig_fim)
            print(f'Tamanho do contig gerado nessa iteração: {tam_contig}')

            if match >= 3 and volta < 50 and tam_contig < 3000 and num_pula < 4:
                volta += 1
                antigo = kmer
                kmer1 = self.definir_kmer(contig_fim, None, False, inicio=True)
                kmer2 = self.definir_kmer(contig_fim, None, False, inicio=False)
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

    def modo_normal(self, diretorio, diretorio2, paired, rev, kmer, contig):
        print("Executando bbduk")
        match = self.executar_bbduk(diretorio, diretorio2, paired, rev, kmer)
        print("Executando cap3")
        contig_fim = self.realizar_montagem(contig, turbo=False)

        return contig_fim, len(contig_fim)

    def executar_bbduk(self, diretorio, diretorio2, paired, rev, kmer):
        if type(kmer) is tuple:
            kmer1, kmer2 = kmer
            if paired:
                print("Analisando kmer 5'")
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado1.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer1}')
                print("Analisando kmer 3'")
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado2.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer2}')
            else:
                print("Analisando kmer 5'")
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer1}')
                print("Analisando kmer 3'")
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer2}')
            return self.contar_reads("resultado1.fasta") + self.contar_reads("resultado2.fasta")
        else:
            if paired:
                os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer}')
            else:
                os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp={rev} mm=f k=25 threads=6 literal={kmer}')

            return self.contar_reads("resultado.fasta")

    def contar_reads(self, arquivo):
        match = 0

        with open(arquivo, 'r') as file:
            for linha in file:
                if linha.startswith('>'):
                    match += 1

        return match

    def realizar_montagem(self, contig, turbo):
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

        print("Iniciando montagem.")
        os.system("cap3 matches.fa -p 99 > consenso")

        contigs = []
        contig_atual = ""

        with open('matches.fa.cap.contigs', 'r') as arquivo:
            for linha in arquivo:
                if linha.startswith('>'):
                    if contig_atual:
                        contigs.append(contig_atual)
                    contig_atual = ""
                else:
                    contig_atual += linha.strip()

        if contig_atual:
            contigs.append(contig_atual)

        maior_contig = max(contigs, key=len, default="")
        print("Montagem feita.")

        return maior_contig

if __name__ == "__main__":
    tela = TelaPython()
    tela.iniciar()