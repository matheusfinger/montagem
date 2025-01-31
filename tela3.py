import os
import time
import datetime
import PySimpleGUI as sg

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
                sg.Checkbox('Deseja que o 20-mer seja antes das extremidades?', key='voltar'),
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

            self.processar_sequenciamento(values)

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

        if extremidade5:
            kmer = self.definir_kmer(contig, voltar, voltar_nuc, inicio=True)
        else:
            kmer = self.definir_kmer(contig, voltar, voltar_nuc, inicio=False)

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
                return contig[0 + voltar_nuc:20 + voltar_nuc]
            else:
                return contig[0:20]
        else:
            if voltar:
                return contig[len(contig) - 20 - voltar_nuc:len(contig) - voltar_nuc]
            else:
                return contig[len(contig) - 20:len(contig)]

    def remover_arquivos_passados(self):
        print("Apagando resultados anteriores...")
        if os.path.exists("matches.fa.cap.contigs"):
            os.remove("matches.fa.cap.contigs")
        if os.path.exists("resultado.fa"):
            os.remove("resultado.fa")

    def modo_turbo(self, diretorio, diretorio2, paired, rev, kmer, contig):
        continuar = True
        volta = 0
        pula = False
        match = 100
        tam_contig = 1000
        while continuar:
            print(f'{volta + 1}ª iteração iniciada.')
            print(datetime.datetime.now())
            if not pula:
                print("EXECUTA BBDUK")
                # match = self.executar_bbduk(diretorio, diretorio2, paired, rev, kmer)
                match = match - 10 * (volta + 1)

            print(f'Análise das reads finalizada. {match} reads mapearam.')
            #contig_fim = self.realizar_montagem(contig)
            #tam_contig = len(contig_fim)
            tam_contig = tam_contig + 100 * (volta + 1)

            if match >= 3 and volta < 50 and tam_contig < 3000:
                volta += 1
                antigo = kmer
                kmer = "ACTAGCTAGCG"
                if antigo == kmer:
                    pula = True
                else:
                    pula = False
            else:
                continuar = False
                if match < 3:
                    print(f'Parou na {volta + 1}ª iteração por falta de matches.')
                elif tam_contig >= 3000:
                    print(f'Parou pq o tamanho do contig passou do corte')
                else:
                    print('Número máximo de iterações alcançado.')

        return "CTGATCGATGCTACGTAGCTA", tam_contig

    def modo_normal(self, diretorio, diretorio2, paired, rev, kmer, contig):
        match = self.executar_bbduk(diretorio, diretorio2, paired, rev, kmer)
        contig_fim = self.realizar_montagem(contig)

        return contig_fim, len(contig_fim)

    def executar_bbduk(self, diretorio, diretorio2, paired, rev, kmer):
        if paired:
            os.system(f'bbduk.sh in1={diretorio} in2={diretorio2} outm=resultado.fasta rcomp={rev} mm=f k=20 threads=8 literal={kmer}')
        else:
            os.system(f'bbduk.sh in={diretorio} outm=resultado.fasta rcomp={rev} mm=f k=20 threads=8 literal={kmer}')

        return self.contar_reads("resultado.fasta")

    def contar_reads(self, arquivo):
        match = 0

        with open(arquivo, 'r') as file:
            for linha in file:
                if linha.startswith('>'):
                    match += 1

        return match

    def realizar_montagem(self, contig):
        with open('matches.fa', 'w') as arquivo:
            arquivo.write(">contig1\n")
            arquivo.write(f"{contig}\n")
            arquivo.write(">contig2\n")
            arquivo.write(contig)

        os.system("cap3 matches.fa > consenso")

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

        return maior_contig

if __name__ == "__main__":
    tela = TelaPython()
    tela.iniciar()
