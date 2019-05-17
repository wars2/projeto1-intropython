# Importa todos os módulos necessários para a execução do programa
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageTk
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Abrir():
    file =  filedialog.askopenfilename(initialdir = "/",title = "Abrir arquivo",filetypes = (("csv files","*.csv"), ("All files", "*.")))
#    arq = open(file, 'r') #Abre o arquivo desejado para ler os dados.
    dados = pd.read_csv(file) #Lê o arquivo com os dados.
#    arq.close() #Fecha o arquivo.
    return dados
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def breitwigner():
    # Valores iniciais para a otimização da curva.
    # E  = Energia
    # gamma Largura a meia altura (FWHM) da distribuição
    # M posição onde se encontra o valor máximo da distribuição
    # a inclinação usada para perceber o efeito do background
    # b valor de intercepção no eixo y usada para perceber o efeito do background
    # A a altura da distribuição de Breit-Wigner
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
    tk.Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    tk.Entry(janela).grid(row=6, column=2)

    tk.Label(janela, text= "Entre com o valor da posição do máximo da distribuição (M):").grid(row=7, column=0, columnspan=2)
    tk.Entry(janela).grid(row=7, column=2)
    
    tk.Label(janela, text= "Entre com o valor do parâmtro a que é a inclinação usada para perceber o efeito do background:").grid(row=8, column=0, columnspan=2)
    tk.Entry(janela).grid(row=8, column=2)

    tk.Label(janela, text= "Entre com o valor do parâmetro b de intercepção no eixo y usada para perceber o efeito do background:").grid(row=9, column=0, columnspan=2)
    tk.Entry(janela).grid(row=9, column=2)

    tk.Label(janela, text= "Entre com o valor do parâmetro A da altura da distribuição de Breit-Wigner:").grid(row=10, column=0, columnspan=2)
    tk.Entry(janela).grid(row=10, column=2)

#    return a*E+b+A*((2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))))/((E**2-M**2)**2+M**2*gamma**2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def gauss():
    tk.Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    tk.Entry(janela).grid(row=6, column=2)

    tk.Label(janela, text= "Entre com o valor da posição do máximo da distribuição:").grid(row=7, column=0, columnspan=2)
    tk.Entry(janela).grid(row=7, column=2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
janela = tk.Tk() # Cria uma janela
janela.title('Curva de fit') # Define o título da janela.

tk.Button(janela, text='Abrir arquivo', font=('Times', '12', 'bold'), command=Abrir).grid(row=0, column=0)
#tk.Button(janela, text='Calcular fit', font=('Times', '12', 'bold')).grid(row=0, column=1)

tk.Label(janela, text= "Escolha o pico para calcular o Fit?").grid(row=1, column=0, columnspan=2)

ck1 = tk.Checkbutton(janela, text= "Rho", command=gauss).grid(row=2, column=0)
ck2 = tk.Checkbutton(janela, text= "Phi", command=gauss).grid(row=2, column=1)
ck3 = tk.Checkbutton(janela, text= "J/Psi", command=gauss).grid(row=3, column=0)
ck4 = tk.Checkbutton(janela, text= "Psi'", command=gauss).grid(row=3, column=1)
ck5 = tk.Checkbutton(janela, text= "Upsilon", command=gauss).grid(row=4, column=0)
ck6 = tk.Checkbutton(janela, text= "Z", command=breitwigner).grid(row=4, column=1)
ck7 = tk.Checkbutton(janela, text= "Todos").grid(row=5, column=0)
#----------------------------------------------------------------------------------------------------------------------------
ds = Abrir()
inv_mass_log = np.log10(ds["M"])
bins=300
weights = []
for a in ds["M"]:
    weights.append(bins/np.log(10)/a)

plt.title('Histograma da massa invariante de dois múons \n')
plt.hist(inv_mass_log, bins=bins, range=(-0.5,2.5), weights=weights, color="cyan")
plt.yscale('log')
plt.xlabel('log10(Massa invariante) [log10(GeV)]')
plt.ylabel('Número de eventos [log10]')
plt.savefig('Hist.png') #Salva o gráfico como uma imagem PNG.
plt.show()

img = Image.open("Hist.png")
tkimage = ImageTk.PhotoImage(img)
tk.Label(janela, image=tkimage).grid(row=1, column=2, rowspan=5)

janela.geometry("800x600")
janela.mainloop() #mantém a janela aberta