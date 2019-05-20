# Importa todos os módulos necessários para a execução do programa
from tkinter import *
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageTk

ds = pd.DataFrame()

janela = Tk() # Cria uma janela
janela.title('Curva de Fit') # Define o título da janela.

a1 = StringVar()
b1 = StringVar()
A1 = StringVar()
gama1 = StringVar()
M1 = StringVar()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def breitwigner(E, gamma, M, a, b, A):
    return a*E+b+A*((2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))))/((E**2-M**2)**2+M**2*gamma**2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Abrir():
    file =  filedialog.askopenfilename(initialdir = "/",title = "Abrir arquivo",filetypes = (("csv files","*.csv"), ("All files", "*.")))
#    arq = open(file, 'r') #Abre o arquivo desejado para ler os dados.
    global ds
    ds = pd.read_csv(file) #Lê o arquivo com os dados.
#    arq.close() #Fecha o arquivo.
    Imp_Graf()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Imp_Graf():
    inv_mass_log = np.log10(ds["M"])
    bins=300
    weights = []
    for i in ds["M"]:
        weights.append(bins/np.log(10)/i)

    plt.title('Histograma da massa invariante de dois múons \n')
    histogram = plt.hist(inv_mass_log, bins=bins, range=(-0.5,2.5), weights=weights, color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:] )
    plt.yscale('log')
    plt.xlabel('log10(Massa invariante) [log10(GeV)]')
    plt.ylabel('Número de eventos [log10]')
    plt.savefig('Hist.png') #Salva o gráfico como uma imagem PNG.
    plt.show()

    img = Image.open("Hist.png")
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=2, rowspan=5)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular():
    gamma = gama1.get()
    M = M1.get()
    a = a1.get()
    b = b1.get()
    A = A1.get()
    #print(gama.get(), M.get(), a.get(), b.get(), A.get())
    
    initials1 = [gamma, M, a, b, A]
    best1, covariance1 = curve_fit(breitwigner, x, y, p0=initials1, sigma=np.sqrt(y))
    error1 = np.sqrt(np.diag(covariance1))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Breitwigner():
    Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    Entry(janela, textvar=gama1).grid(row=6, column=2)
        
    Label(janela, text= "Entre com o valor da posição do máximo da distribuição (M):").grid(row=7, column=0, columnspan=2)
    Entry(janela, textvar=M1).grid(row=7, column=2)

    Label(janela, text= "Entre com o valor do parâmtro a que é a inclinação usada para perceber o efeito do background:").grid(row=8, column=0, columnspan=2)
    Entry(janela, textvar=a1).grid(row=8, column=2)

    Label(janela, text= "Entre com o valor do parâmetro b de intercepção no eixo y usada para perceber o efeito do background:").grid(row=9, column=0, columnspan=2)
    Entry(janela, textvar=b1).grid(row=9, column=2)

    Label(janela, text= "Entre com o valor do parâmetro A da altura da distribuição de Breit-Wigner:").grid(row=10, column=0, columnspan=2)
    Entry(janela, textvar=A1).grid(row=10, column=2)
    
    Button(janela, text="Calcular Fit", command=Calcular).grid(row=11, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Gauss():
    Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    sigma= Entry(janela).grid(row=6, column=2)
    
    Label(janela, text= "Entre com o valor da posição do máximo da distribuição:").grid(row=7, column=0, columnspan=2)
    media = Entry(janela).grid(row=7, column=2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
Button(janela, text='Abrir arquivo', font=('Times', '12', 'bold'), command=Abrir).grid(row=0, column=0)
#tk.Button(janela, text='Calcular fit', font=('Times', '12', 'bold')).grid(row=0, column=1)

Label(janela, text= "Escolha o pico para calcular o Fit?").grid(row=1, column=0, columnspan=2)

ck1 = Checkbutton(janela, text= "Rho", command=Gauss).grid(row=2, column=0)
ck2 = Checkbutton(janela, text= "Phi", command=Gauss).grid(row=2, column=1)
ck3 = Checkbutton(janela, text= "J/Psi", command=Gauss).grid(row=3, column=0)
ck4 = Checkbutton(janela, text= "Psi'", command=Gauss).grid(row=3, column=1)
ck5 = Checkbutton(janela, text= "Upsilon", command=Gauss).grid(row=4, column=0)
ck6 = Checkbutton(janela, text= "Z", command=Breitwigner).grid(row=4, column=1)
ck7 = Checkbutton(janela, text= "Todos (Estamos trabalhano nisto!)").grid(row=5, column=0, columnspan=2)
#----------------------------------------------------------------------------------------------------------------------------


janela.geometry("800x600")
janela.mainloop() #mantém a janela aberta