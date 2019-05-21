# Importa todos os módulos necessários para a execução do programa
from tkinter import *
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageTk
from scipy.optimize import curve_fit

ds = pd.DataFrame()
x = ()
y = ()
inv_mass_log = ()

janela = Tk() # Cria uma janela
janela.title('Curva de Fit') # Define o título da janela.

a1 = StringVar()
b1 = StringVar()
A1 = StringVar()
gama1 = StringVar()
M1 = StringVar()
c1 = IntVar()
c2 = IntVar()
c3 = IntVar()
c4 = IntVar()
c5 = IntVar()
c6 = IntVar()

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def breitwigner(E, gamma, M, a, b, A):
    return a*E+b+A*((2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))))/((E**2-M**2)**2+M**2*gamma**2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Verifica():
    if ((c1.get() == 1 or c2.get() == 1 or c3.get() == 1 or c4.get() == 1 or c5.get() == 1) and c6.get() == 0):
        Gauss()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 1):
        Breitwigner()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        canvas = Canvas(janela).grid(row=6, column=0, rowspan=7, columnspan=3, stick=N+S+E+W)
    else:
        print("Não se pode calcular mais de 1 fit no momento!")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Abrir():
    global ds
    
    file =  filedialog.askopenfilename(initialdir = "/",title = "Abrir arquivo",filetypes = (("csv files","*.csv"),("CSV files", "*.csv")))
    ds = pd.read_csv(file) #Lê o arquivo com os dados.

    Imp_Graf()
        
    Label(janela, text= "Escolha o pico para calcular o Fit?").grid(row=1, column=0, columnspan=2)
    
    Checkbutton(janela, text= "Rho", variable = c1, command=Verifica).grid(row=2, column=0)
    Checkbutton(janela, text= "Phi", variable = c2, command=Verifica).grid(row=2, column=1)
    Checkbutton(janela, text= "J/Psi", variable = c3, command=Verifica).grid(row=3, column=0)
    Checkbutton(janela, text= "Psi'", variable = c4, command=Verifica).grid(row=3, column=1)
    Checkbutton(janela, text= "Upsilon", variable = c5, command=Verifica).grid(row=4, column=0)
    Checkbutton(janela, text= "Z", variable = c6, command=Verifica).grid(row=4, column=1)
    Checkbutton(janela, text= "Todos (Estamos trabalhano nisto!)").grid(row=5, column=0, columnspan=2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Imp_Graf():
    global x
    global y
    global inv_mass_log
    
    inv_mass_log = np.log10(ds["M"])
    bins=300
    weights = []
    for i in ds["M"]:
        weights.append(bins/np.log(10)/i)

    plt.title('Histograma da massa invariante de dois múons \n')
    histogram = plt.hist(inv_mass_log, bins=bins, range=(-0.5,2.5), weights=weights, color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    plt.yscale('log')
    plt.xlabel('log10(Massa invariante) [log10(GeV)]')
    plt.ylabel('Número de eventos [log10]')
    plt.savefig('Hist.png') #Salva o gráfico como uma imagem PNG.
    plt.show()

    img = Image.open("Hist.png")
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=2, rowspan=5)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular1():
#    global gamma
#    global M
#    global a
#    global b
#    global A
    
    gamma = float(gama1.get())
    M = float(M1.get())
    a = float(a1.get())
    b = float(b1.get())
    A = float(A1.get())
    
    limitedmasses = inv_mass_log[(inv_mass_log > 1.85) & (inv_mass_log < 2.05)]
    histogram = plt.hist(limitedmasses, bins=bins, range=(1.85,2.05), color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    
    initials = [gamma, M, a, b, A]
    best, covariance = curve_fit(breitwigner, x, y, p0=initials, sigma=np.sqrt(y))
    error = np.sqrt(np.diag(covariance))
    plt.plot(x, breitwigner(x, *best), 'r-', label='gamma = {}, M = {}'.format(best[0], best[1]))
    plt.yscale("log")
    plt.xlabel('Massa Invariante (log10) [GeV]')
    plt.ylabel('Número de eventos (log10)')
    plt.title('Fit de Breit-Wigner para o pico Z')
    plt.savefig("Hist.png")
    plt.show()
    
    img = Image.open("Hist.png")
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=2, rowspan=5)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Breitwigner():
    canvas = Canvas(janela).grid(row=6, column=0, rowspan=7, columnspan=3, stick=N+S+E+W)
    
    Label(canvas, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    Entry(canvas, textvar=gama1).grid(row=6, column=2)
    
    Label(canvas, text= "Entre com o valor da posição do máximo da distribuição (M):").grid(row=7, column=0, columnspan=2)
    Entry(canvas, textvar=M1).grid(row=7, column=2)

    Label(canvas, text= "Entre com o valor do parâmtro a que é a inclinação usada para perceber o efeito do background:").grid(row=8, column=0, columnspan=2)
    Entry(canvas, textvar=a1).grid(row=8, column=2)

    Label(canvas, text= "Entre com o valor do parâmetro b de intercepção no eixo y usada para perceber o efeito do background:").grid(row=9, column=0, columnspan=2)
    Entry(canvas, textvar=b1).grid(row=9, column=2)

    Label(canvas, text= "Entre com o valor do parâmetro A da altura da distribuição de Breit-Wigner:").grid(row=10, column=0, columnspan=2)
    Entry(canvas, textvar=A1).grid(row=10, column=2)
    
    Button(canvas, text="Calcular Fit", command=Calcular1).grid(row=11, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Gauss():
    canvas = Canvas(janela).grid(row=6, column=0, rowspan=7, columnspan=3, stick=N+S+E+W)
    
    Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=6, column=0, columnspan=2)
    sigma= Entry(janela).grid(row=6, column=2)
    
    Label(janela, text= "Entre com o valor da posição do máximo da distribuição:").grid(row=7, column=0, columnspan=2)
    media = Entry(janela).grid(row=7, column=2)
        
    Button(canvas, text="Calcular Fit", command=Calcular).grid(row=8, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
Button(janela, text='Abrir arquivo', font=('Times', '12', 'bold'), width = 15, command=Abrir).grid(row=0, column=0)

#janela.geometry("800x600")
janela.mainloop() #mantém a janela aberta