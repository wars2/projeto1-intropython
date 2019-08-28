"""
Universidade do Estado do Rio de Janeiro.
Instituto de Física Armando Dias Tavares.
Professoras: Clemencia M. Herrera e Helena B. Malbouisson
Autores: Apóllo V. da Silva, Gabriel M. S. Campos, Hamilton S. Gama Filho, Hugo B. S. Borges,
         Luiz Carlos R. Silva Junior e Miguel A. G. Pereira.
Este programa foi desenvolvido durante as aulas de Introdução à Python lecionada na UERJ.
"""

# Importa todos os módulos necessários para a execução do programa
import platform
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageTk
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from funcoes_encap import funcoes as f

ds = pd.DataFrame()
x = ()
y = ()
inv_mass_log = ()
system = platform.system() # Diz qual SO está sendo utilizado.

# Criação da janela com o tkinter.
janela = Tk() # Cria uma janela
janela.title('Curva de Fit') # Define o título da janela.

if (system == "Linux"):
    janela.wm_iconbitmap('@iconpy.xbm') # Insere o ícone nas janelas do programa no Sistema Operacional Linux.
elif (system == 'Windows'):
    janela.wm_iconbitmap('iconpy.ico') # Insere o ícone nas janelas do programa no Sistema Operacional Windows.
else:
    print("Não implementamos este programa para MAC ou outros sistemas operacionais.")

e1 = StringVar()
e2 = StringVar()
e3 = StringVar()
e4 = StringVar()
e5 = StringVar()
e6 = StringVar()
e7 = StringVar()
e8 = StringVar()

c1 = IntVar()
c2 = IntVar()
c3 = IntVar()
c4 = IntVar()
c5 = IntVar()
c6 = IntVar()

global esperado
esperado = 0
global err_esperado
err_esperado = 0
global esperado1 
esperado1 = 0
global err_esperado1
err_esperado1 = 0

low_lim = -0.5
up_lim = 2.5
bins = 300
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Carregar_img(image): # Método para chamar a imagem gerada no primeiro Imp_Graf().
    global img
    global tkimage
    
    img = Image.open(image)
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=6, rowspan=20)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def convergencia(f, funcao, best, error, n):
    if (f == 0 and f == 1):
        fit = ("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! \n")
    elif (f >= 1 and f <= 19):
        fit = ("O fit convergiu! \n")
    else:
        fit = ("O fit provavelmente esta divergindo... Tente outros valores iniciais! \n")
    
    valor = varfit.get()
    if (valor == options[1]):     # Relativo a escolher a distribuição Duas Gaussianas + Exponencial
        ch2, pval = chisquare(y, funcao(x, *best))
        test = ch2 /(bins - 8)
        discrepancia = np.absolute(best[1] - esperado)
        discrepancia1 = np.absolute(best[4] - esperado1)
        erro_padrao = np.sqrt(error[1]**2 + err_esperado**2)
        erro_padrao1 = np.sqrt(error[4]**2 + err_esperado1**2)

        first = ('A discrepancia do valor obtido no primeiro pico com o valor esperado e: %6.4f \n' %(discrepancia))
        second = ('A discrepancia do valor obtido no segundo pico com o valor esperado e: %6.4f \n' %(discrepancia1))
        third = ('A raiz da soma quadratica do erro do valor obtido no primeiro pico com o erro do valor esperado e: %6.4f \n' %(erro_padrao))
        fourth = ('A raiz da soma quadratica do erro do valor obtido no segundo pico com o erro do valor esperado e: %6.4f \n' %(erro_padrao1))
        fifth = ('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f \n' % (test))
        sixth = ("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pelas discrepancias: \n")
        if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao and discrepancia < 2*erro_padrao1):
            seventh = ("O fit e bom!")
        else:
            seventh = ("O fit e ruim! Tente outros valores.")
        message2 = str(first +"\n"+ second +"\n"+ third +"\n"+ fourth +"\n"+ fifth +"\n"+ sixth +"\n"+ seventh)

    else:
        ch2, pval = chisquare(y, funcao(x, *best))
        test = ch2 /(bins - n)
        discrepancia = np.absolute(best[1] - esperado)
        erro_padrao = np.sqrt(error[1]**2 + err_esperado**2)

        first = ('A discrepancia do valor obtido com o valor esperado e: %6.4f \n' %(discrepancia))
        second = ('A raiz da soma quadratica do erro do valor obtido com o erro do valor esperado e: %6.4f \n' %(erro_padrao))
        third = ('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f \n' % (test))
        fourth = ("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pela discrepancia: \n")
        if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao):
            fifth = ("O fit e bom!")
        else:
            fifth = ("O fit e ruim! Tente outros valores.")
        message2 = str(first +"\n"+ second +"\n"+ third +"\n"+ fourth +"\n"+ fifth)
    return fit, message2
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def apagar():
    canvas = Canvas(janela).grid(row=7, column=0, rowspan=15, columnspan=6, stick=N+S+E+W)
    canvas2 = Canvas(janela, width=120, height=30).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Verificar(): # Método para verificar qual pico foi selecionado.
    global low_lim
    global up_lim
    global bins
    
    if (c1.get() == 1 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Rho
        low_lim = -0.15
        up_lim = -0.05
        bins = 50
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 1 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Phi
        low_lim = -0.09
        up_lim = 0.09
        bins = 50
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 1 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico J/Psi
        low_lim = 0.47
        up_lim = 0.52
        bins = 50
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 1 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Psi'
        low_lim = 0.51
        up_lim = 0.62
        bins = 50
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 1 and c6.get() == 0):
        # Limite para o pico Upsilon
        low_lim = 0.92
        up_lim = 1.02
        bins = 50
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 1):
        # Limite para o pico Z
        low_lim = 1.85
        up_lim = 2.05
        bins = 100
        Escolher_fit()
    else:
        canvas = Canvas(janela).grid(row=5, column=0, rowspan=15, columnspan=6, stick=N+S+E+W)
        canvas2 = Canvas(janela, width=120, height=30).grid(row=0, column=1)
        Carregar_img("Hist0.png")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Abrir():
    global ds
    
    file =  filedialog.askopenfilename(initialdir = "/",title = "Abrir arquivo",filetypes = (("csv files","*.csv"),("CSV files", "*.csv")))
    ds = pd.read_csv(file, engine='python') #Lê o arquivo com os dados.

    Imp_Graf()
        
    Label(janela, text= "Escolha o pico para calcular o ajuste:").grid(row=1, column=0, columnspan=3)
    
    Checkbutton(janela, text= u"\u03c1", variable = c1, width=12, command=Verificar).grid( row=2, column=0)     # u"\u03c1" = código unicode para a letra grega Rho.
    Checkbutton(janela, text= u"\u03a6", variable = c2, width=12, command=Verificar).grid( row=2, column=1)     # u"\u03a6" = código unicode para a letra grega Phi.
    Checkbutton(janela, text= u"J/\u03c8", variable = c3, width=12, command=Verificar).grid(row=2, column=2)   # u"J/\u03c8" = código unicode para a letra grega J/Psi.
    Checkbutton(janela, text= u"\u03c8'", variable = c4, width=12, command=Verificar).grid(row=3, column=0)    # u"\u03c8'" = código unicode para a letra grega Psi'.
    Checkbutton(janela, text= u"\u03a5", variable = c5, width=12, command=Verificar).grid(row=3, column=1)     # u"\u03a5" = código unicode para a letra grega Upsilon.
    Checkbutton(janela, text= "Z", variable = c6, width=12, command=Verificar).grid(row=3, column=2)
    Checkbutton(janela, text= "Todos (Estamos trabalhando nisto!)", state = DISABLED).grid(row=4, column=0, columnspan=3)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Escolher_fit():
    global varfit
    global options
    Label(janela, text="Escolha qual distribuição usar para o ajuste: ").grid(row=5, column=0, columnspan=3)
    options = ["Breit-Wigner",
               "Duas Gaussianas + Exponencial",
               "CrystalBall + Exponencial"]
    varfit = StringVar(janela)
    OptionMenu(janela, varfit, options[0], options[1], options[2]).grid(row=6, column=1)
    varfit.trace("w",Escolher_dist)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Escolher_dist(*args):
    valor = varfit.get()
    if (valor == options[0]):
        apagar()
        Entrar_dados1() # Breit-Wigner
    elif (valor == options[1]):
        apagar()
        Entrar_dados2() # Duas Gaussianas + Exponencial
    elif (valor == options[2]):
        apagar()
        Entrar_dados3() # Crystalball + Exponencial
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Imp_Graf():
    global x
    global y
    global inv_mass_log
    global bins
    
    inv_mass_log = np.log10(np.sqrt(2*ds.pt1*ds.pt2*(np.cosh(ds.eta1-ds.eta2)-np.cos(ds.phi1-ds.phi2))))
    weights = []
    for i in ds["M"]:
        weights.append(bins/np.log(10)/i)
    
    plt.title('Histograma da massa invariante de dois múons \n')
    histogram = plt.hist(inv_mass_log, bins=bins, range=(low_lim, up_lim), weights=weights, color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    plt.yscale('log')
    plt.xlabel('log10(Massa invariante) [log10(GeV)]')
    plt.ylabel('Número de eventos [log10]')
    plt.annotate(r'$\rho$', xy=(-0.11,300000))
    plt.annotate(r'$\Phi$', xy=(0.010,300000))
    plt.annotate(r'$J/\psi$', xy=(0.51,650000))
    plt.annotate(r'$\psi$´', xy=(0.57,100000))
    plt.annotate(r'$\Upsilon$', xy=(0.975,200000))
    plt.annotate('Z', xy=(1.965,30000))
    plt.rcParams['figure.figsize'] = (8,6)  # Redimensiona a o gráfico padrão para o tamanho fixo.
    plt.savefig('Hist0.png') #Salva o gráfico como uma imagem PNG.

    Carregar_img("Hist0.png")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados1(): # Breit-Wigner - 5 graus de liberdade
    Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=7, column=0, columnspan=5)
    e1.set("4")
    Entry(janela, textvar=e1).grid(row=7, column=5)

    Label(janela, text= "Entre com o valor da posição do máximo da distribuição (M):").grid(row=8, column=0, columnspan=5)
    e2.set("91")
    Entry(janela, textvar=e2).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor do parâmtro que é a inclinação usada para perceber o efeito do background:").grid(row=9, column=0, columnspan=5)
    e3.set("-2")
    Entry(janela, textvar=e3).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor do parâmetro b de intercepção no eixo y usada para perceber o efeito do background:").grid(row=10, column=0, columnspan=5)
    e4.set("150")
    Entry(janela, textvar=e4).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor da amplitude da distribuição de Breit-Wigner:").grid(row=11, column=0, columnspan=5)
    e5.set("13000")
    Entry(janela, textvar=e5).grid(row=11, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular1).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados2(): # Duas Gaussianas + Exponencial - 8 graus de liberdade
    Label(janela, text= "Entre com o valor da altura referente a primeira Gaussiana:").grid(row=7, column=0, columnspan=5)
    e1.set("1000")
    Entry(janela, textvar=e1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor da média referente a primeira Gaussiana:").grid(row=8, column=0, columnspan=5)
    e2.set("9.4")
    Entry(janela, textvar=e2).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor do desvio padrão referente a primeira Gaussiana:").grid(row=9, column=0, columnspan=5)
    e3.set("0.1")
    Entry(janela, textvar=e3).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor da altura referente a segunda Gaussiana:").grid(row=10, column=0, columnspan=5)
    e4.set("400")
    Entry(janela, textvar=e4).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor da média referente a segunda Gaussiana:").grid(row=11, column=0, columnspan=5)
    e5.set("10")
    Entry(janela, textvar=e5).grid(row=11, column=5)

    Label(janela, text= "Entre com o valor do desvio padrão referente a segunda Gaussiana:").grid(row=12, column=0, columnspan=5)
    e6.set("0.1")
    Entry(janela, textvar=e6).grid(row=12, column=5)

    Label(janela, text= "Entre com o valor da constante:").grid(row=13, column=0, columnspan=5)
    e7.set("0")
    Entry(janela, textvar=e7).grid(row=13, column=5)

    Label(janela, text= "Entre com o valor da inclinação:").grid(row=14, column=0, columnspan=5)
    e8.set("1")
    Entry(janela, textvar=e8).grid(row=14, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular2).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados3(): # Crystalball + Exponencial - 6 graus de liberdade
    Label(janela, text= "Entre com o valor que define como a função decresce no pico:").grid(row=7, column=0, columnspan=5)
    e1.set("1")
    Entry(janela, textvar=e1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor de n:").grid(row=8, column=0, columnspan=5)
    e2.set("0")
    Entry(janela, textvar=e2).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor da posição do centro do pico:").grid(row=9, column=0, columnspan=5)
    e3.set("3.7")
    Entry(janela, textvar=e3).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor do desvio padrão referente ao pico:").grid(row=10, column=0, columnspan=5)
    e4.set("1")
    Entry(janela, textvar=e4).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor da constante:").grid(row=11, column=0, columnspan=5)
    e5.set("0")
    Entry(janela, textvar=e5).grid(row=11, column=5)
    
    Label(janela, text= "Entre com o valor da inclinação:").grid(row=12, column=0, columnspan=5)
    e6.set("-1")
    Entry(janela, textvar=e6).grid(row=12, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular3).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular1(): # Breit-Wigner
    global bins

    a1 = float(e1.get())
    a2 = float(e2.get())
    a3 = float(e3.get())
    a4 = float(e4.get())
    a5 = float(e5.get())
    
    plt.clf()
    inv_mass = np.sqrt(2*ds.pt1*ds.pt2*(np.cosh(ds.eta1-ds.eta2)-np.cos(ds.phi1-ds.phi2)))
    newlowlim = np.power(10,low_lim)
    newuplim  = np.power(10,up_lim)
    limitedmasses = inv_mass[(inv_mass > newlowlim ) & (inv_mass < newuplim)]
    histogram = plt.hist(limitedmasses, bins=bins, range=(newlowlim,newuplim), color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    
    initials = [a1, a2, a3, a4, a5]
    best, covariance = curve_fit(f.breitwigner, x, y, p0=initials, sigma=np.sqrt(y))
    error = np.sqrt(np.diag(covariance))
    plt.plot(x, f.breitwigner(x, *best), 'r-')
    plt.xlabel('Massa Invariante  [GeV]')
    plt.ylabel('Número de eventos ')
    plt.title('Ajuste de Breit-Wigner')
    plt.rcParams['figure.figsize'] = (8,6)
    plt.savefig("Hist1.png")
    Carregar_img("Hist1.png")
    
    # Diferença entre os valores iniciais e o melhor valor após o 1º curve_fit.
    dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4])]
    # Iteração para convergir para o melhor valor dos parâmetros.
    i = int
    i = 0
    while (dif[0] > 0 and dif[1] > 0 and dif[2] > 0 and dif[3] > 0 and dif[4] > 0 and i <= 20):
        initials = [best[0], best[1], best[2], best[3], best[4]]
        best, covariance = curve_fit(f.breitwigner, x, y, p0=initials, sigma=np.sqrt(y))
        error = np.sqrt(np.diag(covariance))
        # Atualização do valor da diferença, para que haja convergência da iteração.
        dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4])]
        i += 1
        iteration = str("Quantidade de iteração: {:1.1f} \n".format(i))
    
    first = "Valor da largura de decaimento (gamma) = {:4.4f} +- {:4.4f} \n".format(best[0], error[0])
    second = "Valor do pico da distribuicao (M) = {:4.4f} +- {:4.4f} \n".format(best[1], error[1])
    third = "a = {:4.4f} +- {:4.4f} \n".format(best[2], error[2])
    fourth = "b = {:4.4f} +- {:4.4f} \n".format(best[3], error[3])
    fifth = "A = {:4.4f} +- {:4.4f} \n".format(best[4], error[4])
    message1 = str(first +"\n"+ second +"\n"+ third +"\n"+ fourth +"\n"+ fifth)
    
    fit, message2 = convergencia(i, f.breitwigner, best, error, 5)
    
    message = (message1+"\n"+iteration+"\n"+fit+"\n"+message2)
    messagebox.showinfo("ATENÇÃO!", message)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular2(): # Duas Gaussianas + Exponencial
    global bins
    a1 = float(e1.get())
    a2 = float(e2.get())
    a3 = float(e3.get())
    a4 = float(e4.get())
    a5 = float(e5.get())
    a6 = float(e6.get())
    a7 = float(e7.get())
    a8 = float(e8.get())
    plt.clf()
    inv_mass = np.sqrt(2*ds.pt1*ds.pt2*(np.cosh(ds.eta1-ds.eta2)-np.cos(ds.phi1-ds.phi2)))
    newlowlim = np.power(10,low_lim)
    newuplim  = np.power(10,up_lim)
    limitedmasses = inv_mass[(inv_mass > newlowlim ) & (inv_mass < newuplim)]
    histogram = plt.hist(limitedmasses, bins=bins, range=(newlowlim,newuplim), color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    
    initials = [a1, a2, a3, a4, a5, a6, a7, a8]
    best, covariance = curve_fit(f.doublegaussianexpo, x, y, p0=initials, sigma=np.sqrt(y))
    error = np.sqrt(np.diag(covariance))
    plt.plot(x, f.doublegaussianexpo(x, *best), 'r-')
    plt.xlabel('Massa Invariante  [GeV]')
    plt.ylabel('Número de eventos ')
    plt.title('Ajuste de Dupla Gaussiana + Exponencial')
    plt.rcParams['figure.figsize'] = (8,6)
    plt.savefig("Hist2.png")
    Carregar_img("Hist2.png")
        
    # Diferença entre os valores iniciais e o melhor valor após o 1º curve_fit.
    dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4]), np.absolute(best[5] - initials[5]), np.absolute(best[6] - initials[6]), np.absolute(best[7] - initials[7])]
    # Iteração para convergir para o melhor valor dos parâmetros.
    i = int
    i = 0
    while (dif[0] > 0 and dif[1] > 0 and dif[2] > 0 and dif[3] > 0 and dif[4] > 0 and dif[5] > 0 and dif[6] > 0 and dif[7] > 0 and i <= 20):
        initials = [best[0], best[1], best[2], best[3], best[4], best[5], best[6], best[7]]
        best, covariance = curve_fit(f.doublegaussianexpo, x, y, p0=initials, sigma=np.sqrt(y))
        error = np.sqrt(np.diag(covariance))
        # Atualização do valor da diferença, para que haja convergência da iteração.
        dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4]), np.absolute(best[5] - initials[5]), np.absolute(best[6] - initials[6]), np.absolute(best[7] - initials[7])]
        i += 1
        iteration = str("Quantidade de iteração: {:1.1f} \n".format(i))
    
    first = "Valor de a1 = {:4.4f} +- {:4.4f} \n".format(best[0], error[0])
    second = "Valor de mean1 = {:4.4f} +- {:4.4f} \n".format(best[1], error[1])
    third = "Valor de sigma1 = {:4.4f} +- {:4.4f} \n".format(best[2], error[2])
    fourth = "Valor de a2 = {:4.4f} +- {:4.4f} \n".format(best[3], error[3])
    fifth = "Valor de mean2 = {:4.4f} +- {:4.4f} \n".format(best[4], error[4])
    sixth = "Valor de sigma2 = {:4.4f} +- {:4.4f} \n".format(best[5], error[5])
    sixth = "Valor da constante = {:4.4f} +- {:4.4f} \n".format(best[6], error[6])
    seventh = "Valor da inclinação = {:4.4f} +- {:4.4f} \n".format(best[7], error[7])
    message1 = str(first +"\n"+ second +"\n"+ third +"\n"+ fourth +"\n"+ fifth +"\n"+ sixth +"\n"+ seventh)
    
    fit, message2 = convergencia(i, f.doublegaussianexpo, best, error, 8)
    
    message = (message1+"\n"+iteration+"\n"+fit+"\n"+message2)
    messagebox.showinfo("ATENÇÃO!", message)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular3(): # Crystal Ball + Exponencial
    global bins
    a1 = float(e1.get())
    a2 = float(e2.get())
    a3 = float(e3.get())
    a4 = float(e4.get())
    a5 = float(e5.get())
    a6 = float(e6.get())
    
    plt.clf()
    inv_mass = np.sqrt(2*ds.pt1*ds.pt2*(np.cosh(ds.eta1-ds.eta2)-np.cos(ds.phi1-ds.phi2)))
    newlowlim = np.power(10,low_lim)
    newuplim  = np.power(10,up_lim)
    limitedmasses = inv_mass[(inv_mass > newlowlim ) & (inv_mass < newuplim)]
    histogram = plt.hist(limitedmasses, bins=bins, range=(newlowlim,newuplim), color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    
    initials = [a1, a2, a3, a4, a5, a6]
    best, covariance = curve_fit(f.crystalexpo, x, y, p0=initials, sigma=np.sqrt(y), maxfev=50000)
    error = np.sqrt(np.diag(covariance))
    plt.plot(x, f.crystalexpo(x, *best), 'r-')
    plt.xlabel('Massa Invariante  [GeV]')
    plt.ylabel('Número de eventos ')
    plt.title('Ajuste de Crystal_Ball + Exponencial')
    plt.rcParams['figure.figsize'] = (8,6)
    plt.savefig("Hist3.png")
    Carregar_img("Hist3.png")
    
    # Diferença entre os valores iniciais e o melhor valor após o 1º curve_fit.
    dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4]), np.absolute(best[5] - initials[5])]
    # Iteração para convergir para o melhor valor dos parâmetros.
    i = int
    i = 0
    while (dif[0] > 0 and dif[1] > 0 and dif[2] > 0 and dif[3] > 0 and dif[4] > 0 and dif[5] and i <= 20):
        initials = [best[0], best[1], best[2], best[3], best[4], best[5]]
        best, covariance = curve_fit(f.crystalexpo, x, y, p0=initials, sigma=np.sqrt(y))
        error = np.sqrt(np.diag(covariance))
        # Atualização do valor da diferença, para que haja convergência da iteração.
        dif = [np.absolute(best[0] - initials[0]), np.absolute(best[1] - initials[1]), np.absolute(best[2] - initials[2]), np.absolute(best[3] - initials[3]), np.absolute(best[4] - initials[4]), np.absolute(best[5] - initials[5])]
        i += 1
        iteration = str("Quantidade de iteração: {:1.1f} \n".format(i))
    
    first = "Valor de a= {:4.4f} +- {:4.4f} \n".format(best[0], error[0])
    second = "Valor de n = {:4.4f} +- {:4.4f} \n".format(best[1], error[1])
    third = "Valor de mean = {:4.4f} +- {:4.4f} \n".format(best[2], error[2])
    fourth = "Valor de sigma = {:4.4f} +- {:4.4f} \n".format(best[3], error[3])
    fifth = "Valor da constante = {:4.4f} +- {:4.4f} \n".format(best[4], error[4])
    sixth = "Valor da inclinação = {:4.4f} +- {:4.4f} \n".format(best[5], error[5])
    
    message1 = str(first +"\n"+ second +"\n"+ third +"\n"+ fourth +"\n"+ fifth +"\n"+ sixth)
    
    fit, message2 = convergencia(i, f.crystalexpo, best, error, 6)
    
    message = (message1+"\n"+iteration+"\n"+fit+"\n"+message2)
    messagebox.showinfo("ATENÇÃO!", message)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
Button(janela, text='Abrir arquivo', font=('Times', '12', 'bold'), width = 12, command=Abrir).grid(row=0, column=0)

#janela.geometry("300x150")
janela.mainloop() #mantém a janela aberta