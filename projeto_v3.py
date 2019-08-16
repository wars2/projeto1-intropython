# Importa todos os módulos necessários para a execução do programa
from tkinter import *
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageTk
from scipy.optimize import curve_fit
from scipy.stats import norm
from tkinter import messagebox

ds = pd.DataFrame()
x = ()
y = ()
inv_mass_log = ()
################################################################################################################################
############################################### DEFINIÇÕES DAS FUNÇÕES DOS AJUSTES #############################################
################################################################################################################################
def expo(x, const, slope):
    """ Uma curva exponencial. 
  parametros: const, slope """
    return np.exp(const + slope*x)

def line(x, intercept, slope):
    """ Polinômio do primeiro grau. """
    return slope*x + intercept

def breitwigner(E, gamma, M, a, b, A):
    ''' E (é a energia)
        gamma (a largura total do meio no máximo da distribuição)
        M (valor onde ocorre o máximo da distribuição)
        a (inclinação que é usada para pereber o efeito de backgrund)
        b (intercepção em y, que é usada para perceber o efeito de background)
        A ("amplitude" da distribuição de Breit-Wigner) '''
    return a*E+b+A*( (2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))) )/((E**2-M**2)**2+M**2*gamma**2)

def gaussian(x, a, x0, sigma):
    ''' a (altura do pico)
        x0 (média, ordena a posição do centro do pico)
        sigma (desvio padrão, controla a largura da curva)
    '''
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def doublegaussian(x, a1, x01, sigma1, a2, x02, sigma2):
    '''Duas gaussianas somadas.
    '''
    return a1*np.exp(-(x-x01)**2/(2*sigma1**2))+a2*np.exp(-(x-x02)**2/(2*sigma2**2))

def crystalball(x, a, n, xb, sig):
    x = x+0j
    #N, a, n, xb, sig = parametros
    if a < 0:
        a = -a
    if n < 0:
        n = -n
    aa = abs(a)
    A = (n/aa)**n * np.exp(- aa**2 / 2)
    B = n/aa - aa
    C = n/aa / (n-1.) * np.exp(-aa**2/2.)
    D = np.sqrt(pi/2.) * (1. + erf(aa/np.sqrt(2.)))
    N = 1. / (sig * (C+D))
    total = 0.*x
    total += ((x-xb)/sig  > -a) * N * np.exp(- (x-xb)**2/(2.*sig**2))
    total += ((x-xb)/sig <= -a) * N * A * (B - (x-xb)/sig)**(-n)
    try:
      return total.real

    except:
      return total
    return total

def crystalexpo(x, a, n, xb, sig, const, slope):#6 parâmetros
    ''' Crystal-ball somada com exponencial.
    '''
    return (crystalball(x, a, n, xb, sig)+expo(x, const, slope))

def doublecrystal(x, a1, n1, xb1, sig1, a2, n2, xb2, sig2): #8 parâmetros
    ''' Duas crystal-ball somadas.
    '''
    return (crystalball(x, a1, n1, xb1, sig1)+crystalball(x, a2, n2, xb2, sig2))

def gaussianexpo(x, a, x0, sigma, const, slope): #5 parâmetros
    ''' Gaussiana somada com exponencial.
    '''
    return(gaussian(x, a, x0, sigma)+expo(x, const, slope))

def doublegaussianexpo(x, a1, x01, sigma1, a2, x02, sigma2, const, slope): #8 parâmetros
    ''' Soma de duas guassianas e uma exponencial.
    '''
    return (doublegaussian(x, a1, x01, sigma1, a2, x02, sigma2)+expo(x, const, slope))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
'''def Calcular_Breitwigner:
    print(" gamma (a largura total do meio no maximo da distribuicao),\n M(valor onde ocorre o maximo da distribuicao),\n a (inclinacao que e usada para pereber o efeito de backgrund),\n b (intercepcao em y, que e usada para perceber o efeito de background),\n A (amplitude da distribuicao de Breit-Wigner). \nRecomendamos os valores: 4 91 -2 150 13000 para um fit bem sucedido.")
    initials_z =  [float(x) for x in input('Preencha com os parametros da Breit-Wigner (gamma, M, a, b, A) dando apenas um espaco entre eles: ').split()]
    # Vamos importar o módulo que é usado na otimização, executar a otimização e calcular as incertezas dos parâmetros otimizados.
    best_z, covariance = curve_fit(breitwigner, x, y, p0=initials_z, sigma=np.sqrt(y))
    error_z = np.sqrt(np.diag(covariance))

    # Vamos imprimir os valores e incertezas obtidos com a otimização.
    print("Valores com incertezas")
    print("")
    first = "Valor da largura de decaimento (gamma) = {} +- {}".format(best_z[0], error_z[0])
    second = "Valor do pico da distribuicao (M) = {} +- {}".format(best_z[1], error_z[1])
    third = "a = {} +- {}".format(best_z[2], error_z[2])
    fourth = "b = {} +- {}".format(best_z[3], error_z[3])
    fifth = "A = {} +- {}".format(best_z[4], error_z[4])
    print(first)
    print(second)
    print(third)
    print(fourth)
    print(fifth)

    # Diferença entre os valores iniciais e o melhor valor após o 1º curve_fit.
    dif_z = [np.absolute(best_z[0] - initials_z[0]), np.absolute(best_z[1] - initials_z[1]), np.absolute(best_z[2] - initials_z[2]), np.absolute(best_z[3] - initials_z[3]), np.absolute(best_z[4] - initials_z[4])]

    # Iteração para convergir para o melhor valor dos parâmetros.
    while (dif_z[0] > 0 and dif_z[1] > 0 and dif_z[2] > 0 and dif_z[3] > 0 and dif_z[4] > 0 and i <= 20):
        initials_z = [best_z[0], best_z[1], best_z[2], best_z[3], best_z[4]]
        best_z, covariance = curve_fit(breitwigner, x, y, p0=initials_z, sigma=np.sqrt(y))
        error_z = np.sqrt(np.diag(covariance))
        first = "Valor da largura de decaimento (gamma) = {} +- {}".format(best_z[0], error_z[0])
        second = "Valor do pico da distribuicao (M) = {} +- {}".format(best_z[1], error_z[1])
        third = "a = {} +- {}".format(best_z[2], error_z[2])
        fourth = "b = {} +- {}".format(best_z[3], error_z[3])
        fifth = "A = {} +- {}".format(best_z[4], error_z[4])
        print(first)
        print(second)
        print(third)
        print(fourth)
        print(fifth)
        # Atualização do valor da diferença, para que haja convergência da iteração.
        dif_z = [np.absolute(best_z[0] - initials_z[0]), np.absolute(best_z[1] - initials_z[1]), np.absolute(best_z[2] - initials_z[2]), np.absolute(best_z[3] - initials_z[3]), np.absolute(best_z[4] - initials_z[4])]
        i += 1
        print("Iteracao numero: ", i)

    # Critério de avaliação dos ajustes: Iterações, chi2 normalizado e compatibilidade.
    print("Numero de iteracoes: ", i)
    print("Baseado nos numeros de iteracoes: ")
    if (i == 0 and i == 1):
        print("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! ")
    elif (i >= 1 and i <= 19):
        print("O fit convergiu!")
    else:
        print("O fit provavelmente esta divergindo... Tente outros valores iniciais!")
    
    ch2, pval = chisquare(y, breitwigner(x, *best_z))
    test = ch2 /(bins - 5)
    discrepancia = np.absolute(best_z[1] - esperado)
    erro_padrao = np.sqrt(error_z[1]**2 + err_esperado**2)
    
    print('A discrepancia do valor obtido com o valor esperado e: %6.4f' %(discrepancia))
    print('O erro padrao do erro do valor obtido com o erro do valor esperado e: %6.4f' %(erro_padrao))
    print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
    print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pela discrepancia: ")
    if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao):
        print ("O fit e bom!")
    else:
        print("O fit e ruim! Tente outros valores.")
    
    # Plotando o ajuste.
    plt.plot(x, breitwigner(x, *best_z), 'r-')
    plt.xlabel('Massa Invariante [GeV]')
    plt.ylabel('Numero de Eventos')
    plt.title('Boson Z: Ajuste com Breit-Wigner')
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular_DuasGaussianasExpo:
    print(" a (altura do pico), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1000 9.4 0.1 400 10 0.1 0 1 para um fit bem sucedido.")
    initials_upsilon =  [float(x) for x in input('Preencha com os parametros das duas Gaussianas e da exponencial (a1, mean1, sigma1, a2, mean2, sigma2, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
    best_upsilon, covariance = curve_fit(doublegaussianexpo, x, y, p0=initials_upsilon, sigma=np.sqrt(y))
    error_upsilon = np.sqrt(np.diag(covariance))

    print("Valores com incertezas")
    print("")
    primeiro = "Valor de a1 = {} +- {}".format(best_upsilon[0], error_upsilon[0])
    segundo = "Valor de mean1 = {} +- {}".format(best_upsilon[1], error_upsilon[1])
    terceiro = "Valor de sigma1 = {} +- {}".format(best_upsilon[2], error_upsilon[2])
    quarto = "Valor de a2 = {} +- {}".format(best_upsilon[3], error_upsilon[3])
    quinto = "Valor de mean2 = {} +- {}".format(best_upsilon[4], error_upsilon[4])
    sexto = "Valor de sigma2 = {} +- {}".format(best_upsilon[5], error_upsilon[5])
    setimo = "Valor de constante = {} +- {}".format(best_upsilon[6], error_upsilon[6])
    oitavo = "Valor da inclinacao = {} +- {}".format(best_upsilon[7], error_upsilon[7])
    print(primeiro)
    print(segundo)
    print(terceiro)
    print(quarto)
    print(quinto)
    print(sexto)
    print(setimo)
    print(oitavo)

    dif_upsilon = [np.absolute(best_upsilon[0] - initials_upsilon[0]), np.absolute(best_upsilon[1] - initials_upsilon[1]), np.absolute(best_upsilon[2] - initials_upsilon[2]), np.absolute(best_upsilon[3] - initials_upsilon[3]), np.absolute(best_upsilon[4] - initials_upsilon[4]), np.absolute(best_upsilon[5] - initials_upsilon[5]), np.absolute(best_upsilon[6] - initials_upsilon[6]), np.absolute(best_upsilon[7] - initials_upsilon[7])]

    while (dif_upsilon[0] > 0 and dif_upsilon[1] > 0 and dif_upsilon[2] > 0 and dif_upsilon[3] > 0 and dif_upsilon[4] > 0 and dif_upsilon[5] > 0 and dif_upsilon[6] > 0 and dif_upsilon[7] > 0 and i <= 20):
        initials_upsilon = [best_upsilon[0], best_upsilon[1], best_upsilon[2], best_upsilon[3], best_upsilon[4], best_upsilon[5], best_upsilon[6], best_upsilon[7]]
        best_upsilon, covariance = curve_fit(doublegaussianexpo, x, y, p0=initials_upsilon, sigma=np.sqrt(y))
        error_upsilon = np.sqrt(np.diag(covariance))
        primeiro = "Valor de a1 = {} +- {}".format(best_upsilon[0], error_upsilon[0])
        segundo = "Valor de mean1 = {} +- {}".format(best_upsilon[1], error_upsilon[1])
        terceiro = "Valor de sigma1 = {} +- {}".format(best_upsilon[2], error_upsilon[2])
        quarto = "Valor de a2 = {} +- {}".format(best_upsilon[3], error_upsilon[3])
        quinto = "Valor de mean2 = {} +- {}".format(best_upsilon[4], error_upsilon[4])
        sexto = "Valor de sigma2 = {} +- {}".format(best_upsilon[5], error_upsilon[5])
        setimo = "Valor de constante = {} +- {}".format(best_upsilon[6], error_upsilon[6])
        oitavo = "Valor da inclinacao = {} +- {}".format(best_upsilon[7], error_upsilon[7])
        print(primeiro)
        print(segundo)
        print(terceiro)
        print(quarto)
        print(quinto)
        print(sexto)
        print(setimo)
        print(oitavo)
        dif_upsilon = [np.absolute(best_upsilon[0] - initials_upsilon[0]), np.absolute(best_upsilon[1] - initials_upsilon[1]), np.absolute(best_upsilon[2] - initials_upsilon[2]), np.absolute(best_upsilon[3] - initials_upsilon[3]), np.absolute(best_upsilon[4] - initials_upsilon[4]), np.absolute(best_upsilon[5] - initials_upsilon[5]), np.absolute(best_upsilon[6] - initials_upsilon[6]), np.absolute(best_upsilon[7] - initials_upsilon[7])]
        i += 1
        print("Iteracao numero: ", i)

    print("Numero de iteracoes: ", i)
    print("Baseado nos numeros de iteracoes: ")
    if (i == 1):
        print("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! ")
    elif (i >= 1 and i <= 19):
        print("O fit convergiu!")
    else:
        print("O fit provavelmente esta divergindo... Tente outros valores iniciais!")
    
    ch2, pval = chisquare(y, doublegaussianexpo(x, *best_upsilon))
    test = ch2 /(bins - 8)
    discrepancia = np.absolute(best_upsilon[1] - esperado)
    discrepancia1 = np.absolute(best_upsilon[4] - esperado1)
    erro_padrao = np.sqrt(error_upsilon[1]**2 + err_esperado**2)
    erro_padrao1 = np.sqrt(error_upsilon[4]**2 + err_esperado1**2)
    
    print('A discrepancia do valor obtido no primeiro pico com o valor esperado e: %6.4f' %(discrepancia))
    print('A discrepancia do valor obtido no segundo pico com o valor esperado e: %6.4f' %(discrepancia1))
    print('O erro padrao do erro do valor obtido no primeiro pico com o erro do valor esperado e: %6.4f' %(erro_padrao))
    print('O erro padrao do erro do valor obtido no segundo pico com o erro do valor esperado e: %6.4f' %(erro_padrao1))
    print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
    print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pelas discrepancias: ")
    if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao and discrepancia < 2*erro_padrao1):
        print ("O fit e bom!")
    else:
        print("O fit e ruim! Tente outros valores.")
    
    plt.plot(x, doublegaussianexpo(x, *best_upsilon), 'r-')
    plt.xlabel('Massa Invariante [GeV]')
    plt.ylabel('Numero de Eventos')
    plt.title('Upsilon: Ajuste com Duas Gaussianas + Exponencial')
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular_GaussianaExpo:
    print(" a (define como a funcao decresce no pico), \n n (), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1.6 0.9 3.1 1 0 -1 para um fit bem sucedido.")
    initials_jpsi =  [float(x) for x in input('Preencha com os parametros da Crystal-Ball e da exponencial (a, n, mean, sigma, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
    best_jpsi, covariance = curve_fit(crystalexpo, x, y, p0=initials_jpsi, sigma=np.sqrt(y))
    error_jpsi = np.sqrt(np.diag(covariance))

    print("Valores com incertezas")
    print("")
    primeiro = "Valor de a = {} +- {}".format(best_jpsi[0], error_jpsi[0])
    segundo = "Valor de n = {} +- {}".format(best_jpsi[1], error_jpsi[1])
    terceiro = "Valor de mean = {} +- {}".format(best_jpsi[2], error_jpsi[2])
    quarto = "Valor de sigma = {} +- {}".format(best_jpsi[3], error_jpsi[3])
    quinto = "Valor da constante = {} +- {}".format(best_jpsi[4], error_jpsi[4])
    sexto = "Valor da inclinacao = {} +- {}".format(best_jpsi[5], error_jpsi[5])
    print(primeiro)
    print(segundo)
    print(terceiro)
    print(quarto)
    print(quinto)
    print(sexto)
    dif_jpsi = [np.absolute(best_jpsi[0] - initials_jpsi[0]), np.absolute(best_jpsi[1] - initials_jpsi[1]), np.absolute(best_jpsi[2] - initials_jpsi[2]), np.absolute(best_jpsi[3] - initials_jpsi[3]), np.absolute(best_jpsi[4] - initials_jpsi[4]), np.absolute(best_jpsi[5] - initials_jpsi[5])]

    while (dif_jpsi[0] > 0 and dif_jpsi[1] > 0 and dif_jpsi[2] > 0 and dif_jpsi[3] > 0 and dif_jpsi[4] > 0 and dif_jpsi[5] and i <= 20):
        initials_jpsi = [best_jpsi[0], best_jpsi[1], best_jpsi[2], best_jpsi[3], best_jpsi[4], best_jpsi[5]]
        best_jpsi, covariance = curve_fit(crystalexpo, x, y, p0=initials_jpsi, sigma=np.sqrt(y))
        error_jpsi = np.sqrt(np.diag(covariance))
        primeiro = "Valor de a = {} +- {}".format(best_jpsi[0], error_jpsi[0])
        segundo = "Valor de n = {} +- {}".format(best_jpsi[1], error_jpsi[1])
        terceiro = "Valor de mean = {} +- {}".format(best_jpsi[2], error_jpsi[2])
        quarto = "Valor de sigma = {} +- {}".format(best_jpsi[3], error_jpsi[3])
        quinto = "Valor da constante = {} +- {}".format(best_jpsi[4], error_jpsi[4])
        sexto = "Valor da inclinacao = {} +- {}".format(best_jpsi[5], error_jpsi[5])
        print(primeiro)
        print(segundo)
        print(terceiro)
        print(quarto)
        print(quinto)
        print(sexto)
        dif_jpsi = [np.absolute(best_jpsi[0] - initials_jpsi[0]), np.absolute(best_jpsi[1] - initials_jpsi[1]), np.absolute(best_jpsi[2] - initials_jpsi[2]), np.absolute(best_jpsi[3] - initials_jpsi[3]), np.absolute(best_jpsi[4] - initials_jpsi[4]), np.absolute(best_jpsi[5] - initials_jpsi[5])]
        i += 1
        print("Iteracao numero: ", i)

    print("Numero de iteracoes: ", i)
    print("Baseado nos numeros de iteracoes: ")
    if (i == 1):
        print("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! ")
    elif (i >= 1 and i <= 19):
        print("O fit convergiu!")
    else:
        print("O fit provavelmente esta divergindo... Tente outros valores iniciais!")
    
    ch2, pval = chisquare(y, crystalexpo(x, *best_jpsi))
    test = ch2 /(bins - 6)
    discrepancia = np.absolute(best_jpsi[2] - esperado)
    erro_padrao = np.sqrt(error_jpsi[2]**2 + err_esperado**2)
    
    print('A discrepancia do valor obtido com o valor esperado e: %6.4f' %(discrepancia))
    print('O erro padrao do erro do valor obtido com o erro do valor esperado e: %6.4f' %(erro_padrao))
    print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
    print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pela discrepancia: ")
    if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao):
        print ("O fit e bom!")
    else:
        print("O fit e ruim! Tente outros valores.")
    
    plt.plot(x, crystalexpo(x, *best_jpsi), 'r-')
    plt.xlabel('Massa Invariante [GeV]')
    plt.ylabel('Numero de Eventos')
    plt.title('J/Psi: Ajuste com Gaussian + Exponencial')
    plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular_CrystalballExpo:
    print(" a (define como a funcao decresce no pico), \n n (), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1 0 3.7 1 0 -1 para um fit bem sucedido.")
    initials_psiprime =  [float(x) for x in input('Preencha com os parametros da Crystal-Ball e da exponencial (a, n, mean, sigma, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
    best_psiprime, covariance = curve_fit(crystalexpo, x, y, p0=initials_psiprime, sigma=np.sqrt(y))
    error_psiprime = np.sqrt(np.diag(covariance))

    print("Valores com incertezas")
    print("")
    primeiro = "Valor de a = {} +- {}".format(best_psiprime[0], error_psiprime[0])
    segundo = "Valor de n = {} +- {}".format(best_psiprime[1], error_psiprime[1])
    terceiro = "Valor de mean = {} +- {}".format(best_psiprime[2], error_psiprime[2])
    quarto = "Valor de sigma = {} +- {}".format(best_psiprime[3], error_psiprime[3])
    quinto = "Valor da constante = {} +- {}".format(best_psiprime[4], error_psiprime[4])
    sexto = "Valor da inclinacao = {} +- {}".format(best_psiprime[5], error_psiprime[5])
    print(primeiro)
    print(segundo)
    print(terceiro)
    print(quarto)
    print(quinto)
    print(sexto)

    dif_psiprime = [np.absolute(best_psiprime[0] - initials_psiprime[0]), np.absolute(best_psiprime[1] - initials_psiprime[1]), np.absolute(best_psiprime[2] - initials_psiprime[2]), np.absolute(best_psiprime[3] - initials_psiprime[3]), np.absolute(best_psiprime[4] - initials_psiprime[4]), np.absolute(best_psiprime[5] - initials_psiprime[5])]

    while (dif_psiprime[0] > 0 and dif_psiprime[1] > 0 and dif_psiprime[2] > 0 and dif_psiprime[3] > 0 and dif_psiprime[4] > 0 and dif_psiprime[5] > 0 and i <= 20):
        initials_psiprime = [best_psiprime[0], best_psiprime[1], best_psiprime[2], best_psiprime[3], best_psiprime[4], best_psiprime[5]]
        best_psiprime, covariance = curve_fit(crystalexpo, x, y, p0=initials_psiprime, sigma=np.sqrt(y))
        error_psiprime = np.sqrt(np.diag(covariance))
        primeiro = "Valor de a = {} +- {}".format(best_psiprime[0], error_psiprime[0])
        segundo = "Valor de n = {} +- {}".format(best_psiprime[1], error_psiprime[1])
        terceiro = "Valor de mean = {} +- {}".format(best_psiprime[2], error_psiprime[2])
        quarto = "Valor de sigma = {} +- {}".format(best_psiprime[3], error_psiprime[3])
        quinto = "Valor da constante = {} +- {}".format(best_psiprime[4], error_psiprime[4])
        sexto = "Valor da inclinacao = {} +- {}".format(best_psiprime[5], error_psiprime[5])
        print(primeiro)
        print(segundo)
        print(terceiro)
        print(quarto)
        print(quinto)
        print(sexto)
        dif_psiprime = [np.absolute(best_psiprime[0] - initials_psiprime[0]), np.absolute(best_psiprime[1] - initials_psiprime[1]), np.absolute(best_psiprime[2] - initials_psiprime[2]), np.absolute(best_psiprime[3] - initials_psiprime[3]), np.absolute(best_psiprime[4] - initials_psiprime[4]), np.absolute(best_psiprime[5] - initials_psiprime[5])]
        i += 1
        print("Iteracao número: ", i)

    print("Numero de iteracoes: ", i)
    print("Baseado nos numeros de iteracoes: ")
    if (i == 1):
        print("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! ")
    elif (i >= 1 and i <= 19):
        print("O fit convergiu!")
    else:
        print("O fit provavelmente esta divergindo... Tente outros valores iniciais!")
    
    ch2, pval = chisquare(y, crystalexpo(x, *best_psiprime))
    test = ch2 /(bins - 6)
    discrepancia = np.absolute(best_psiprime[2] - esperado)
    erro_padrao = np.sqrt(error_psiprime[2]**2 + err_esperado**2)

    print('A discrepancia do valor obtido com o valor esperado e: %6.4f' %(discrepancia))
    print('O erro padrao do erro do valor obtido com o erro do valor esperado e: %6.4f' %(erro_padrao))
    print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
    print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pela discrepancia: ")
    if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao):
        print ("O fit e bom!")
    else:
        print("O fit e ruim! Tente outros valores.")
    
    plt.plot(x, crystalexpo(x, *best_psiprime), 'r-')
    plt.xlabel('Massa Invariante [GeV]')
    plt.ylabel('Numero de Eventos')
    plt.title('Psi-prime: Ajuste com CrystalBall + Exponencial')
    plt.show()'''
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
# Criação da janela com o tkinter.
janela = Tk() # Cria uma janela
janela.title('Curva de Fit') # Define o título da janela.

a1 = StringVar()
b1 = StringVar()
A1 = StringVar()
gama1 = StringVar()
M1 = StringVar()
sigma1 = StringVar()
media1 = StringVar()
c1 = IntVar()
c2 = IntVar()
c3 = IntVar()
c4 = IntVar()
c5 = IntVar()
c6 = IntVar()
bins = 300
low_lim = -0.5
up_lim = 2.5
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Carregar_img(): # Método para chamar a imagem gerada no primeiro Imp_Graf().
    global img
    global tkimage
    
    img = Image.open("Hist0.png")
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=6, rowspan=20)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def apagar():
    canvas = Canvas(janela).grid(row=7, column=0, rowspan=15, columnspan=6, stick=N+S+E+W)
    canvas2 = Canvas(janela, width=120, height=30).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Verificar(): # Método para verificar qual pico foi selecionado.
    global low_lim
    global up_lim
    
    if (c1.get() == 1 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Rho
        low_lim = -0.15
        up_lim = -0.05
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 1 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Phi
        low_lim = -0.09
        up_lim = 0.09
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 1 and c4.get() == 0 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico J/Psi
        low_lim = 0.47
        up_lim = 0.52
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 1 and c5.get() == 0 and c6.get() == 0):
        # Limite para o pico Psi'
        low_lim = 0.51
        up_lim = 0.62
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 1 and c6.get() == 0):
        # Limite para o pico Upsilon
        low_lim = 0.96
        up_lim = 0.99
        Escolher_fit()
    elif (c1.get() == 0 and c2.get() == 0 and c3.get() == 0 and c4.get() == 0 and c5.get() == 0 and c6.get() == 1):
        # Limite para o pico Z
        low_lim = 1.85
        up_lim = 2.05
        Escolher_fit()
    else:
        canvas = Canvas(janela).grid(row=5, column=0, rowspan=15, columnspan=6, stick=N+S+E+W)
        canvas2 = Canvas(janela, width=120, height=30).grid(row=0, column=1)
        low_lim = -0.5
        up_lim = 2.5
        Carregar_img()
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
    Label(janela, text="Escolha qual distribuição usar para o ajuste: ").grid(row=5, column=0)
    options = ["Breit-Wigner",
               "Duas Gaussianas + Exponencial",
               "Gaussiana + Exponencial",
               "CrystalBall + Exponencial"]
    varfit = StringVar(janela)
    OptionMenu(janela, varfit, options[0], options[1], options[2], options[3]).grid(row=6, column=1)
#    varfit.set(options[0])
    varfit.trace("w",Escolher_dist)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Escolher_dist(*args):
    valor = varfit.get()
    if (valor == options[0]):
        apagar()
        Entrar_dados1() # Breit-Wigner
    elif (varfit.get() == options[1]):
        apagar()
        Entrar_dados2() # Duas Gaussianas + Exponencial
        messagebox.showinfo("Atenção!", "Implementar esta parte!")
    elif (valor == options[2]):
        apagar()
        Entrar_dados3() # Gaussiana + Exponencial
        messagebox.showinfo("Atenção!", "Implementar esta parte!")
    elif (valor == options[3]):
        apagar()
        Entrar_dados4() # Crystalball + Exponencial
        messagebox.showinfo("Atenção!", "Implementar esta parte!")
#    else:
#        messagebox.showinfo("Atenção!", "Selecione uma distribuição para calcular o ajuste!")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Imp_Graf():
    global x
    global y
    global inv_mass_log
    global bins
    
    inv_mass_log = np.log10(ds["M"])
#    bins=300
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
    plt.annotate(r'$J/\psi$', xy=(0.49,700000))
    plt.annotate(r'$\psi$´', xy=(0.57,100000))
    plt.annotate(r'$\Upsilon$', xy=(0.975,200000))
    plt.annotate('Z', xy=(1.965,30000))
    plt.rcParams['figure.figsize'] = (8,6)  # Redimensiona a o gráfico padrão para o tamanho fixo.
    plt.savefig('Hist0.png') #Salva o gráfico como uma imagem PNG.
    
    Carregar_img()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Calcular1():
    gamma = float(gama1.get())
    M = float(M1.get())
    a = float(a1.get())
    b = float(b1.get())
    A = float(A1.get())
    plt.clf()
    inv_mass = ds["M"]
    newlowlim = np.power(10,low_lim)
    newuplim  = np.power(10,up_lim)
    limitedmasses = inv_mass[(inv_mass > newlowlim ) & (inv_mass < newuplim)]
    histogram = plt.hist(limitedmasses, bins=bins, range=(newlowlim,newuplim), color="cyan")
    y = histogram[0]
    x = 0.5*( histogram[1][0:-1] + histogram[1][1:])
    
    initials = [gamma, M, a, b, A]
    best, covariance = curve_fit(breitwigner, x, y, p0=initials, sigma=np.sqrt(y))
    error = np.sqrt(np.diag(covariance))
    plt.plot(x, breitwigner(x, *best), 'r-')
    #plt.yscale("log")
    plt.xlabel('Massa Invariante  [GeV]')
    plt.ylabel('Número de eventos ')
    plt.title('Ajuste de Breit-Wigner')
    plt.rcParams['figure.figsize'] = (8,6)
    plt.savefig("Hist.png")
    #plt.show()
    img = Image.open("Hist.png")
    tkimage = ImageTk.PhotoImage(img)
    Label(janela, image=tkimage).grid(row=1, column=6, rowspan=20)
        
    first = "Valor de gamma (FWMH) = {:4.4f} +- {:4.4f} \n".format(best[0], error[0])
    second = "Valor onde a distribuição M é máxima: = {:4.4f} +- {:4.4f} \n".format(best[1], error[1])
    third = "a = {:4.4f} +- {:4.4f} \n".format(best[2], error[2])
    fourth = "b = {:4.4f} +- {:4.4f} \n".format(best[3], error[3])
    fifth = "A = {:4.4f} +- {:4.4f} \n".format(best[4], error[4])
    message = str(first + second + third + fourth + fifth)
    messagebox.showinfo("ATENÇÃO!", message)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados1(): # Breit-Wigner
    Label(janela, text= "Entre com o valor da largura a meia altura (FWHM) do pico:").grid(row=7, column=0, columnspan=5)
    Entry(janela, textvar=gama1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor da posição do máximo da distribuição (M):").grid(row=8, column=0, columnspan=5)
    Entry(janela, textvar=M1).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor do parâmtro a que é a inclinação usada para perceber o efeito do background:").grid(row=9, column=0, columnspan=5)
    Entry(janela, textvar=a1).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor do parâmetro b de intercepção no eixo y usada para perceber o efeito do background:").grid(row=10, column=0, columnspan=5)
    Entry(janela, textvar=b1).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor do parâmetro A da altura da distribuição de Breit-Wigner:").grid(row=11, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=11, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular1).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados2(): # Duas Gaussianas + Exponencial
    Label(janela, text= "Entre com o valor ...:").grid(row=7, column=0, columnspan=5)
    Entry(janela, textvar=gama1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor ...:").grid(row=8, column=0, columnspan=5)
    Entry(janela, textvar=M1).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=9, column=0, columnspan=5)
    Entry(janela, textvar=a1).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=10, column=0, columnspan=5)
    Entry(janela, textvar=b1).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=11, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=11, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=12, column=0, columnspan=5)
    Entry(janela, textvar=a1).grid(row=12, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=13, column=0, columnspan=5)
    Entry(janela, textvar=b1).grid(row=13, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=14, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=14, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular1).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados3(): # Gaussiana + Exponencial
    Label(janela, text= "Entre com o valor ...:").grid(row=7, column=0, columnspan=5)
    Entry(janela, textvar=gama1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor ...:").grid(row=8, column=0, columnspan=5)
    Entry(janela, textvar=M1).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=9, column=0, columnspan=5)
    Entry(janela, textvar=a1).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=10, column=0, columnspan=5)
    Entry(janela, textvar=b1).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=11, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=11, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular1).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def Entrar_dados4(): # Crystalball + Exponencial
    Label(janela, text= "Entre com o valor ...:").grid(row=7, column=0, columnspan=5)
    Entry(janela, textvar=gama1).grid(row=7, column=5)
    
    Label(janela, text= "Entre com o valor ...:").grid(row=8, column=0, columnspan=5)
    Entry(janela, textvar=M1).grid(row=8, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=9, column=0, columnspan=5)
    Entry(janela, textvar=a1).grid(row=9, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=10, column=0, columnspan=5)
    Entry(janela, textvar=b1).grid(row=10, column=5)

    Label(janela, text= "Entre com o valor ...:").grid(row=11, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=11, column=5)
    
    Label(janela, text= "Entre com o valor ...:").grid(row=12, column=0, columnspan=5)
    Entry(janela, textvar=A1).grid(row=12, column=5)
    
    Button(janela, text="Calcular Fit", font=('Times', '12', 'bold'), width = 12, command=Calcular1).grid(row=0, column=1)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
Button(janela, text='Abrir arquivo', font=('Times', '12', 'bold'), width = 12, command=Abrir).grid(row=0, column=0)

#janela.geometry("800x600")
janela.mainloop() #mantém a janela aberta