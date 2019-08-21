import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.special import erf

##############################################################################
##################### DEFINIÇÕES DAS FUNÇÕES DOS AJUSTES #####################
##############################################################################
class funcoes:
    
    global expo
    def expo(x, const, slope):
        """ Uma curva exponencial. 
            parametros: const, slope.
        """
        return np.exp(const + slope*x)

    global line
    def line(x, intercept, slope):
        """ Polinômio do primeiro grau. 
        """
        return slope*x + intercept
    
    def breitwigner(E, gamma, M, a, b, A):
        ''' E (é a energia)
            gamma (a largura total do meio no máximo da distribuição)
            M (valor onde ocorre o máximo da distribuição)
            a (inclinação que é usada para pereber o efeito de backgrund)
            b (intercepção em y, que é usada para perceber o efeito de background)
            A ("amplitude" da distribuição de Breit-Wigner) '''
        return a*E+b+A*( (2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))) )/((E**2-M**2)**2+M**2*gamma**2)

    global gaussian
    def gaussian(x, a, x0, sigma):
        ''' a (altura do np.pico)
            x0 (média, ordena a posição do centro do np.pico)
            sigma (desvio padrão, controla a largura da curva)
        '''
        return a*np.exp(-(x-x0)**2/(2*sigma**2))

    global doublegaussian
    def doublegaussian(x, a1, x01, sigma1, a2, x02, sigma2):
        '''Duas gaussianas somadas.
        '''
        return a1*np.exp(-(x-x01)**2/(2*sigma1**2))+a2*np.exp(-(x-x02)**2/(2*sigma2**2))

    global crystalball
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
        D = np.sqrt(np.pi/2.) * (1. + erf(aa/np.sqrt(2.)))
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
    
    def histograma():
        # Leitura do arquivo que contém os dados.
        ds = pd.read_csv('DoubleMuRun2011A.csv')
        print(ds.head())

        # Fórmula da massa invariante, informação que será analisada no projeto.
        invariant_mass = np.sqrt(2*ds.pt1*ds.pt2*(np.cosh(ds.eta1-ds.eta2)-np.cos(ds.phi1-ds.phi2) ))
        print('Os primeiros cinco valores calculados (em unidades GeV)')
        print(invariant_mass[0:5])

        # Limitando o ajuste próximo ao pico do histograma.
        # Cada escolha seleciona uma ressonância distinta.
        global escolha 
        escolha = 0
        global esperado
        esperado = 0
        global err_esperado
        err_esperado = 0
        global esperado1 
        esperado1 = 0
        global err_esperado1
        err_esperado1 = 0
        while(escolha>4 or escolha<1):
            escolha = int(input("Escolha 1, 2, 3 ou 4: Enter 1 --> Z, Enter 2 --> Upsilon, Enter 3 --> J/Psi ou Enter 4 --> Psi':  "))
            if escolha == 1:
                lowerlimit = 70
                upperlimit = 110
                esperado = 91.1876
                err_esperado = 0.0021

            elif escolha == 2:
                lowerlimit = 9.2
                upperlimit = 10.2
                esperado = 9.46030
                err_esperado = 0.00026
                esperado1 = 10.02326
                err_esperado1 = 0.00031

            elif escolha == 3:
                lowerlimit = 2.95
                upperlimit = 3.2
                esperado = 3.096916
                err_esperado = 0.000011
        
            else:
                lowerlimit = 3.55
                upperlimit = 3.78
                esperado = 3.686093
                err_esperado = 0.00004
        
        # Definindo o número de bins que o histograma terá.
        global bins
        bins = int(input('Insira a binagem desejada: '))

        # Selecionando os valores de massa invariante que estão dentro das limitações.
        limitedmasses = invariant_mass[(invariant_mass > lowerlimit) & (invariant_mass < upperlimit)]

        # Criando um histograma com os valores selecionados.
        histogram = plt.hist(limitedmasses, bins=bins, range=(lowerlimit,upperlimit))

        # No eixo y, o número de eventos para cada bin (pode ser obtido a partir da variável histograma).
        # No eixo x, centro das classes.
        global y
        y = histogram[0]
        global x
        x = 0.5*( histogram[1][0:-1] + histogram[1][1:] )
        return escolha, esperado, err_esperado, histogram, x, y
    
    def initials():
        global initials_z
        global initials_upsilon 
        global initials_psiprime 
        global initials_jpsi 
        if escolha == 1:
            print(" gamma (a largura total do meio no maximo da distribuicao),\n M(valor onde ocorre o maximo da distribuicao),\n a (inclinacao que e usada para pereber o efeito de backgrund),\n b (intercepcao em y, que e usada para perceber o efeito de background),\n A (amplitude da distribuicao de Breit-Wigner). \nRecomendamos os valores: 4 91 -2 150 13000 para um fit bem sucedido.")
            initials_z =  [float(x) for x in input('Preencha com os parametros da Breit-Wigner (gamma, M, a, b, A) dando apenas um espaco entre eles: ').split()]
            return initials_z
        elif escolha  == 2: 
            print(" a (altura do pico), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1000 9.4 0.1 400 10 0.1 0 1 para um fit bem sucedido.")
            initials_upsilon =  [float(x) for x in input('Preencha com os parametros das duas Gaussianas e da exponencial (a1, mean1, sigma1, a2, mean2, sigma2, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
            return initials_upsilon
        elif escolha == 3: 
            print(" a (define como a funcao decresce no pico), \n n (), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1.6 0.9 3.1 1 0 -1 para um fit bem sucedido.")
            initials_jpsi =  [float(x) for x in input('Preencha com os parametros da Crystal-Ball e da exponencial (a, n, mean, sigma, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
            return initials_jpsi
        else:
            print(" a (define como a funcao decresce no pico), \n n (), \n mean (mean, ordena a posicao do centro do pico), \n sigma (desvio padrao, controla a largura da curva), \n const e inclinacao (np.exp(const + inclinacao*x)). \nRecomendamos os valores: 1 0 3.7 1 0 -1 para um fit bem sucedido.")
            initials_psiprime =  [float(x) for x in input('Preencha com os parametros da Crystal-Ball e da exponencial (a, n, mean, sigma, constante, inclinacao) dando apenas um espaco entre eles: ').split()]
            return initials_psiprime

    def fit(funcao, valor):
        best, covariance = curve_fit(funcao, x, y, p0=valor, sigma=np.sqrt(y))
        error = np.sqrt(np.diag(covariance))
        return best, covariance, error

    def print_valores(best, error):
        print("Valores com incertezas")
        print("")
        if escolha == 1:
            first = "Valor da largura de decaimento (gamma) = {} +- {}".format(best[0], error[0])
            second = "Valor do pico da distribuicao (M) = {} +- {}".format(best[1], error[1])
            third = "a = {} +- {}".format(best[2], error[2])
            fourth = "b = {} +- {}".format(best[3], error[3])
            fifth = "A = {} +- {}".format(best[4], error[4])
            print(first)
            print(second)
            print(third)
            print(fourth)
            print(fifth)
        elif escolha == 2:
            primeiro = "Valor de a1 = {} +- {}".format(best[0], error[0])
            segundo = "Valor de mean1 = {} +- {}".format(best[1], error[1])
            terceiro = "Valor de sigma1 = {} +- {}".format(best[2], error[2])
            quarto = "Valor de a2 = {} +- {}".format(best[3], error[3])
            quinto = "Valor de mean2 = {} +- {}".format(best[4], error[4])
            sexto = "Valor de sigma2 = {} +- {}".format(best[5], error[5])
            setimo = "Valor de constante = {} +- {}".format(best[6], error[6])
            oitavo = "Valor da inclinacao = {} +- {}".format(best[7], error[7])
            print(primeiro)
            print(segundo)
            print(terceiro)
            print(quarto)
            print(quinto)
            print(sexto)
            print(setimo)
            print(oitavo)
        elif escolha == 3:
            primeiro = "Valor de a = {} +- {}".format(best[0], error[0])
            segundo = "Valor de n = {} +- {}".format(best[1], error[1])
            terceiro = "Valor de mean = {} +- {}".format(best[2], error[2])
            quarto = "Valor de sigma = {} +- {}".format(best[3], error[3])
            quinto = "Valor da constante = {} +- {}".format(best[4], error[4])
            sexto = "Valor da inclinacao = {} +- {}".format(best[5], error[5])
            print(primeiro)
            print(segundo)
            print(terceiro)
            print(quarto)
            print(quinto)
            print(sexto)
        else:
            primeiro = "Valor de a = {} +- {}".format(best[0], error[0])
            segundo = "Valor de n = {} +- {}".format(best[1], error[1])
            terceiro = "Valor de mean = {} +- {}".format(best[2], error[2])
            quarto = "Valor de sigma = {} +- {}".format(best[3], error[3])
            quinto = "Valor da constante = {} +- {}".format(best[4], error[4])
            sexto = "Valor da inclinacao = {} +- {}".format(best[5], error[5])
            print(primeiro)
            print(segundo)
            print(terceiro)
            print(quarto)
            print(quinto)
            print(sexto)


    def convergencia(f, funcao, best, error, n):
        print("Numero de iteracoes: ", f)
        print("Baseado nos numeros de iteracoes: ")
        if (f == 0 and f == 1):
            print("O fit ficou bom? Legal! \nNao ficou? Tente outros valores iniciais! ")
        elif (f >= 1 and f <= 19):
            print("O fit convergiu!")
        else:
            print("O fit provavelmente esta divergindo... Tente outros valores iniciais!")

        if escolha == 2:     
            ch2, pval = chisquare(y, funcao(x, *best))
            test = ch2 /(bins - 8)
            discrepancia = np.absolute(best[1] - esperado)
            discrepancia1 = np.absolute(best[4] - esperado1)
            erro_padrao = np.sqrt(error[1]**2 + err_esperado**2)
            erro_padrao1 = np.sqrt(error[4]**2 + err_esperado1**2)

            print('A discrepancia do valor obtido no primeiro pico com o valor esperado e: %6.4f' %(discrepancia))
            print('A discrepancia do valor obtido no segundo pico com o valor esperado e: %6.4f' %(discrepancia1))
            print('A raiz da soma quadratica do erro do valor obtido no primeiro pico com o erro do valor esperado e: %6.4f' %(erro_padrao))
            print('A raiz da soma quadratica do erro do valor obtido no segundo pico com o erro do valor esperado e: %6.4f' %(erro_padrao1))
            print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
            print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pelas discrepancias: ")
            if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao and discrepancia < 2*erro_padrao1):
                print ("O fit e bom!")
            else:
                print("O fit e ruim! Tente outros valores.")

        else:
            ch2, pval = chisquare(y, funcao(x, *best))
            test = ch2 /(bins - n)
            discrepancia = np.absolute(best[1] - esperado)
            erro_padrao = np.sqrt(error[1]**2 + err_esperado**2)

            print('A discrepancia do valor obtido com o valor esperado e: %6.4f' %(discrepancia))
            print('A raiz da soma quadratica do erro do valor obtido com o erro do valor esperado e: %6.4f' %(erro_padrao))
            print('O valor do chi2 dividido pelo numero de graus de liberdade e: %6.3f' % (test))
            print("Baseado no valor de chi2 dividido pelo numero de graus de liberdade e pela discrepancia: ")

            if (test >=0 and test <= 10 and discrepancia < 2*erro_padrao):
                print ("O fit e bom!")
            else:
                print("O fit e ruim! Tente outros valores.")

    def plot(funcao, best, titulo):
        plt.plot(x, funcao(x, *best), 'r-')
        plt.xlabel('Massa Invariante [GeV]')
        plt.ylabel('Numero de Eventos')
        plt.title(titulo)
        plt.show()
