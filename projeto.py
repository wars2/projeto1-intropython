import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ds = pd.read_csv('DoubleMuRun2011A.csv') #Lê o arquivo com os dados.

#Parte de interação com o usuário
lowerlimit = int(input("Entre com o valor inteiro para o limite inferior do histograma:"))
upperlimit = int(input("Entre com o valor inteiro para o limite superior do histograma:"))
bins = int(input("Entre com um valor inteiro para o número de subdivisões do histograma:"))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
invariant_mass = np.sqrt(2*ds.pt1*ds.pt2*((np.cosh(ds.eta1-ds.eta2))-(np.cos(ds.phi1-ds.phi2)))) # Calcula o valor da massa invariante para um par de múons.
limitedmasses = invariant_mass[(invariant_mass > lowerlimit) & (invariant_mass < upperlimit)] # Seleciona as massas invariantes dentro dos limites escolhidos.

histogram = plt.hist(limitedmasses, bins=bins, range=(lowerlimit, upperlimit)) # Plota o histograma da massa invariante de acordo com os valores de entrada.
y = histogram[0] # Frequências referentes ao valor de bins usado
x = 0.5*( histogram[1][0:-1] + histogram[1][1:] ) # valor médio entre os limites de cada bin
plt.xlabel('Invariant mass [GeV]')
plt.ylabel('Number of event')
plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

gama1 = int(input("Entre com um valor inteiro para GAMA, ou seja, a largura a meia altura (FWHM) para o primeiro pico:"))
M1 =  int(input("Entre com um valor inteiro para M, ou seja, a posição da massa invariante onde se tem a máxima frequência para o primeiro pico:"))
print("")

gama2 = int(input("Entre com um valor inteiro para GAMA, ou seja, a largura a meia altura (FWHM) para o segundo pico:"))
M2 =  int(input("Entre com um valor inteiro para M, ou seja, a posição da massa invariante onde se tem a máxima frequência para o segundo pico:"))
print("")

gama3 = int(input("Entre com um valor inteiro para GAMA, ou seja, a largura a meia altura (FWHM) para o terceiro pico:"))
M3 =  int(input("Entre com um valor inteiro para M, ou seja, a posição da massa invariante onde se tem a máxima frequência para o terceiro pico:"))
print("")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def breitwigner(E, gamma, M, a, b, A):
    return a*E+b+A*((2*np.sqrt(2)*M*gamma*np.sqrt(M**2*(M**2+gamma**2)))/(np.pi*np.sqrt(M**2+np.sqrt(M**2*(M**2+gamma**2)))))/((E**2-M**2)**2+M**2*gamma**2)

# Valores iniciais para a otimização da curva.
# E  = Energia
# gamma Largura a meia altura (FWHM) da distribuição
# M posição onde se encontra o valor máximo da distribuição
# a inclinação usada para perceber o efeito do background
# b valor de intercepção no eixo y usada para perceber o efeito do background
# A a altura da distribuição de Breit-Wigner

initials1 = [gama1, M1, -2, 200, 13000]
initials2 = [gama2, M2, -2, 200, 13000]
initials3 = [gama3, M3, -2, 200, 13000]

# Let's import the module that is used in the optimization, run the optimization
# and calculate the uncertainties of the optimized parameters.
from scipy.optimize import curve_fit
from scipy.stats import exponnorm

best1, covariance1 = curve_fit(breitwigner, x, y, p0=initials1, sigma=np.sqrt(y))
error1 = np.sqrt(np.diag(covariance1))
    
best2, covariance2 = curve_fit(breitwigner, x, y, p0=initials2, sigma=np.sqrt(y))
error2 = np.sqrt(np.diag(covariance2))

best3, covariance3 = curve_fit(breitwigner, x, y, p0=initials3, sigma=np.sqrt(y))
error3 = np.sqrt(np.diag(covariance3))

# Let's print the values and uncertainties that are got from the optimization.
print("Valores e incertezas para otimização do fit:")
print("")
first = "Valor de gama = {} +- {}".format(best1[0], error1[0])
second = "Valor de M = {} +- {}".format(best1[1], error1[1])
third = "a = {} +- {}".format(best1[2], error1[2])
fourth = "b = {} +- {}".format(best1[3], error1[3])
fifth = "A = {} +- {}".format(best1[4], error1[4])
print(first)
print(second)
print(third)
print(fourth)
print(fifth)

plt.hist(limitedmasses, bins=bins, range=(lowerlimit, upperlimit))
plt.plot(x, breitwigner(x, *best1), 'r-', label='gamma1 = {}, M1 = {}'.format(best1[0], best1[1]))
plt.plot(x, breitwigner(x, *best2), 'k-', label='gamma2 = {}, M2 = {}'.format(best2[0], best2[1]))
plt.plot(x, breitwigner(x, *best3), 'g-', label='gamma3 = {}, M3 = {}'.format(best3[0], best3[1]))
plt.xlabel('Invariant mass [GeV]')
plt.ylabel('Number of event')
plt.title('Histogram and Breit-Wigner fit')
plt.legend()
plt.show()