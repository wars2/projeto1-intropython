import numpy as np

##############################################################################
##################### DEFINIÇÕES DAS FUNÇÕES DOS AJUSTES #####################
##############################################################################
class funcoes:

    def expo(x, const, slope):
        """ Uma curva exponencial. 
            parametros: const, slope.
        """
        return np.exp(const + slope*x)

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

    
