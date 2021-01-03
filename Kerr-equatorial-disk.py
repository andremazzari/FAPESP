from numpy import empty,zeros,arcsin,tan,amin,linspace,cbrt,sqrt,cos,arctan2,pi,arccos,arctan,sin,array,append
from scipy.special import ellipj, ellipkinc, ellipk
import matplotlib.pyplot as plt
from matplotlib import colors
from time import time

t0 = time()

#parametros do problema
M = 1
a = 0.1
ro = 10000*M
thetao = 1.4835
reh = M + sqrt(M**2 - a**2)
rehm = M - sqrt(M**2 - a**2)
rstatic = 2*M #limite estatico no plano equatorial
DM = 1.15e-3

Z1 = 1 + (1 - a**2)**(1/3)*((1 + a)**(1/3) + (1 - a)**(1/3))
Z2 = sqrt(3*a**2 + Z1**2)
risco = 3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2*Z2))

print(risco)

def R(s, ro, Rcase, roots, rij):
    if Rcase == 1:
        k = ((rij[2][1]*rij[3][0])/(rij[2][0]*rij[3][1])).real
        X = sqrt((rij[2][0]*rij[3][1]).real)*s/2
        X -= ellipkinc(arcsin(sqrt((((ro - roots[3])*rij[2][0])/((ro - roots[2])*rij[3][0])).real)), k)
        r = rij[2][0]*roots[3] - rij[3][0]*roots[2]*ellipj(X, k)[0]**2
        r /= rij[2][0] - rij[3][0]*ellipj(X, k)[0]**2
    elif Rcase == 2:
        A = sqrt((rij[2][1]*rij[3][1]).real)
        B = sqrt((rij[2][0]*rij[3][0]).real)
        
        x3o = ((A*(ro - roots[0]) - B*(ro - roots[1]))/(A*(ro - roots[0]) + B*(ro - roots[1]))).real
        k3 = (((A + B)**2 - (roots[1] - roots[0])**2)/(4*A*B)).real
        X3 = sqrt(A*B)*s - ellipkinc(arccos(x3o), k3)
        r = B*roots[1] - A*roots[0] + (A*roots[0] + B*roots[1])*ellipj(X3, k3)[1]
        r /= B - A + (A + B)*ellipj(X3, k3)[1]
    else:
        C = sqrt((rij[2][0]*rij[3][1]).real)
        D = sqrt((rij[2][1]*rij[3][0]).real)
        
        x4o = (ro - roots[1].real)/roots[1].imag
        g0 = sqrt((4*roots[1].imag**2 - (C - D)**2)/((C + D)**2 - 4*roots[1].imag**2))
        k4 = 4*C*D/(C + D)**2
        X4 = ellipkinc(arctan(x4o) + arctan(g0), k4) - (C + D)*s/2
        sc = ellipj(X4, k4)[0]/ellipj(X4, k4)[1]
        r = roots[1].imag*((sc - g0)/(1 + g0*sc)) + roots[1].real
        
    return r.real
#Funções para o modelo do disco de acreção
def FuncAlpha(a,r):
    return 1 + (a/r)**2 + 2*a**2/r**3

def FuncBeta(a,r):
    return 1 + a/r**(3/2)
    
def FuncL(a,r):
    return 1 - 3/r + 2*a/r**(3/2)
    
def FuncD(a,r):
    return 1 - 2/r + (a/r)**2
    
def FuncE(a,r):
    return 1 + 4*(a/r)**2 - 4*a**2/r**3 + 3*(a/r)**4
    
def Temperature(a,r): #Kelvin 
    #rom = (2e3)*(M**(-2/3)*DM**(2/3))*FuncAlpha(a,r)**(2/3)*FuncBeta(a,r)**(-8/15)*FuncD(a,r)**(-1/3)*FuncE(a,r)**(-1/3)*FuncL(a,r)**(2/3)
    #rmi = (40)*(M**(-2/3)*DM**(16/20))*FuncAlpha(a,r)**(20/21)*FuncBeta(a,r)**(-36/21)*FuncD(a,r)**(-8/21)*FuncE(a,r)**(-10/21)*FuncL(a,r)**(16/21)
    rom = 100
    rmi = 4
    if r > rom:#outer region
        return (8e7)*(M**(-1/2)*DM**(3/10))*r**(-3/4)*FuncAlpha(a,r)**(-1/10)*FuncBeta(a,r)**(-1/5)*FuncD(a,r)**(-3/20)*FuncE(a,r)**(1/20)*FuncL(a,r)**(3/10)
    elif r <= rom and r > rmi: #Middle region
        return (3e8)*(M**(-3/5)*DM**(2/5))*r**(-9/10)*FuncBeta(a,r)**(-2/5)*FuncD(a,r)**(-1/5)*FuncL(a,r)**(2/5)
    else:#inner region
        return (4e7)*(M**(-1/4))*r**(-3/8)*FuncAlpha(a,r)**(-1/2)*FuncBeta(a,r)**(1/2)*FuncE(a,r)**(1/20)
    
def Density(a,r):
    #rom = (2e3)*(M**(-2/3)*DM**(2/3))*FuncAlpha(a,r)**(2/3)*FuncBeta(a,r)**(-8/15)*FuncD(a,r)**(-1/3)*FuncE(a,r)**(-1/3)*FuncL(a,r)**(2/3)
    #rmi = (40)*(M**(-2/3)*DM**(16/20))*FuncAlpha(a,r)**(20/21)*FuncBeta(a,r)**(-36/21)*FuncD(a,r)**(-8/21)*FuncE(a,r)**(-10/21)*FuncL(a,r)**(16/21)
    rom = 100
    rmi = 4
    if r > rom:#outer region
        return (8e1)*(M**(-5/4)*DM**(11/20))*r**(-15/8)*FuncAlpha(a,r)**(-17/20)*FuncBeta(a,r)**(3/10)*FuncD(a,r)**(-11/40)*FuncE(a,r)**(17/40)*FuncL(a,r)**(11/20)
    elif r <= rom and r > rmi: #Middle region
        return (10)*(M**(-11/10)*DM**(2/5))*r**(-33/20)*FuncAlpha(a,r)**(-1)*FuncBeta(a,r)**(3/5)*FuncD(a,r)**(-1/5)*FuncE(a,r)**(1/2)*FuncL(a,r)**(2/5)
    else:#inner region
        return (1e-4)*(M**(-1/4)*DM**(-2))*r**(3/2)*FuncAlpha(a,r)**(-4)*FuncBeta(a,r)**(6)*FuncD(a,r)*FuncE(a,r)**(2)*FuncL(a,r)**(-2)

def MagneticField(a,r):
    rom = 100
    rmi = 4
    if r > rom:#outer region
        return (7e8)*(M**(-7/8)*DM**(17/40))*r**(-21/16)*FuncAlpha(a,r)**(-19/40)*FuncBeta(a,r)**(1/20)*FuncD(a,r)**(-17/80)*FuncE(a,r)**(19/80)*FuncL(a,r)**(17/40)
    elif r <= rom and r > rmi: #Middle region
        return (1e9)*(M**(-17/20)*DM**(2/5))*r**(-51/40)*FuncAlpha(a,r)**(-1/2)*FuncBeta(a,r)**(1/10)*FuncD(a,r)**(-1/5)*FuncE(a,r)**(1/4)*FuncL(a,r)**(2/5)
    else:#inner region
        return (7e7)*(M**(-1/2))*r**(-3/4)*FuncAlpha(a,r)**(-1)*FuncBeta(a,r)*FuncE(a,r)**(1/2)

def FreeBoundEmission(a,r):
    Const = 1e25
    if r < risco:
        return 0
    else:
        return Const*Density(a,r)*Temperature(a,r)**(-1/2)
    
def FreeFreeEmission(a,r):
    return FreeBoundEmission(a,r)*Temperature(a,r)/(8e5)

def CyclotronEmission(a,r):
    Const = 1
    if r < risco:
        return 0
    else:
        if Temperature(a,r) > 1e10:#Syncrotron
            return Const*(1e-9)*(Temperature(a,r)*MagneticField(a,r))**2
        else: #Cyclotron'''
            return Const*Temperature(a,r)*MagneticField(a,r)**2

def Emission(r): #perfil de emissao do disco de acreção
    '''
    if r < 20*M and r > 6*M:
        E = 1
    else: 
        E = 0
    '''
    E = (FreeBoundEmission(a,r) + FreeFreeEmission(a,r) + CyclotronEmission(a,r))
    
    if E > 1e15:
        E /= 1e17
    

        
    return E

def Energy(r, theta, Lz):
    flag = 1
    Energy = Emission(r)
    if Energy == 0:
        return 0
    else:
        if flag == 0:
            rho = r**2 + (a*cos(theta))**2
            
            gtt = -(1 - 2*M*r/rho)
            gpp = (r**2 + a**2 + 2*M*r*a**2/rho)*sin(theta)**2
            gtp = - 2*M*r*a*sin(theta)**2/rho
            
            omega = sqrt(M/r**3)
            
            ut = 1/sqrt(-gtt - 2*omega*gtp - omega**2*gpp)
            
            return Energy/(ut*(1 - omega*Lz))**4
        else:
            f = sqrt(r**3 - 3*M*r**2  + 2*a*r**(3/2))/(r**(3/2) + a - Lz)
            return Energy*f**4

maxLz = 50*sin(thetao)

NumLz = 151 #Deve ser um numero impar
ListaLz = linspace(-maxLz, maxLz, NumLz)
delta = 2*maxLz/(NumLz - 1)
ObsInt = zeros([NumLz,NumLz], float)
Xcoord = empty([NumLz,NumLz],float)
Ycoord = empty([NumLz,NumLz], float)
g = 0

for Lz in ListaLz:
    h = 0
    Q = -(a*cos(thetao))**2 + (Lz/tan(thetao))**2 + 0.00001
    
    while h <= (NumLz-1)/2:
        #coordenadas dos pontos calculados
        if h == 0:#y=0
            Xcoord[int((NumLz-1)/2)][g] = -Lz/sin(thetao)
            Ycoord[int((NumLz-1)/2)][g] = 0
        else:
            Xcoord[int((NumLz-1)/2) + h][g] = -Lz/sin(thetao)
            Ycoord[int((NumLz-1)/2) + h][g] = sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
            Xcoord[int((NumLz-1)/2) - h][g] = -Lz/sin(thetao)
            Ycoord[int((NumLz-1)/2) - h][g] = -sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
            
        if Q <= 0: #trajetoria não cruza o plano equatorial
            h += 1
            Q += 2*delta*sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
            continue
        
        #parametros do potencial angular
        deltaT = 0.5*(1 - (Q + Lz**2)/a**2)
        up = deltaT + sqrt(deltaT**2 + Q/a**2)
        um = deltaT - sqrt(deltaT**2 + Q/a**2)
        
        if up<=0 or um == 0:
            print("erro 2:", up,um,Q, Lz)
        if up>1:
            if up < 1.00000001:
                up = 1.0
            else:
                print("erro 3:", up,um,Q,Lz)
                
        theta1 = arccos(sqrt(up))
        theta4 = arccos(-sqrt(up))
        if thetao < theta1 or thetao > theta4:
            h += 1
            Q += 2*delta*sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
            print("aviso 1", Q,Lz)
            continue
        if abs(cos(thetao)/sqrt(up)) > 1:
            print("erro 4:",cos(thetao)/sqrt(up), thetao, theta1, theta4)
        
        G = (2/sqrt(-um*a**2))*ellipk(up/um)
        Go = -(1/sqrt(-um*a**2))*ellipkinc(arcsin(cos(thetao)/sqrt(up)), up/um)
        
        #Parametros para o calculo das raizes do potencial radial
        A = a**2 - Q - Lz**2
        B = 2*M*(Q + (Lz - a)**2)
        C = -a**2*Q
        
        P = -A**2/12 - C
        S = -(A/3)*((A/6)**2 - C) - B**2/8
        
        if ((P/3)**3 + (S/2)**2) >= 0:
            omegap = cbrt(-S/2 + sqrt((P/3)**3 + (S/2)**2))
            omegam = cbrt(-S/2 - sqrt((P/3)**3 + (S/2)**2))
            e0 = omegap + omegam - A/3
        else:
            x = -S/2
            y = sqrt(-((P/3)**3 + (S/2)**2))
            rp = sqrt(x**2 + y**2)
            thetap = arctan2(y,x)
            imax = 0
            for i in [1,2]:
                if cos(thetap/3 + (2/3)*pi*i) >= cos(thetap/3 +(2/3)*pi*imax):
                    imax = i
            e0 = 2*cbrt(rp)*cos(thetap/3 + (2/3)*pi*imax) - A/3
        
        if e0 <= 0:
            print("erro 1", e0)
            
        z = sqrt(e0/2) 
        #calcula as raizes do potencial radial em cada caso
        Roots = empty([4], complex)
        if (-A/2 - z**2 - B/(4*z)) >= 0: #quatro raizes reais
            Rcase = 1
        
            Roots[0] = -z - sqrt(-A/2 - z**2 + B/(4*z))
            Roots[1] = -z + sqrt(-A/2 - z**2 + B/(4*z))
            Roots[2] = z - sqrt(-A/2 - z**2 - B/(4*z))
            Roots[3] = z + sqrt(-A/2 - z**2 - B/(4*z))
            
        elif (-A/2 - z**2 + B/(4*z)) >= 0: #duas raizes reais
            Rcase = 2
        
            Roots[0] = -z - sqrt(-A/2 - z**2 + B/(4*z))
            Roots[1] = -z + sqrt(-A/2 - z**2 + B/(4*z))
            Roots[2] = complex(z , - sqrt(-(-A/2 - z**2 - B/(4*z))))
            Roots[3] = complex(z , sqrt(-(-A/2 - z**2 - B/(4*z))))
            
        else: #nenhuma raiz real
            Rcase = 3
        
            Roots[0] = complex(-z , - sqrt(-(-A/2 - z**2 + B/(4*z))))
            Roots[1] = complex(-z , sqrt(-(-A/2 - z**2 + B/(4*z))))
            Roots[2] = complex(z , - sqrt(-(-A/2 - z**2 - B/(4*z))))
            Roots[3] = complex(z , sqrt(-(-A/2 - z**2 - B/(4*z))))
        
        #calcula as diferenças entre as raízes
        rij = empty([4,4], complex)
        for i in [0,1,2,3]:
            for j in [0,1,2,3]:
                rij[i][j] = Roots[i] - Roots[j]
                
        #calcula o parametro máximo em cada caso
        smr = 0 #parametro em que teriamos um turning point radial 
        if Rcase == 1:
            k = ((rij[2][1]*rij[3][0])/(rij[2][0]*rij[3][1])).real
            smax = ellipkinc(arcsin(sqrt((((ro - Roots[3])*rij[2][0])/((ro - Roots[2])*rij[3][0])).real)), k)
            if Roots[3].real < reh:
                smax -= ellipkinc(arcsin(sqrt((((reh - Roots[3])*rij[2][0])/((reh - Roots[2])*rij[3][0]))).real), k)
            else:
                smr = 2/sqrt((rij[2][0]*rij[3][1]).real)*smax
                smax += ellipkinc(arcsin(sqrt((((ro - Roots[3])*rij[2][0])/((ro - Roots[2])*rij[3][0]))).real), k)
            smax *= 2/sqrt((rij[2][0]*rij[3][1]).real)
        elif Rcase == 2:
            A = sqrt((rij[2][1]*rij[3][1]).real)
            B = sqrt((rij[2][0]*rij[3][0]).real)
            
            x3o = ((A*(ro - Roots[0]) - B*(ro - Roots[1]))/(A*(ro - Roots[0]) + B*(ro - Roots[1]))).real
            x3eh = ((A*(reh - Roots[0]) - B*(reh - Roots[1]))/(A*(reh - Roots[0]) + B*(reh - Roots[1]))).real
            k3 = (((A + B)**2 - (Roots[1] - Roots[0])**2)/(4*A*B)).real
            
            smax = ellipkinc(arccos(x3o), k3)
            smax -= ellipkinc(arccos(x3eh), k3)
            smax /= sqrt(A*B)
        else:
            C = sqrt((rij[2][0]*rij[3][1]).real)
            D = sqrt((rij[2][1]*rij[3][0]).real)
            
            x4o = (ro - Roots[1].real)/Roots[1].imag
            x4eh = (reh - Roots[1].real)/Roots[1].imag
            g0 = sqrt((4*Roots[1].imag**2 - (C - D)**2)/((C + D)**2 - 4*Roots[1].imag**2))
            k4 = 4*C*D/(C + D)**2
            
            smax = ellipkinc(arctan(x4o) + arctan(g0), k4)
            smax -= ellipkinc(arctan(x4eh) + arctan(g0), k4)
            smax *= 2/(C + D)
        
        #Hemisferio sul
        p = 0
        s = p*G - Go
        while s < smax:
            r = R(s, ro, Rcase, Roots, rij)
            ObsInt[int((NumLz-1)/2) - h][g] += Energy(r,pi/2,Lz)
            p += 1
            s = p*G - Go
        if h != 0: #Hemisferio norte
            p = 1
            s = p*G + Go
            while s < smax:
                r = R(s, ro, Rcase, Roots, rij)
                ObsInt[int((NumLz-1)/2) + h][g] += Energy(r,pi/2,Lz)
                p += 1
                s = p*G + Go
                

        h += 1
        if Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2 < 0:
            print(Q,Lz,h, Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
        Q += 2*delta*sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
    g += 1

#plt.pcolormesh(Xcoord, Ycoord, ObsInt)
plt.pcolormesh(Xcoord, Ycoord, ObsInt, norm=colors.SymLogNorm(linthresh = 1))
plt.title("a = " + str(a) + " ; \u03B8 = 85")
print(time() - t0)
plt.figure()
x = linspace(reh,250*reh,1000)
E = [Emission(r) for r in x]
plt.yscale('log')
plt.title("Perfil de emissão do disco de acreção")
plt.xlabel("r/M")
plt.ylabel("Intensidade Emitida")
plt.plot(x,E)


        
        