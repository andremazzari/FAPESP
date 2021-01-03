from numpy import sqrt,sin,linspace,cos,zeros,empty,tan,array,arange
import matplotlib.pyplot as plt
from time import time
from math import isnan
from Kerr_ODEs import *

t0 = time()

#parametros do problema
M = 1
a = 0.9
ro = 500*M
thetao = 1.4835
reh = M + sqrt(M**2 - a**2)
rehm = M - sqrt(M**2 - a**2)
rstatic = 2*M #limite estatico no plano equatorial
sigmao = ro**2 + (a*cos(thetao))**2
Deltao = ro**2 - 2*M*ro + a**2

def Kerr_solver(x0, u0, dt):
    t = 0
    x1 = x0[0]
    x2 = x0[1]
    x3 = x0[2]
    u1 = u0[0]
    u2 = u0[1]
    u3 = u0[2]
    
    r = []
    hemisfer = 0 #0:north 1:south
    EH = 0
    Dradius = 20*M
    while x1 <= x0[0]:
        x = array([t,x1,x2,x3], float)
        u = array([-1,u1,u2,u3], float)
        
        dx, du, EHRK = Runge_Kutta_45_step(x,u,dt)
        
        if EHRK == 1:
            EH = 1
            break

        x1 += dx[1]
        x2 += dx[2]
        x3 += dx[3]
        u1 += du[1]
        u2 += du[2]
        u3 += du[3]
        
        if hemisfer == 0 and cos(x[2]) < 0:#passou para hemisferio sul
            r.append(x1 - dx[1]/2) #faz a média das posições
            hemisfer = 1
        elif hemisfer == 1 and cos(x[2]) > 0:#passou para hemiferio norte
            r.append(x1 - dx[1]/2) #faz a média das posições
            hemisfer = 0
        
        if x1 <= (reh + 1e-2): #Horizonte de eventos
            EH = 1
            break
        
        if u1 > 0 and x1 > Dradius:
            break
        
        t += dt
        
    Lz, Q, H = Kerr_conserved_quantities([t,x1,x2,x3], [-1,u1,u2,u3])
    
    return r, EH, Q, H

def Emission(r): #perfil de emissao do disco de acreção
    if r < 20*M and r > 6*M:
        return 1
    else: 
        return 0

maxLz = 30*sin(thetao) #valor maximo de Lz

NumLz = 81 #Deve ser um numero impar
ListaLz = linspace(-maxLz, maxLz, NumLz)
delta = 2*maxLz/(NumLz - 1)
ObsInt = zeros([NumLz,NumLz], float)
Xcoord = empty([NumLz,NumLz],float)
Ycoord = empty([NumLz,NumLz], float)
g = 0

NanQ = []
NanL = []
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
        
        #Hemisferio sul
        R0 = ((ro**2 + a**2) - a*Lz)**2 - Deltao*((Lz - a)**2 + Q)
        Theta0 = Q - cos(thetao)**2*(-(a)**2 + (Lz/sin(thetao))**2)
        if abs(R0) < 1e-10:
            R0 = 0
        if abs(Theta0) < 1e-10 or Theta0 < 0:
            Theta0 = 0
        #momentos iniciais
        ur0 = -gamma_11([0,ro,thetao,0])*sqrt(R0)/sigmao
        utheta0 = gamma_22([0,ro,thetao,0])*sqrt(Theta0)/sigmao
        uphi0 = Lz
        
        H = Kerr_conserved_quantities([0,ro,thetao,0], [-1,ur0,utheta0,uphi0])[2]
        
        flag = 0
        dt = 5
        while flag == 0:
            r, EH, Qf, Hf  = Kerr_solver([ro, thetao, 0], [ur0, utheta0, uphi0], dt)
            if isnan(Hf) == True:
                if dt <= 0.2:
                    NanQ.append(Q)
                    NanL.append(Lz)
                    break
                else:
                    dt /= 5
            if abs(1 - Hf/H) < 0.001:
                flag = 1
            else:
                dt /= 5
            
        for i in arange(len(r)):
            ObsInt[int((NumLz-1)/2) - h][g] += Emission(r[i])
        if h != 0: #Hemisferio norte
            utheta0 = -gamma_22([0,ro,thetao,0])*sqrt(Theta0)/sigmao
            H = Kerr_conserved_quantities([0,ro,thetao,0], [-1,ur0,utheta0,uphi0])[2]
            flag = 0
            dt = 5
            while flag == 0:
                r, EH, Qf, Hf  = Kerr_solver([ro, thetao, 0], [ur0, utheta0, uphi0], dt)
                if isnan(Hf) == True:
                    if dt <= 0.2:
                        NanQ.append(-Q)
                        NanL.append(Lz)
                        break
                    else:
                        dt /= 5
                if abs(1 - Hf/H) < 0.001:
                    flag = 1
                else:
                    dt /= 5
            for i in arange(len(r)):
                ObsInt[int((NumLz-1)/2) + h][g] += Emission(r[i])
                

        h += 1
        if Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2 < 0:
            print(Q,Lz,h, Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
        Q += 2*delta*sqrt(Q + (a*cos(thetao))**2 - (Lz/tan(thetao))**2)
        print(h,g)
    g += 1
    
plt.pcolormesh(Xcoord, Ycoord, ObsInt)
plt.title("a = 0.9; \u03B8 = 85")
print("NanLz: ",NanL)
print("NanQ: ", NanQ)
print(len(NanL))
print(time() - t0)
