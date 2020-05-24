#Com atenuação, NAO esta adptado para emissao nao-uniforme
#DiskV1 representa uma primeira versão do disco de acreção, mas que estava apresentando algum problema que ainda nao entendi
#DiskV2 parece funcionar
from numpy import sqrt,arctan,cos,sin,linspace,zeros,arange,searchsorted,array,pi,ones,copy,tan,empty,arctan2,meshgrid,tile,tanh,flip,sign,arccos,concatenate
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time

t0 = time()

M = 1 #Massa do buraco negro
d = 100*M #Distância do plano à origem
E = 1 #Energia do raio de luz
N = 200 #numero de divisoes na integração
Nd = 501 #Numero de divisoes no grid para plotar o disco de acreção
Np = 300 #Numero de parametros de impacto utilizados
Nf = 100 #Numero de divisoes no ajusto fino do intervalo de integração
NumAngulos = 70 #Numero de angulos que serão utilizados
alpha = 1.39626 #Angulo de inclinação do disco


#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y=d/ds(u) , Y[2]=phi, Y[3]=z=d/ds(phi)
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*E*Y[0]**3, Y[3], 2*Y[0]*Y[1]*b*E]

#define evento de atravessar o horizonte de eventos
def EventHorizon(s, Y):
    return 1/Y[0] - 2*M
EventHorizon.terminal = True #define este evento como terminal

#versao 2 do perfil de acreção
def Emission(r, phi, theta):
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)
    
    ypa = y*sin(alpha)**2 + z*sin(alpha)*cos(alpha)
    zpa = y*sin(alpha)*cos(alpha) + z*cos(alpha)**2
    rpa = sqrt(ypa**2 + zpa**2)
    
    ype = y - ypa
    zpe = z - zpa
    rpe = sqrt(x**2 + ype**2 + zpe**2)
    
    SinPhi = rpe/r
    CosPhi = rpa/r
    
    if abs(r*SinPhi) > 2*(r*CosPhi)**2 + 2:
        return 1
    else:
        return 0
    
def Absortion(r, phi, theta):
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)
    
    ypa = y*sin(alpha)**2 + z*sin(alpha)*cos(alpha)
    zpa = y*sin(alpha)*cos(alpha) + z*cos(alpha)**2
    rpa = sqrt(ypa**2 + zpa**2)
    
    ype = y - ypa
    zpe = z - zpa
    rpe = sqrt(x**2 + ype**2 + zpe**2)
    
    SinPhi = rpe/r
    CosPhi = rpa/r
    
    Ta = 0.9 #Taxa de absorção
    
    if abs(r*SinPhi) > 2*(r*CosPhi)**2 + 2:
        return Ta
    else:
        return 0

#Evento em que o raio entra e sai do circulo com r=10sqrt(2)M
def LimitInterval(s, Y):
    return 1/Y[0] - 10*sqrt(2)*M

def IntegrationInterval(i, state, theta):
    x = linspace(sol.t[i-1], sol.t[i], Nf)
    if state==0: #Entrando no disco
        for s in x:
            if Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta)!=0:
                return s
    else: #Saindo do disco
        sp = x[0]
        for s in x:
            if Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta)==0:
                return sp
            sp = s

#Lista de parametros de impacto que serao usados
ListaParametros = linspace(0, 10*M, Np)
#lista de angulos que serão utilizados
Angulos = linspace(0,pi/2,NumAngulos)

#Array onde será salvo a intensidade observada em cada ponto de observação
ObsInt = zeros((4*NumAngulos-3,len(ListaParametros)))
MaxObsInt = 0 #Intensidade maxima observada

k = 0 #indice radial(parametro de impacto)
for b in ListaParametros:
    #valores iniciais
    r0 = sqrt(b**2 + d**2)
    u0 = 1/r0
    phi0 = arctan(b/d)
    y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
    z0 = u0**2*b*E
    
    sol = solve_ivp(F, [0, 150], [u0, y0, phi0, z0], events=(EventHorizon, LimitInterval), dense_output=True, max_step=0.1)
    
    #Encontra os limites da regiao de interesse da trajetoria do raio
    if len(sol.t_events[0]) != 0:
        LimitPos = array([searchsorted(sol.t, sol.t_events[1][0]), len(sol.t) - 1])
    else:
        LimitPos = searchsorted(sol.t, [sol.t_events[1][0], sol.t_events[1][1]])
        
    #realiza a integração para cada um dos angulos
    m = 0 #indice angular
    for theta in Angulos:
        
        if Emission(1/sol.sol(sol.t[int(LimitPos[0])])[0], sol.sol(sol.t[int(LimitPos[0])])[2], theta)!=0:
            state = 1
            si = [sol.t[int(LimitPos[0])]] #TALVEZ PRECISE DE AJUSTE FINO AQUI
            sf = []
        else:
            state = 0
            si = []
            sf = []
        
        #define os intervalos de integração
        for i in arange(LimitPos[0], LimitPos[1]):
            if Emission(1/sol.sol(sol.t[i])[0], sol.sol(sol.t[i])[2], theta) != state:
                if state==0:
                    si.append(IntegrationInterval(i,state,theta))
                    state = 1
                else:
                    sf.append(IntegrationInterval(i,state,theta))
                    state = 0
                    
        #Caso ainda nao tenha nenhum intervalo de integração, verifica se perdeu algo no final
        if len(si)==0:
            if len(sol.t_events[0])!=0:
                sp = sol.t_events[0][0]
            else:
                sp = sol.t_events[1][1]
            x = linspace(sol.t[LimitPos[1]-1], sp, Nf)
            for s in x:
                if Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta)!=0:
                    si.append(s)
                    sf.append(sp)
                    break
        
        if len(si)!=len(sf):
            if len(sol.t_events[0])!=0: #Entrou no horizonte de eventos
                if Emission(1/sol.sol(sol.t_events[0][0])[0], sol.sol(sol.t_events[0][0])[2], theta)!=0:
                    sf.append(sol.t_events[0][0] - 0.1)
                else:
                    x = linspace(sol.t[LimitPos[1]-1], sol.t_events[0][0], Nf)
                    sp = x[0]
                    for s in x:
                        if Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta)==0:
                            sf.append(sp)
                            break
                        sp = s
            elif Emission(1/sol.sol(sol.t_events[1][1])[0], sol.sol(sol.t_events[1][1])[2], theta)!=0:
                sf.append(sol.t_events[1][1])
            else:
                sf.append(IntegrationInterval(LimitPos[1],state,theta))
        if len(si)!=len(sf):
            print("erro2", k, b)
            print(si)
            print(sf)
            k += 1
            continue
        
        rp = 1*M
        #realiza a integração
        for i in flip(arange(len(si))):
            h = (sf[i] - si[i])/N
            
            s = sf[i] #começa pelo fim
            
            for j in range(N):
                if (j!=0) or (i != (len(si) - 1)): #Não é o primeiro ponto
                    #fator de red-shift
                    g = (1 - 2*M/rp)/(1 - 2*M*sol.sol(s)[0])
                    
                    Delta = Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta) - Absortion(1/sol.sol(s)[0], sol.sol(s)[2], theta)*pow(g, 3/2)*ObsInt[m][k]
                
                    #atualiza o red-shift
                    ObsInt[m][k] *= g**2
    
                    ObsInt[m][k] += h*pow(g,1/2)*Delta*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
                else:#caso seja o primeiro ponto
                    ObsInt[m][k] += h*Emission(1/sol.sol(s)[0], sol.sol(s)[2], theta)*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
                
                rp = 1/sol.sol(s)[0]
                s -= h
                
        #realiza o red-shift para o infinito
        ObsInt[m][k] *= (1 - 2*M/rp)**2
        
        #segundo quadrante é igual
        if theta!=0 and theta!=pi/2:
            ObsInt[2*NumAngulos-m-2][k] = ObsInt[m][k]
            
        #realiza o mesmo para os terceiro e quarto quadrantes
        
        l = 2*NumAngulos - 2 + m
        
        if Emission(1/sol.sol(sol.t[int(LimitPos[0])])[0], sol.sol(sol.t[int(LimitPos[0])])[2], -theta)!=0:
            state = 1
            si = [sol.t[int(LimitPos[0])]] #TALVEZ PRECISE DE AJUSTE FINO AQUI
            sf = []
        else:
            state = 0
            si = []
            sf = []
        
        #define os intervalos de integração
        for i in arange(LimitPos[0], LimitPos[1]):
            if Emission(1/sol.sol(sol.t[i])[0], sol.sol(sol.t[i])[2], -theta) != state:
                if state==0:
                    si.append(IntegrationInterval(i,state,-theta))
                    state = 1
                else:
                    sf.append(IntegrationInterval(i,state,-theta))
                    state = 0
                    
        #Caso ainda nao tenha nenhum intervalo de integração, verifica se perdeu algo no final
        if len(si)==0:
            if len(sol.t_events[0])!=0:
                sp = sol.t_events[0][0]
            else:
                sp = sol.t_events[1][1]
            x = linspace(sol.t[LimitPos[1]-1], sp, Nf)
            for s in x:
                if Emission(1/sol.sol(s)[0], sol.sol(s)[2], -theta)!=0:
                    si.append(s)
                    sf.append(sp)
                    break
        
        if len(si)!=len(sf):
            if len(sol.t_events[0])!=0: #Entrou no horizonte de eventos
                if Emission(1/sol.sol(sol.t_events[0][0])[0], sol.sol(sol.t_events[0][0])[2], -theta)!=0:
                    sf.append(sol.t_events[0][0] - 0.1)
                else:
                    x = linspace(sol.t[LimitPos[1]-1], sol.t_events[0][0], Nf)
                    sp = x[0]
                    for s in x:
                        if Emission(1/sol.sol(s)[0], sol.sol(s)[2], -theta)==0:
                            sf.append(sp)
                            break
                        sp = s
            elif Emission(1/sol.sol(sol.t_events[1][1])[0], sol.sol(sol.t_events[1][1])[2], -theta)!=0:
                sf.append(sol.t_events[1][1])
            else:
                sf.append(IntegrationInterval(LimitPos[1],state, -theta))
        if len(si)!=len(sf):
            print("erro2", k, b)
            print(si)
            print(sf)
            k += 1
            continue
        
        rp = 1*M
        #realiza a integração
        for i in flip(arange(len(si))):
            h = (sf[i] - si[i])/N
            
            s = sf[i] #começa pelo fim
            
            for j in range(N):
                if (j!=0) or (i != (len(si) - 1)): #Não é o primeiro ponto
                    #fator de red-shift
                    g = (1 - 2*M/rp)/(1 - 2*M*sol.sol(s)[0])
                    
                    Delta = Emission(1/sol.sol(s)[0], sol.sol(s)[2], -theta) - Absortion(1/sol.sol(s)[0], sol.sol(s)[2], -theta)*pow(g, 3/2)*ObsInt[l][k]
                
                    #atualiza o red-shift
                    ObsInt[l][k] *= g**2
    
                    ObsInt[l][k] += h*pow(g,1/2)*Delta*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
                else:#caso seja o primeiro ponto
                    ObsInt[l][k] += h*Emission(1/sol.sol(s)[0], sol.sol(s)[2], -theta)*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
                
                rp = 1/sol.sol(s)[0]
                s -= h
                
        #realiza o red-shift para o infinito
        ObsInt[l][k] *= (1 - 2*M/rp)**2
                
        #quarto quadrante é igual
        if theta!=pi/2:
            pos = l % (2*NumAngulos-2)
            ObsInt[4*NumAngulos-4-pos][k] = ObsInt[l][k]
        
            
        m += 1
        
    print(k)
    k += 1

#plota o perfil do disco de acreção
AccretionDisk = empty([Nd,Nd], int)
for i in range(Nd):
    for j in range(Nd):
        x = (20/Nd)*i - 10
        y = (20/Nd)*j - 10
        AccretionDisk[j][i] = Emission(sqrt(x**2 + y**2), arctan2(y,x), pi/2)
fig, ax = plt.subplots()
plt.pcolormesh(linspace(-10,10,Nd), linspace(-10,10,Nd), AccretionDisk)
circle1 = plt.Circle((0, 0), 2*M, color='black')
ax.add_patch(circle1)
plt.axis("square")
plt.show()

Angulos= concatenate((linspace(0,pi/2,NumAngulos),linspace(pi/2,pi,NumAngulos)[1:NumAngulos-1],linspace(pi,3/2*pi,NumAngulos),linspace(3/2*pi,2*pi,NumAngulos)[1:NumAngulos]))

#Agora plota a imagem do buraco negro
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

r, th = meshgrid(ListaParametros, Angulos)

plt.pcolormesh(th, r, ObsInt)
plt.colorbar(label='Observed Intensity')

print(time() - t0)