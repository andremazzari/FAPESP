from numpy import sqrt,arctan,cos,sin,linspace,zeros,arange,searchsorted,array,pi,ones,copy,tan,empty,arctan2,meshgrid,tile,tanh,flip,log,exp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from time import time

t0 = time()

M = 1 #Massa do buraco negro
d = 100*M #Distância do plano à origem
E = 1 #Energia do raio de luz
N = 1000 #Numero de divisoes na integração trapeoidal
Nd = 501 #Numero de divisoes no grid para plotar o disco de acreção
Np = 1000 #Numero de parametros de impacto utilizados
Nf = 100 #Numero de divisoes no ajusto fino do intervalo de integração


#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y=d/ds(u) , Y[2]=phi, Y[3]=z=d/ds(phi)
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*E*Y[0]**3, Y[3], 2*Y[0]*Y[1]*b*E]

#define evento de atravessar o horizonte de eventos
def EventHorizon(s, Y):
    return 1/Y[0] - 2*M
EventHorizon.terminal = True #define este evento como terminal

#perfil de emissão
def Emission(r, phi):
    if r < 6*M:
        return 1
    else:
        return 0

#perfil de absroção
def Absortion(r, phi):
    a = log(9)/(4*M)
    Ta = 0.9*exp(-a*(r - 2*M)) #Taxa de absorção
    if r < 6*M:
        return 1
    else:
        return 0

#Evento em que o raio entra e sai do circulo com r=10sqrt(2)M
def LimitInterval(s, Y):
    return 1/Y[0] - 10*sqrt(2)*M

def IntegrationInterval(i, state):
    x = linspace(sol.t[i-1], sol.t[i], Nf)
    if state==0: #Entrando no disco
        for s in x:
            if Emission(1/sol.sol(s)[0], sol.sol(s)[2])!=0:
                return s
    else: #Saindo do disco
        sp = x[0]
        for s in x:
            if Emission(1/sol.sol(s)[0], sol.sol(s)[2])==0:
                return sp
            sp = s

#Lista de parametros de impacto que serao usados
ListaParametros = linspace(0, 10*M, Np)

#Intensidade observada em cada parametro de impacto
ObsInt = zeros(len(ListaParametros))
MaxObsInt = 0 #Intensidade maxima observada
bmax = 0
k = 0
for b in ListaParametros:
    #valores iniciais
    r0 = sqrt(b**2 + d**2)
    u0 = 1/r0
    phi0 = arctan(b/d)
    y0 = (1/(r0**2))*sqrt(1-(1-((2*M)/r0))*((b**2)/(r0**2)))
    z0 = u0**2*b*E
    
    sol = solve_ivp(F, [0, 150], [u0, y0, phi0, z0], events=(EventHorizon, LimitInterval), dense_output=True, max_step=0.1)
    
    if len(sol.t_events[0]) != 0:
        LimitPos = array([searchsorted(sol.t, sol.t_events[1][0]), len(sol.t) - 1])
    else:
        LimitPos = searchsorted(sol.t, [sol.t_events[1][0], sol.t_events[1][1]])
    
    if Emission(1/sol.sol(sol.t[int(LimitPos[0])])[0], sol.sol(sol.t[int(LimitPos[0])])[2])!=0:
        print("erro1", k, b)
        k += 1
        continue
    else:
        state = 0
        si = []
        sf = []
    
    #define os intervalos de integração
    for i in arange(LimitPos[0], LimitPos[1]):
        if Emission(1/sol.sol(sol.t[i])[0], sol.sol(sol.t[i])[2]) != state:
            if state==0:
                si.append(IntegrationInterval(i,state))
                state = 1
            else:
                sf.append(IntegrationInterval(i,state))
                state = 0
                
    #Caso ainda nao tenha nenhum intervalo de integração, verifica se perdeu algo no final
    if len(si)==0:
        if len(sol.t_events[0])!=0:
            sp = sol.t_events[0][0]
        else:
            sp = sol.t_events[1][1]
        x = linspace(sol.t[LimitPos[1]-1], sp, Nf)
        for s in x:
            if Emission(1/sol.sol(s)[0], sol.sol(s)[2])!=0:
                si.append(s)
                sf.append(sp)
                break
    
    if len(si)!=len(sf):
        if len(sol.t_events[0])!=0:#Entrou no Horizonte de Eventos
            if Emission(1/sol.sol(sol.t_events[0][0])[0], sol.sol(sol.t_events[0][0])[2])!=0:
                sf.append(sol.t_events[0][0] - 0.1)
            else:
                x = linspace(sol.t[LimitPos[1]-1], sol.t_events[0][0], Nf)
                sp = x[0]
                for s in x:
                    if Emission(1/sol.sol(s)[0], sol.sol(s)[2])==0:
                        sf.append(sp)
                        break
                    sp = s
        elif Emission(1/sol.sol(sol.t_events[1][1])[0], sol.sol(sol.t_events[1][1])[2])!=0:
            sf.append(sol.t_events[1][1])
        else:
            sf.append(IntegrationInterval(LimitPos[1],state))
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
                
                Delta = Emission(1/sol.sol(s)[0], sol.sol(s)[2]) - Absortion(1/sol.sol(s)[0], sol.sol(s)[2])*pow(g, 3/2)*ObsInt[k]
            
                #atualiza o red-shift
                ObsInt[k] *= g**2

                ObsInt[k] += h*pow(g,1/2)*Delta*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
            else:#caso seja o primeiro ponto
                ObsInt[k] += h*Emission(1/sol.sol(s)[0], sol.sol(s)[2])*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)
            
            rp = 1/sol.sol(s)[0]
            s -= h
            
    #realiza o red-shift para o infinito
    ObsInt[k] *= (1 - 2*M/rp)**2
    
    if ObsInt[k] > MaxObsInt:
        MaxObsInt = ObsInt[k]
        bmax = b
    print(k, b, ObsInt[k])
    k += 1


#plota o perfil do disco de acreção
AccretionDisk = empty([Nd,Nd], int)
for i in range(Nd):
    for j in range(Nd):
        x = (20/Nd)*i - 10
        y = (20/Nd)*j - 10
        AccretionDisk[j][i] = Emission(sqrt(x**2 + y**2), arctan2(y,x))
fig, ax = plt.subplots()
plt.pcolormesh(linspace(-10,10,Nd), linspace(-10,10,Nd), AccretionDisk)
circle1 = plt.Circle((0, 0), 2*M, color='black')
ax.add_patch(circle1)
plt.axis("square")
plt.show()

#Plota o perfil de intensidade observado
fig, ax = plt.subplots()
ax.plot(ListaParametros, ObsInt,'-b' , color='blue', label='Observed profile')
plt.xlim(0,10*M)
plt.ylim(0, 1.1*MaxObsInt)
plt.xlabel("b/M")
plt.ylabel("Observed Intensity")

#Plota a imagem do buraco negro observada
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
azm = linspace(0, 2*pi, 100)
r, th = meshgrid(ListaParametros, azm)
z = tile(ObsInt, (r.shape[0], 1))

plt.pcolormesh(th, r, z)
plt.clim(0, 1.1*MaxObsInt)
plt.colorbar(label='Observed Intensity')

print(MaxObsInt, bmax)
print(time() - t0)