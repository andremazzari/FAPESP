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
Ta = 0 #Taxa de absorção
Ng = 10 #Numero de pontos em cada integração gaussiana
Nd = 501 #Numero de divisoes no grid para plotar o disco de acreção
Np = 600 #Numero de parametros de impacto utilizados
Nf = 100 #Numero de divisoes no ajusto fino do intervalo de integração
NumAngulos = 80 #Numero de angulos que serão utilizados
alpha = pi/2 #Angulo de inclinação do disco

#funçao para integração gaussiana
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

#define a função que representa o lado direito do sistema
#Y[0]=u , Y[1]=y=d/ds(u) , Y[2]=phi, Y[3]=z=d/ds(phi)
def F(s, Y):
    return [Y[1], (2 - (3*Y[0]**2 - 7*M*Y[0]**3)*b**2)*E*Y[0]**3, Y[3], 2*Y[0]*Y[1]*b*E]

#define evento de atravessar o horizonte de eventos
def EventHorizon(s, Y):
    return 1/Y[0] - 2*M
EventHorizon.terminal = True #define este evento como terminal

#Versao 1 do perfil de acreçao(com algum problema!!!)
def Disk1(r, phi, theta):
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)
    
    flag = 0
    if abs(z*sin(alpha) + y*cos(alpha)) > 1e-7:
        cosB = 1/sqrt(1 + (x/(z*sin(alpha) + y*cos(alpha)))**2)
        sinB = x*cosB/(z*sin(alpha) + y*cos(alpha))
    elif abs(x)>1e-7:
        print(abs(z*sin(alpha) + y*cos(alpha)))
        print(x)
        if ((z*sin(alpha) + y*cos(alpha))/x) > 0:
            sinB = 1/sqrt(1 + ((z*sin(alpha) + y*cos(alpha))/x)**2)
        else:
            sinB = -1/sqrt(1 + ((z*sin(alpha) + y*cos(alpha))/x)**2)
        cosB = sinB*(z*sin(alpha) + y*cos(alpha))/x
    else:
        flag = 1
        rp = r
        phip = phi
    
    if flag == 0:
        if cosB>1 or sinB>1:
            print("erro5",b,theta)
    
        yp = -x*cos(alpha)*sinB + y*(cosB + (sin(alpha)**2)*(1 - cosB)) + z*sin(alpha)*cos(alpha)*(1 - cosB)
        zp = x*sin(alpha)*sinB + y*sin(alpha)*cos(alpha)*(1 - cosB) + z*(cosB + (cos(alpha)**2)*(1 - cosB))

        rp = sqrt(yp**2 + zp**2)
        phip = arctan2(yp, zp) - alpha
    
    if abs(rp*sin(phip)) > 2*(rp*cos(phip))**2 + 2:
        return 1
    else:
        return 0

#versao 2 do perfil de acreção
def DiskV2(r, phi, theta):
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
    
    if abs(r*CosPhi) < 5*M:
        return 1
    else:
        return 0

#Evento em que o raio entra e sai do circulo com r=10sqrt(2)M
def LimitInterval(s, Y):
    return 1/Y[0] - 10*sqrt(2)*M

def IntegrationInterval(i, state, theta):
    x = linspace(sol.t[i-1], sol.t[i], Nf)
    if state==0: #Entrando no disco
        for s in x:
            if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)!=0:
                return s
    else: #Saindo do disco
        sp = x[0]
        for s in x:
            if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)==0:
                return sp
            sp = s

def Integrand(s, theta):
    if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)==0:
        print("errro3",k,b,s,i,j,theta,m)
    return DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)*sqrt((1/(1 - 2*M*sol.sol(s)[0]))*(sol.sol(s)[1]/(sol.sol(s)[0]**2))**2 + (sol.sol(s)[3]/sol.sol(s)[0])**2)*(1 - 2*M*sol.sol(s)[0])**2

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
        
        if DiskV2(1/sol.sol(sol.t[int(LimitPos[0])])[0], sol.sol(sol.t[int(LimitPos[0])])[2], theta)!=0:
            state = 1
            si = [sol.t[int(LimitPos[0])]] #TALVEZ PRECISE DE AJUSTE FINO AQUI
            sf = []
        else:
            state = 0
            si = []
            sf = []
        
        #define os intervalos de integração
        for i in arange(LimitPos[0], LimitPos[1]):
            if DiskV2(1/sol.sol(sol.t[i])[0], sol.sol(sol.t[i])[2], theta) != state:
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
                if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)!=0:
                    si.append(s)
                    sf.append(sp)
                    break
        
        if len(si)!=len(sf):
            if len(sol.t_events[0])!=0:
                if DiskV2(1/sol.sol(sol.t_events[0][0])[0], sol.sol(sol.t_events[0][0])[2], theta)!=0:
                    sf.append(sol.t_events[0][0])
                else:
                    x = linspace(sol.t[LimitPos[1]-1], sol.t_events[0][0], Nf)
                    sp = x[0]
                    for s in x:
                        if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], theta)==0:
                            sf.append(sp)
                            break
                        sp = s
            elif DiskV2(1/sol.sol(sol.t_events[1][1])[0], sol.sol(sol.t_events[1][1])[2], theta)!=0:
                sf.append(sol.t_events[1][1])
            else:
                sf.append(IntegrationInterval(LimitPos[1],state,theta))
        if len(si)!=len(sf):
            print("erro2", k, b)
            print(si)
            print(sf)
            k += 1
            continue
        
        #realiza a integração
        for i in flip(arange(len(si))):
            x, w = gaussxw(Ng)
            xp = 0.5*(sf[i] - si[i])*x + 0.5*(sf[i] + si[i])
            wp = 0.5*(sf[i] - si[i])*w
            
            for j in range(Ng):
                ObsInt[m][k] = (1 - Ta)*ObsInt[m][k]
                ObsInt[m][k] += wp[j]*Integrand(xp[j], theta)
        
        #segundo quadrante é igual
        if theta!=0 and theta!=pi/2:
            ObsInt[2*NumAngulos-m-2][k] = ObsInt[m][k]
            
        #realiza o mesmo para os terceiro e quarto quadrantes
        
        l = 2*NumAngulos - 2 + m
        
        if DiskV2(1/sol.sol(sol.t[int(LimitPos[0])])[0], sol.sol(sol.t[int(LimitPos[0])])[2], -theta)!=0:
            state = 1
            si = [sol.t[int(LimitPos[0])]] #TALVEZ PRECISE DE AJUSTE FINO AQUI
            sf = []
        else:
            state = 0
            si = []
            sf = []
        
        #define os intervalos de integração
        for i in arange(LimitPos[0], LimitPos[1]):
            if DiskV2(1/sol.sol(sol.t[i])[0], sol.sol(sol.t[i])[2], -theta) != state:
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
                if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], -theta)!=0:
                    si.append(s)
                    sf.append(sp)
                    break
        
        if len(si)!=len(sf):
            if len(sol.t_events[0])!=0:
                if DiskV2(1/sol.sol(sol.t_events[0][0])[0], sol.sol(sol.t_events[0][0])[2], -theta)!=0:
                    sf.append(sol.t_events[0][0])
                else:
                    x = linspace(sol.t[LimitPos[1]-1], sol.t_events[0][0], Nf)
                    sp = x[0]
                    for s in x:
                        if DiskV2(1/sol.sol(s)[0], sol.sol(s)[2], -theta)==0:
                            sf.append(sp)
                            break
                        sp = s
            elif DiskV2(1/sol.sol(sol.t_events[1][1])[0], sol.sol(sol.t_events[1][1])[2], -theta)!=0:
                sf.append(sol.t_events[1][1])
            else:
                sf.append(IntegrationInterval(LimitPos[1],state, -theta))
        if len(si)!=len(sf):
            print("erro2", k, b)
            print(si)
            print(sf)
            k += 1
            continue
        
        #realiza a integração
        for i in flip(arange(len(si))):
            x, w = gaussxw(Ng)
            xp = 0.5*(sf[i] - si[i])*x + 0.5*(sf[i] + si[i])
            wp = 0.5*(sf[i] - si[i])*w
            
            for j in range(Ng):
                ObsInt[l][k] = (1 - Ta)*ObsInt[l][k]
                ObsInt[l][k] += wp[j]*Integrand(xp[j], -theta)
                
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
        AccretionDisk[j][i] = DiskV2(sqrt(x**2 + y**2), arctan2(y,x), pi/2)
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