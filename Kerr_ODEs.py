from numpy import cos,sin,sqrt,array

M = 1
a = 0.9
reh = M + sqrt(M**2 - a**2)

def gamma_11(x):
    sigma = x[1]**2 + (a*cos(x[2]))**2
    Delta = x[1]**2 + -2*M*x[1] + a**2
    return sigma/Delta

def gamma_12(x):
    return 0

def gamma_13(x):
    return 0

def gamma_21(x):
    return 0

def gamma_22(x):
    sigma = x[1]**2 + (a*cos(x[2]))**2
    return sigma

def gamma_23(x):
    return 0

def gamma_31(x):
    return 0

def gamma_32(x):
    return 0

def gamma_33(x):
    sigma = x[1]**2 + (a*cos(x[2]))**2
    Delta = x[1]**2 + -2*M*x[1] + a**2
    return ((x[1]**2 + a**2)**2 - Delta*(a*sin(x[2]))**2)*(sin(x[2])**2/sigma)

def inv_gamma11(x):
    return 1/gamma_11(x)

def inv_gamma12(x):
    return 0

def inv_gamma13(x):
    return 0

def inv_gamma21(x):
    return 0

def inv_gamma22(x):
    return 1/gamma_22(x)

def inv_gamma23(x):
    return 0

def inv_gamma31(x):
    return 0

def inv_gamma32(x):
    return 0

def inv_gamma33(x):
    return 1/gamma_33(x)

def D1inv_gamma11(x):
    return 2*(-x[1]*(-2*M*x[1] + a**2 + x[1]**2) + (-M + x[1])*(a**2*cos(x[2])**2 + x[1]**2))/(a**2*cos(x[2])**2 + x[1]**2)**2

def D2inv_gamma11(x):
    return 4*a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(2*x[2])/(a**2*cos(2*x[2]) + a**2 + 2*x[1]**2)**2

def D3inv_gamma11(x):
    return 0

def D1inv_gamma12(x):
    return 0

def D2inv_gamma12(x):
    return 0
    
def D3inv_gamma12(x):
    return 0

def D1inv_gamma13(x):
    return 0

def D2inv_gamma13(x):
    return 0
    
def D3inv_gamma13(x):
    return 0

def D1inv_gamma21(x):
    return 0

def D2inv_gamma21(x):
    return 0
    
def D3inv_gamma21(x):
    return 0

def D1inv_gamma22(x):
    return -2*x[1]/(a**2*cos(x[2])**2 + x[1]**2)**2

def D2inv_gamma22(x):
    return 4*a**2*sin(2*x[2])/(a**2*cos(2*x[2]) + a**2 + 2*x[1]**2)**2

def D3inv_gamma22(x):
    return 0

def D1inv_gamma23(x):
    return 0

def D2inv_gamma23(x):
    return 0
    
def D3inv_gamma23(x):
    return 0

def D1inv_gamma31(x):
    return 0

def D2inv_gamma31(x):
    return 0
    
def D3inv_gamma31(x):
    return 0

def D1inv_gamma32(x):
    return 0

def D2inv_gamma32(x):
    return 0
    
def D3inv_gamma32(x):
    return 0

def D1inv_gamma33(x):
    return -(2*x[1]*(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2) + 2*(a**2*cos(x[2])**2 + x[1]**2)*(a**2*(M - x[1])*sin(x[2])**2 + 2*x[1]*(a**2 + x[1]**2)))/((a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2*sin(x[2])**2)

def D2inv_gamma33(x):
    return 2*(a**2*(a**2*cos(x[2])**2 + x[1]**2)*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 + a**2*(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)*sin(x[2])**2 + (a**2*cos(x[2])**2 + x[1]**2)*(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2))*cos(x[2])/((a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2*sin(x[2])**3)

def D3inv_gamma33(x):
    return 0

def beta1(x):
    return 0

def beta2(x):
    return 0

def beta3(x):
    Delta = x[1]**2 + -2*M*x[1] + a**2
    return -2*M*x[1]*a/((x[1]**2 + a**2)**2 - Delta*(a*sin(x[2]))**2)

def D1beta1(x):
    return 0

def D2beta1(x):
    return 0

def D3beta1(x):
    return 0

def D1beta2(x):
    return 0

def D2beta2(x):
    return 0

def D3beta2(x):
    return 0

def D1beta3(x):
    return 2*M*a*(-a**4*cos(x[2])**2 + a**2*x[1]**2*cos(x[2])**2 + a**2*x[1]**2 + 3*x[1]**4)/(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2

def D2beta3(x):
    return -4*M*a**3*x[1]*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])*cos(x[2])/(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2

def D3beta3(x):
    return 0

def alpha(x):
    sigma = x[1]**2 + (a*cos(x[2]))**2
    Delta = x[1]**2 + -2*M*x[1] + a**2
    gtt = - (1 - 2*M*x[1]/sigma)
    
    alpha = (2*M*x[1]*a)/((x[1]**2 + a**2)**2 - Delta*(a*sin(x[2]))**2)
    alpha *= 2*M*x[1]*a*sin(x[2])**2
    alpha /= sigma
    alpha -= gtt
    alpha = sqrt(alpha)
    return alpha

def D1alpha(x):
    return M*(4*M*a**2*x[1]**3*cos(x[2])**2 - 4*M*a**2*x[1]**3 - a**6*cos(x[2])**2 - 2*a**4*x[1]**2*cos(x[2])**2 + a**4*x[1]**2 - a**2*x[1]**4*cos(x[2])**2 + 2*a**2*x[1]**4 + x[1]**6)/(sqrt((a**2*cos(x[2])**2 + x[1]**2)*(2*M*x[1] - a**2 - x[1]**2)/(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2))*(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2)

def D2alpha(x):
    return 2*M*a**2*x[1]*(2*M*a**2*x[1] + 2*M*x[1]**3 - a**4 - 2*a**2*x[1]**2 - x[1]**4)*sin(x[2])*cos(x[2])/(sqrt((a**2*cos(x[2])**2 + x[1]**2)*(2*M*x[1] - a**2 - x[1]**2)/(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2))*(a**2*(-2*M*x[1] + a**2 + x[1]**2)*sin(x[2])**2 - (a**2 + x[1]**2)**2)**2)

def D3alpha(x):
    return 0

def Fx1(x,u):
    if inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3] < 0:
        print("erro aqui!:",inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    return (inv_gamma11(x)*u[1] + inv_gamma12(x)*u[2] + inv_gamma13(x)*u[3])/u0 - beta1(x)

def Fx2(x,u):
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    return (inv_gamma21(x)*u[1] + inv_gamma22(x)*u[2] + inv_gamma23(x)*u[3])/u0 - beta2(x)

def Fx3(x,u):
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    return (inv_gamma31(x)*u[1] + inv_gamma32(x)*u[2] + inv_gamma33(x)*u[3])/u0 - beta3(x)

def Fu1(x,u):
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    F = -alpha(x)*u0*D1alpha(x) + u[1]*D1beta1(x) + u[2]*D1beta2(x) + u[3]*D1beta3(x)
    F -= u[1]*(u[1]*D1inv_gamma11(x) + u[2]*D1inv_gamma12(x) + u[3]*D1inv_gamma13(x))/(2*u0)
    F -= u[2]*(u[1]*D1inv_gamma21(x) + u[2]*D1inv_gamma22(x) + u[3]*D1inv_gamma23(x))/(2*u0)
    F -= u[3]*(u[1]*D1inv_gamma31(x) + u[2]*D1inv_gamma32(x) + u[3]*D1inv_gamma33(x))/(2*u0)
    return F

def Fu2(x,u):
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    F = -alpha(x)*u0*D2alpha(x) + u[1]*D2beta1(x) + u[2]*D2beta2(x) + u[3]*D2beta3(x)
    F -= u[1]*(u[1]*D2inv_gamma11(x) + u[2]*D2inv_gamma12(x) + u[3]*D2inv_gamma13(x))/(2*u0)
    F -= u[2]*(u[1]*D2inv_gamma21(x) + u[2]*D2inv_gamma22(x) + u[3]*D2inv_gamma23(x))/(2*u0)
    F -= u[3]*(u[1]*D2inv_gamma31(x) + u[2]*D2inv_gamma32(x) + u[3]*D2inv_gamma33(x))/(2*u0)
    return F

def Fu3(x,u):
    u0 = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma21(x)*u[2]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma23(x)*u[2]*u[3] + inv_gamma31(x)*u[3]*u[1] + inv_gamma32(x)*u[3]*u[2] + inv_gamma33(x)*u[3]*u[3])
    u0 /= alpha(x)
    F = -alpha(x)*u0*D3alpha(x) + u[1]*D3beta1(x) + u[2]*D3beta2(x) + u[3]*D3beta3(x)
    F -= u[1]*(u[1]*D3inv_gamma11(x) + u[2]*D3inv_gamma12(x) + u[3]*D3inv_gamma13(x))/(2*u0)
    F -= u[2]*(u[1]*D3inv_gamma21(x) + u[2]*D3inv_gamma22(x) + u[3]*D3inv_gamma23(x))/(2*u0)
    F -= u[3]*(u[1]*D3inv_gamma31(x) + u[2]*D3inv_gamma32(x) + u[3]*D3inv_gamma33(x))/(2*u0)
    return F

def Runge_Kutta_4_step(x,u,dt):
    k1x = dt*array([dt, Fx1(x,u), Fx2(x,u), Fx3(x,u)], float)
    k1u = dt*array([0, Fu1(x,u), Fu2(x,u), Fu3(x,u)], float)
    
    k2x = dt*array([dt, Fx1(x + k1x/2,u + k1u/2), Fx2(x + k1x/2,u + k1u/2), Fx3(x + k1x/2,u + k1u/2)], float)
    k2u = dt*array([0, Fu1(x + k1x/2,u + k1u/2), Fu2(x + k1x/2,u + k1u/2), Fu3(x + k1x/2,u + k1u/2)], float)
    
    k3x = dt*array([dt, Fx1(x + k2x/2,u + k2u/2), Fx2(x + k2x/2,u + k2u/2), Fx3(x + k2x/2,u + k2u/2)], float)
    k3u = dt*array([0, Fu1(x + k2x/2,u + k2u/2), Fu2(x + k2x/2,u + k2u/2), Fu3(x + k2x/2,u + k2u/2)], float)
    
    k4x = dt*array([dt, Fx1(x + k3x,u + k3u), Fx2(x + k3x,u + k3u), Fx3(x + k3x,u + k3u)], float)
    k4u = dt*array([0, Fu1(x + k3x/2,u + k3u), Fu2(x + k3x,u + k3u), Fu3(x + k3x,u + k3u)], float)
    
    dx = (1/6)*(k1x + 2*k2x + 2*k3x + k4x)
    du = (1/6)*(k1u + 2*k2u + 2*k3u + k4u)
    
    return dx, du

def Runge_Kutta_45_step(x,u,dt):
    k1x = dt*array([dt, Fx1(x,u), Fx2(x,u), Fx3(x,u)], float)
    k1u = dt*array([0, Fu1(x,u), Fu2(x,u), Fu3(x,u)], float)
    
    if (x + k1x/4)[1] <= (reh + 1e-2): #Horizonte de eventos
        return 0,0,1
        
    k2x = dt*array([dt, Fx1(x + k1x/4,u + k1u/4), Fx2(x + k1x/4,u + k1u/4), Fx3(x + k1x/4,u + k1u/4)], float)
    k2u = dt*array([0, Fu1(x + k1x/4,u + k1u/4), Fu2(x + k1x/4,u + k1u/4), Fu3(x + k1x/4,u + k1u/4)], float)
    
    if (x + (3/32)*k1x + k2x*(9/32))[1] <= (reh + 1e-2): #Horizonte de eventos
        return 0,0,1
    
    k3x = dt*array([dt, Fx1(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32)), Fx2(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32)), Fx3(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32))], float)
    k3u = dt*array([0, Fu1(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32)), Fu2(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32)), Fu3(x + (3/32)*k1x + k2x*(9/32),u + (3/32)*k1u + k2u*(9/32))], float)
    
    if (x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x)[1] <= (reh + 1e-2): #Horizonte de eventos
        return 0,0,1
    
    k4x = dt*array([dt, Fx1(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u), Fx2(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u), Fx3(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u)], float)
    k4u = dt*array([0, Fu1(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u), Fu2(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u), Fu3(x + (1932/2197)*k1x - (7200/2197)*k2x + (7296/2197)*k3x,u + (1932/2197)*k1u - (7200/2197)*k2u + (7296/2197)*k3u)], float)
    
    if (x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x)[1] <= (reh + 1e-2): #Horizonte de eventos
        return 0,0,1
    
    k5x = dt*array([dt, Fx1(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u), Fx2(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u), Fx3(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u)], float)
    k5u = dt*array([0, Fu1(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u), Fu2(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u), Fu3(x + (439/216)*k1x - 8*k2x + (3680/513)*k3x - (845/4104)*k4x,u + (439/216)*k1u - 8*k2u + (3680/513)*k3u - (845/4104)*k4u)], float)
    
    if (x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x)[1] <= (reh + 1e-2): #Horizonte de eventos
        return 0,0,1
    
    k6x = dt*array([dt, Fx1(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u), Fx2(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u), Fx3(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u)], float)
    k6u = dt*array([0, Fu1(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u), Fu2(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u), Fu3(x - (8/27)*k1x + 2*k2x - (3544/2565)*k3x + (1859/4104)*k4x - (11/40)*k5x,u - (8/27)*k1u + 2*k2u - (3544/2565)*k3u + (1859/4104)*k4u - (11/40)*k5u)], float)
    
    dx = (16/135)*k1x + (6656/12825)*k3x + (28561/56430)*k4x - (9/50)*k5x + (2/55)*k6x
    du = (16/135)*k1u + (6656/12825)*k3u + (28561/56430)*k4u - (9/50)*k5u + (2/55)*k6u
    
    return dx, du, 0
    
    

def Kerr_conserved_quantities(x,u):
    E = -u[0]
    #componente z do momento angular
    Lz = u[3]
    #constante de carter
    Q = u[2]**2 + cos(x[2])**2*(-(a*E)**2 + (Lz/sin(x[2]))**2)
    #Hamiltoniano
    H = alpha(x)*sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma33(x)*u[3]*u[3] + 2*inv_gamma12(x)*u[1]*u[2] + 2*inv_gamma13(x)*u[1]*u[3] + 2*inv_gamma23(x)*u[2]*u[3]) - beta1(x)*u[1] - beta2(x)*u[2] - beta3(x)*u[3]
    
    return Lz, Q, H

def inv_gamma(x,i,j):
    if i == 1:
        if j == 1:
            return inv_gamma11(x)
        elif j == 2:
            return inv_gamma12(x)
        elif j == 3:
            return inv_gamma13(x)
    elif i == 2:
        if j == 1:
            return inv_gamma21(x)
        elif j == 2:
            return inv_gamma22(x)
        elif j == 3:
            return inv_gamma23(x)
    elif i == 3:
        if j == 1:
            return inv_gamma31(x)
        elif j == 2:
            return inv_gamma32(x)
        elif j == 3:
            return inv_gamma33(x)

def Delta_Hamiltonian_term_x(x,u,up,i):
    if i == 1:
        result = inv_gamma11(x)*(up + u[1]) + 2*inv_gamma12(x)*u[2] + 2*inv_gamma13(x)*u[3]
    elif i == 2:
        result = 2*inv_gamma21(x)*u[1] + inv_gamma22(x)*(up + u[2]) + 2*inv_gamma23(x)*u[3]
    elif i == 3:
        result = 2*inv_gamma31(x)*u[1] + 2*inv_gamma32(x)*u[2] + inv_gamma33(x)*(up + u[3])
    
    U = inv_gamma11(x)*u[1]**2 + inv_gamma22(x)*u[2]**2 + inv_gamma33(x)*u[3]**2 + 2*(inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma23(x)*u[2]*u[3])
    
    if i == 1:
        Up = inv_gamma11(x)*up**2 + inv_gamma22(x)*u[2]**2 + inv_gamma33(x)*u[3]**2 + 2*(inv_gamma12(x)*up*u[2] + inv_gamma13(x)*up*u[3] + inv_gamma23(x)*u[2]*u[3])
    elif i == 2:
        Up = inv_gamma11(x)*u[1]**2 + inv_gamma22(x)*up**2 + inv_gamma33(x)*u[3]**2 + 2*(inv_gamma12(x)*u[1]*up + inv_gamma13(x)*u[1]*u[3] + inv_gamma23(x)*up*u[3])
    elif i == 3:
        Up = inv_gamma11(x)*u[1]**2 + inv_gamma22(x)*u[2]**2 + inv_gamma33(x)*up**2 + 2*(inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*up + inv_gamma23(x)*u[2]*up)
        
    result /= sqrt(Up) + sqrt(U)
    
    if i == 1:
        result -= beta1(x)
    elif i == 2:
        result -= beta2(x)
    elif i == 3:
        result -= beta3(x)
        
    return result

def Delta_Hamiltonian_term_u(x,u,xp,i):
    if i == 1:
        xp = array([x[0],xp,x[2],x[3]], float)
    elif i == 2:
        xp = array([x[0],x[1],xp,x[3]], float)
    elif i == 3:
        xp = array([x[0],x[1],x[2],xp], float)
    
    result = sqrt(inv_gamma11(x)*u[1]*u[1] + inv_gamma22(x)*u[2]*u[2] + inv_gamma33(x)*u[3]*u[3] + 2*(inv_gamma12(x)*u[1]*u[2] + inv_gamma13(x)*u[1]*u[3] + inv_gamma23(x)*u[2]*u[3]))
    
    result += sqrt(inv_gamma11(xp)*u[1]*u[1] + inv_gamma22(xp)*u[2]*u[2] + inv_gamma33(xp)*u[3]*u[3] + 2*(inv_gamma12(xp)*u[1]*u[2] + inv_gamma13(xp)*u[1]*u[3] + inv_gamma23(xp)*u[2]*u[3]))
    
    temp = result
    
    if i == 1:
        result *= -0.5*D1alpha(x)
    elif i == 2:
        result *= -0.5*D2alpha(x)
    elif i == 3:
        result *= -0.5*D3alpha(x)
        
    if i == 1:
        temp2 = D1inv_gamma11(x)*u[1]*u[1] + D1inv_gamma22(x)*u[2]*u[2] + D1inv_gamma33(x)*u[3]*u[3] + 2*(D1inv_gamma12(x)*u[1]*u[2] + D1inv_gamma13(x)*u[1]*u[3] + D1inv_gamma23(x)*u[2]*u[3])
    elif i == 2:
        temp2 = D2inv_gamma11(x)*u[1]*u[1] + D2inv_gamma22(x)*u[2]*u[2] + D2inv_gamma33(x)*u[3]*u[3] + 2*(D2inv_gamma12(x)*u[1]*u[2] + D2inv_gamma13(x)*u[1]*u[3] + D2inv_gamma23(x)*u[2]*u[3])
    elif i == 3:
        temp2 = temp2 = D2inv_gamma11(x)*u[1]*u[1] + D2inv_gamma22(x)*u[2]*u[2] + D2inv_gamma33(x)*u[3]*u[3] + 2*(D2inv_gamma12(x)*u[1]*u[2] + D2inv_gamma13(x)*u[1]*u[3] + D2inv_gamma23(x)*u[2]*u[3])
    
    result += -0.5*(alpha(xp) + alpha(x))*temp2/temp
    
    if i == 1:
        result += D1beta1(x)*u[1] + D1beta2(x)*u[2] + D1beta3(x)*u[3]
    elif i == 2:
        result += D2beta1(x)*u[1] + D2beta2(x)*u[2] + D2beta3(x)*u[3]
    elif i == 3:
        result += D3beta1(x)*u[1] + D3beta2(x)*u[2] + D3beta3(x)*u[3]
        
    return result

def x1_Hamilton_eq(x,u,xp,up,dt):
    temp = Delta_Hamiltonian_term_x([x[0],xp[1],x[2],x[3]],u,up[1],1)
    temp += Delta_Hamiltonian_term_x(x,[u[0],u[1],up[2],up[3]],up[1],1)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],xp[3]],[u[0],u[1],up[2],up[3]],up[1],1)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],xp[3]],[u[0],u[1],u[2],up[3]],up[1],1)
    temp += Delta_Hamiltonian_term_x(x,u,up[1],1)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],x[2],xp[3]],[u[0],u[1],u[2],up[3]],up[1],1)
    temp *= 1/6
    
    return (xp[1] - x[1])/dt - temp

def x2_Hamilton_eq(x,u,xp,up,dt): 
    temp = Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],x[3]],[u[0],up[1],u[2],u[3]],up[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],x[3]],u,up[2],2)
    temp += Delta_Hamiltonian_term_x(x,u,up[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],xp[3]],[u[0],up[1],u[2],up[3]],up[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],x[3]],[u[0],up[1],u[2],u[3]],up[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],xp[3]],[u[0],up[1],up[2],up[3]],up[2],2)
    temp *= 1/6
    
    return (xp[2] - x[2])/dt - temp

def x3_Hamilton_eq(x,u,xp,up,dt): 
    temp = Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],xp[3]],[u[0],up[1],up[2],u[3]],up[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],xp[3]],[u[0],u[1],up[2],u[3]],up[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],x[3]],[u[0],u[1],up[2],u[3]],up[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],x[2],xp[3]],u,up[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],x[3]],[u[0],up[1],up[2],u[3]],up[3],3)
    temp += Delta_Hamiltonian_term_x(x,u,up[3],3)
    temp *= 1/6
    
    return (xp[3] - x[3])/dt - temp

def u1_Hamilton_eq(x,u,xp,up,dt):
    temp = Delta_Hamiltonian_term_u(x,u,xp[1],1)
    temp += Delta_Hamiltonian_term_u([x[0],x[1],xp[2],xp[3]],[u[0],u[1],up[2],up[3]],xp[1],1)
    temp += Delta_Hamiltonian_term_u([x[0],x[1],xp[2],xp[3]],[u[0],up[1],up[2],up[3]],xp[1],1)
    temp += Delta_Hamiltonian_term_u([x[0],x[1],x[2],xp[3]],[u[0],u[1],u[2],up[3]],xp[1],1)
    temp += Delta_Hamiltonian_term_u(x,[u[0],up[1],u[2],u[3]],xp[1],1)
    temp += Delta_Hamiltonian_term_u([x[0],x[1],x[2],xp[3]],[u[0],up[1],u[2],up[3]],xp[1],1)
    temp *= 1/6
    
    return (xp[1] - x[1])/dt - temp

def u2_Hamilton_eq(x,u,xp,up,dt): 
    temp = Delta_Hamiltonian_term_x([x[0],xp[1],x[2],x[3]],[u[0],up[1],u[2],u[3]],xp[2],2)
    temp += Delta_Hamiltonian_term_x(x,u,xp[2],2)
    temp += Delta_Hamiltonian_term_x(x,[u[0],u[1],up[2],u[3]],xp[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],xp[3]],[u[0],up[1],u[2],up[3]],xp[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],x[3]],[u[0],up[1],up[2],u[3]],xp[2],2)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],x[2],xp[3]],[u[0],up[1],up[2],up[3]],xp[2],2)
    temp *= 1/6
    
    return (xp[2] - x[2])/dt - temp

def u3_Hamilton_eq(x,u,xp,up,dt): 
    temp = Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],x[3]],[u[0],up[1],up[2],u[3]],xp[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],x[3]],[u[0],u[1],up[2],u[3]],xp[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],x[1],xp[2],x[3]],[u[0],u[1],up[2],up[3]],xp[3],3)
    temp += Delta_Hamiltonian_term_x(x,u,xp[3],3)
    temp += Delta_Hamiltonian_term_x([x[0],xp[1],xp[2],x[3]],[u[0],up[1],up[2],up[3]],xp[3],3)
    temp += Delta_Hamiltonian_term_x(x,[u[0],u[1],u[2],up[3]],xp[3],3)
    temp *= 1/6
    
    return (xp[3] - x[3])/dt - temp