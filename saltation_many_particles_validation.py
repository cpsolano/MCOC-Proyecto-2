from matplotlib.pylab import *
#from mpl_toolkits.mplot3d import Axes3D
import random , math

#unidades SI
_m= 1.
_kg = 1.
_s = 1.
_mm = 1e-3*_m
_gr = 1e-3*_kg
#------------------

vfx = 5.0*_m/_s    #m/s 
vfy = 0.1 *_m/_s   #m/s


#Parametros
g = 9.81*_m/_s**2
d = 15 * _mm

rho_particula = 2650.*_kg/(_m**3)
rho_agua = 1000.*_kg/(_m**3)

Cd = 0.47   #Constante de drag
Cl = 0.2    #Constante de lifting
Cm = 0.5    #Constante de peso virtual
R = (rho_particula/rho_agua) -1
alpha = 1/(1+R+Cm)
ustar = 0.18*_m/_s # sqrt(tauw/rho_agua)

A = pi*(d/2)**2           #Area transversal de particula
V = (4./3.)*pi*(d/2)**3   #Volumen de particula
m = rho_particula*V       #masa de la particula, grano de arena
#------------------

#Tiempo
dt = 0.001*_s   #paso de tiempo
tmax = 1 *_s   #tiempo maximo de simulacion
ti = 0.*_s      #tiempo actual
t = arange(0,tmax,dt)
Nt = len(t)
#Nt = int32(2*tmax/ dt) #numero de espacios de tiempo, tiene que ser entero porque hare un for
#------------------

#Generacion de particulas
n = 20   #Numero de particulas

x01 = zeros((n,2))   #matriz posicion de las particulas
v01 = zeros((n,2))   #matriz velocidad de las particulas
for i in range(n):
	x01[i][0:2] = array([random.random(),random.random()])*10*d*5  #valores iniciales de la posicion
	v01[i][0:2] = array([random.random(),0.])/2  # y la velocidad de la particula

print"Condiciones iniciales:"
print "Posiociones ="
print x01
print "Velocidades ="
print v01
#------------------

ihat= array([1,0])
jhat= array([0,1])

#Fuerzas constantes
W = array([0, -m*g])            #Vector peso de la particula
fB = array([0, rho_agua*V*g])   #Vetcor Empuje
#------------------

#Funciones y fuerzas
norm = lambda v: sqrt(dot(v,v))   #Funcion norma

#dist = lambda x1,x2: sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2)   #Distancia entre 2 particulas (lo mismo que funcion norm)

def velocity_field(x):   #Perfil logaritmico
	z=x[1]/d
	if z>1./30.:
		uf= ustar*log(30.*z)/0.41
	else:
		uf=0

	return array([uf,0])

vfx = velocity_field([0, 4*d])[0]

k_penal= 1000*0.5*Cd*rho_agua*A*norm(vfx)**2/(1*_mm)   #Constante de resorte

#def UtUb(xy,vi):   #Velocidad relativa sobre y bajo la particula
#	xy = xy[1]
#	if xy > d/2:
#		return array([velocity_field(xy+d/2)-vi[0],velocity_field(xy-d/2)-vi[0]])
#	else:
#		return array([velocity_field(d)-vi[0],velocity_field(0)-vi[0]])


angulo = lambda x1,x2: math.atan2(x1[1]-x2[1],x1[0]-x2[0]) #Angulo de choque entre 2 particulas (Error cuando se divide por 0)

#def lifting(xi,vi):   #Fuerza de lifting
#	Utb = UtUb(xi,vi) 
#	return 3./4.*alfa*Cl*(norm(Utb[0])**2 - norm(Utb[1])**2)
	#return 0.5*rho_agua*(vi[0]-vfx(xi[1]))**2*A*Cl

def particula(z,t):

	zp = zeros(len(z))
	for i in range(n):
		xi = z[4*i:(4*i+2)]
		vi = z[4*i+2:(4*i+4)]

		vf = velocity_field(xi)
		vf_top = norm(velocity_field(xi + (d/2)*jhat)) #falta
		vf_bot = norm(velocity_field(xi - (d/2)*jhat)) #falta
		vrel = vf - vi
		vrel = vf-vi
		fD = (0.5*Cd*alpha*rho_agua*norm(vrel)*A)*vrel  #formula wiki
		#fD = alpha*(R*(d*g/(ustar**2))-(3./4.)*Cd*(vrel)*norm(vrel)) # formula PM
		fL = (0.5*Cl*alpha*rho_agua*(vf_top**2 -vf_bot**2)*A)*jhat #formula wiki
		#fL = alpha*(3/4*CL*(norm(vf_top)**2 - norm(vf_bot)**2))
		Fi = W + fD + fB + fL

		if xi[1] < d/2.:   #Cuando centro de la particula se encuentra a una distancia de d/2 se produce el choque contra el piso
			Fi[1] += -k_penal*(xi[1]-d/2)

		for j in range(n):   #Se revisa si hay colision con alguna de las otras particulas
			if i > j:
				x2 = z[j*4:2+j*4]    #Posicion y velocidad de la 2da particula
				xdif = norm(xi-x2)   #Distancia entre los centros de las particulas

				if xdif < d:         #Revisamos si hay choque revisando si la diferencia de posiciones es menor a 1 diametro
					Fct = abs(-k_penal*(xdif-d))   #La misma idea de choque contra el suelo
					theta = angulo(xi,x2)          #Angulo con que chocan las particulas
					if xi[0] > x2[0]:              #Revisamos posicion de la particula con respecto a la que choco para ver si
						Fcx = Fct*cos(theta)       # acelera o desacelera en el eje x e y
					else:
						Fcx = -Fct*cos(theta)
					if xi[1] >= x2[1]:
						Fcy = Fct*cos(theta)
					else:
						Fcy = -Fct*cos(theta)
					Fc = array([Fcx, Fcy])
					zp[2+i*4:4+i*4] += Fc/m
					zp[2+j*4:4+j*4] -= Fc/m
					#print "choque",t, p1, p2

		zp[i*4:2+i*4] += vi      #Guardamos la derivada de la posicion con respecto al tiempo
		zp[2+i*4:4+i*4] += Fi/m  #Guardamos la derivada de la velocidad con respecto al tiempo

	return zp   #Retornamos el vector derivada


#Simulacion
from scipy.integrate import odeint

z0 = zeros(4*n)      #Vector de condiciones iniciales para la integracion. Guardamos de la forma
for p in range(n):   # z0 = [x1,y1,vx1,vy1,x2,y2,vx2,vy2,x3,y3,vx3,vy3,....,xn,yn,vxn,vyn]
    z0[p*4:2+p*4] += x01[p]
    z0[2+p*4:4+p*4] += v01[p]


z = odeint(particula, z0, t)
#------------------

#Graficos
figure()
ax=gca()
for p in range(n):
    x1 = z[:,p*4:2+p*4]/d
    plot(x1[:,0],x1[:,1],label="p"+str(p+1))
   
ax.axhline(d/2,color="k", linestyle="--")
#plot([0,t],[0,0],label="piso")
#plot(x2[:,0],x2[:,1],label="x2")
#ylim([0,50*_mm])
plt.title("Posicion particulas plano XY")
plt.legend()


#figure()
#for p in range(n):
 #   subplot(2,n,p+1)
 #   x1 = z[:,p*4:2+p*4]
  #  plot(t,x1[:,0],label="x"+str(p+1))
 #   plot(t,x1[:,1],label="y"+str(p+1))
 #   plt.title("Particula "+str(p+1)+'\nPosicion')
 #   plt.legend()

#for p in range(n):
 #   subplot(2,n,p+n+1)
 #   v1 = z[:,2+p*4:4+p*4]
 #   plot(t,v1[:,0],label="vx"+str(p+1))
 #   plot(t,v1[:,1],label="vy"+str(p+1))
 #   plt.title("Velocidad")
 #   plt.legend()

#fig = plt.figure()
#ax = Axes3D(fig)
#for p in range(n):
#	x1 = z[:,p*4:2+p*4]
#	ax.plot(t, x1[:,0],x1[:,1]) 

show()
#------------------

