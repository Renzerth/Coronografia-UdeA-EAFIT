
# coding: utf-8

# ## Simulador
# Simulador del plano del detector en función de la separación angular de dos fuentes puntuales lejanas observadas con un sistema telescópico.
#
# Por: Edgar Rueda, mayo 2018

# In[1]:

from numpy import loadtxt,linspace,meshgrid,sqrt,pi,zeros,sum,ceil,max
from pylab import figure,title,imshow,grid,xticks,yticks,show
from scipy.special import jv
from scipy.misc import imresize
from Functions import zoomC
#get_ipython().magic('matplotlib inline')


# En esta versión se tiene en cuenta el sistema telescópico: Lente f1 + lente f2 + lente f
#
# Basado en el criterio de resolución
# $$\Delta x_c = \theta f = 1.22 \frac{\lambda}{D} f$$
#
# La intensidad del campo se describe como una función proporcional a una Bessel, es decir, se supone corresponde a una onda plana limitada por una pupila y enfocada en el plano focal de la lente.
#
# $$\frac{D^2}{4 R^4}J_1^2\bigg(\frac{\pi D R}{\lambda f}\bigg)$$
#
# El aumento angular $PA$ en un sistema telescópico de dos lentes (objetivo $f_1$ y ocular $f_2$) se define como:
#
# $$PA = \frac{f_1}{f_2}$$
#
# La nueva distancia angular $\alpha_f$ a partir de la distancia angular inicial $\alpha_i$ será:
#
# $$\alpha_f = PA \alpha_i$$
#
# Y su separación en el detector luego de ser enfocados por una lente de distancia focal $f$ y observado con un zoom $Z$ será:
#
# $$\Delta x = Z f\alpha_f$$

# ### Parámetros y separación en el detector

# In[2]:

# parámetros físicos
f1 = 0.2 # focal lente objetivo (m)
f2 = 0.2 # focal lente ocular (m)
f = 0.3 # focal lente formador imagen (m)
D = 0.008 # diámetro apertura lente (m)
wl = 532.0e-9 # longitud de onda (m)
n = 2 #veces el criterio de difracción
ang = n*(1.22)*wl/D # ángulo en términos del criterio de difracción de Rayleigh (rad)

pix = 3.62e-6 # tamaño pixel detector (m)
pixSLM = 26.0e-6 # tamaño pixel modulador (m)
pixN = 500 # muestreo
levels = 256 # resolución detector (bit)
Respix = 2 # muestreo de cada pixel (para mejorar la visualización, no es parámetro físico)
ratio = 1 # razón de intensidades entre el planeta y la estrella

Zoom = 20. # zoom del plano del detector del sistema microscópico

PA = f1/f2 # aumento angular (rad)
ang2 = PA*ang # distancia angular final (rad)

dx = ang2*f # Separación espacial en el plano del dectector (m)

print('PARÁMETROS POR DEFECTO')
print('##########')
print('Focal lente objetivo (m) = ',f1)
print('Focal lente ocular (m) = ',f2)
print('Focal lente formador imagen (m) = ',f)
print('Diámetro apertura lente (m) = ',D)
print('longitud de onda (m) = ',wl)
print('Veces el criterio de difracción = ',n)
print('Tamaño pixel detector (m) = ',pix)
print('Tamaño pixel SLM (m) = ',pixSLM)
print('Muestreo = ',pixN)
print('Número de niveles (bits) = ',levels)
print('Razón intensidad Planeta/estrella = ',ratio)
print('Zoom = ',Zoom)

hago = input('\n¿Desea usar los parámetros por defecto?(s/n)')
if hago == 'n':
    print('Lectura del archivo con parámetros\n')
    print('Recuerde que deben estar en la segunda fila en el orden de aparición de arriba.\n')
    Para = loadtxt('parametros.txt',delimiter=',',skiprows=1)
    f1 = Para[0] # focal lente objetivo (m)
    f2 = Para[1] #focal lente ocular (m)
    f = Para[2] # focal lente formador imagen (m)
    D = Para[3] # diámetro apertura lente (m)
    wl = Para[4] # longitud de onda (m)
    n = Para[5] #veces el criterio de difracción
    ang = n*(1.22)*wl/D # ángulo en términos del criterio de difracción de Rayleigh (rad)

    pix = Para[6] # tamaño pixel detector (m)
    pixSLM = Para[7] # tamaño pixel modulador (m)
    pixN = int(Para[8]) # muestreo
    levels = Para[9] # resolución detector (bit)
    Respix = 2 # muestreo de cada pixel (para mejorar la visualización, no es parámetro físico)
    ratio = Para[10] # razón de intensidades entre el planeta y la estrella

    Zoom = Para[11] # zoom del plano del detector del sistema microscópico
else:
    print('Se usan los parámetros por defecto\n')

print('PARÁMETROS USADOS')
print('##########')
print('Focal lente objetivo (m) = ',f1)
print('Focal lente ocular (m) = ',f2)
print('Focal lente formador imagen (m) = ',f)
print('Diámetro apertura lente (m) = ',D)
print('longitud de onda (m) = ',wl)
print('Veces el criterio de difracción = ',n)
print('Tamaño pixel detector (m) = ',pix)
print('Tamaño pixel SLM (m) = ',pixSLM)
print('Muestreo = ',pixN)
print('Número de niveles (bits) = ',levels)
print('Razón intensidad Planeta/estrella = ',ratio)
print('Zoom = ',Zoom)

print('\nPara una separación angular de %.2e rad (%.2f segundos de arco), o %.f veces el criterio de difracción:'%       (ang,ang*180*3600/(pi),n))
print('Separación de spots = %.2f um (%.f pixeles del detector)' %(Zoom*dx*1e6,Zoom*dx/pix))
#pixN = int(3*dx/pix) # número de pixeles observados; se modifica para que se vean los dos spots



# ### Plano Fase espiral
#
# Se presenta gráficamente como se verían los spots en el plano de la máscara espiral SPP, del coronógrafo.

# In[3]:

dxSLM = ang*f1
pixNSLM = 20 # número de píxeles por lado a observar
RespixSLM = 10 # resolución muestreo por pixel del SLM
N = RespixSLM*pixNSLM # Tamaño matriz simulación spots reales
x = linspace(-pixNSLM*pixSLM,pixNSLM*pixSLM,N) # espacio x-coord
X,Y = meshgrid(x,x) # espacio coordenado 2D
R = sqrt((X+dxSLM/2)**2 + Y**2) # coordenada radial
Star = (0.25*D**2/R**2)*jv(1,pi*D*R/(wl*f1))**2 # Patron estrella, sobre el eje
#Star =zoomC(int(Zoom*len(Star)),len(Star),imresize(Star,Zoom))
#dx = Zoom*dx
R = sqrt((X-dxSLM/2)**2 + Y**2) # coordenada radial desplazada para planeta
planet = (0.25*D**2/R**2)*jv(1,pi*D*R/(wl*f1))**2 # patrón planeta fuera del eje
#planet = zoomC(int(Zoom*len(Star)),len(Star),imresize(planet,Zoom))
PlanoSLM = Star + ratio*planet # Intensidad en el plano detector; planeta con intensidad máxima menor a la estrella


fig = figure(figsize=(12,12))
title(r'Intensidad que llega al SLM; la malla denota los pixeles de tamaño = %.2f um' % (pixSLM*1e6),fontsize=16)
ax = fig.gca()
ax.set_xticks(linspace(0, N, pixNSLM+1))
ax.set_yticks(linspace(0, N, pixNSLM+1))
imshow(PlanoSLM**0.5,cmap='gray',)
grid(color='red',lw=1)
#plt.colorbar()
#plt.show()



# ### Plano Detector Cámara
# Se presenta la imagen observada en el detector luego del aumento programado. La malla corresponde a los pixeles del detector.

# In[4]:

N = Respix*pixN # Tamaño matriz simulación spots reales
x = linspace(-pixN*pix,pixN*pix,N) # espacio x-coord
X,Y = meshgrid(x,x) # espacio coordenado 2D
R = sqrt((X+dx/2)**2 + Y**2) # coordenada radial
Star = (0.25*D**2/R**2)*jv(1,pi*D*R/(wl*f))**2 # Patron estrella, sobre el eje
#Star =zoomC(int(Zoom*len(Star)),len(Star),imresize(Star,Zoom))
#dx = Zoom*dx
R = sqrt((X-dx/2)**2 + Y**2) # coordenada radial desplazada para planeta
planet = (0.25*D**2/R**2)*jv(1,pi*D*R/(wl*f))**2 # patrón planeta fuera del eje
#planet = zoomC(int(Zoom*len(Star)),len(Star),imresize(planet,Zoom))
Plano = Star + ratio*planet # Intensidad en el plano detector; planeta con intensidad máxima menor a la estrella
Plano = zoomC(int(Zoom*len(Star)),len(Star),imresize(Plano,Zoom))

# Simulación de visualización con detector segun resolución pixeles y resolución rango dinámico
Detector = zeros((pixN,pixN))
for ii in range(pixN):
    for jj in range(pixN):
        Aux = sum(Plano[ii*Respix:(ii+1)*Respix,jj*Respix:(jj+1)*Respix])
        Detector[ii,jj] = Aux
DetBit = ceil(levels*Detector/max(Detector))


fig = figure(figsize=(16,16))
title('Intensidad que llega al detector; la malla denota los pixeles de tamaño = %.2f um' % (pix*1e6), fontsize=16)
ax = fig.gca()
ax.set_xticks(linspace(0, N, pixN+1))
ax.set_yticks(linspace(0, N, pixN+1))
imshow(Plano**0.5,cmap='gray',)
grid(color='red',lw=0.5)
#plt.colorbar()
#plt.show()



figure(figsize=(16,16))
title('Simulación de como se observaría en el detector para %g bits de profundidad'% (levels),fontsize=16)
xticks(linspace(0, pixN, pixN+1))
yticks(linspace(0, pixN, pixN+1))
imshow(DetBit,cmap='gray',interpolation='nearest')
#plt.colorbar()
show()


# In[ ]:
