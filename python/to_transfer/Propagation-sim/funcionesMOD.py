"""Some useful functions

function(InfLim,SupLim,step,Total_Width,M)
M boolean --> True: Returns spatial coords. and function value, 
              False: Returns only function value
              """
import numpy as np

def rect(a,b,dx,w,M=False):    
    s = []
    X = []
    for x in np.arange(a,b,dx):
        X.append(x)
        if abs(x/w) < 0.5: s.append(1.)
        #elif abs(x/w) == 0.5: s.append(1./2.)
        else: s.append(0.) 
    
    if M == True: return X,s
    if M == False: return s

def triangle(a,b,dx,w,M=False):
    s = []
    X = []
    for x in np.arange(a,b,dx):
        X.append(x)
        if abs(x/w) <= 1.: s.append(1.-abs(x/w))
        else: s.append(0.)
    
    if M == True: return X,s
    if M == False: return s


def gaussian(a,b,dx,w,M=False):
    s = []
    X = []
    for x in np.arange(a,b,dx):
        X.append(x)
        s.append(np.exp(-np.pi*(x/w)*(x/w)))
    if M == True: return X,s
    if M == False: return s

def circ(a,b,dx,w,M=False):    
    s,X,Y = [],[],[]
    L = abs((b-a)/dx)  #shape de la matriz     
    #if int(L)==L: L=int(L)
    #if int(L)-L < 0.0 : 
     #   L = L+1  #Para corregir el error en decimales y calcular bien el int(L) Ej.: L=40.99999997 (41.0), int(L)=40 
    #else: L=int(L)
    L = int(round(L)) #round() approximates to the nearest
    i=0
    j=0
    for x in np.arange(a,b,dx):
        if i==L: break #This is to prevent that the number of samples 
                       # exceed its real number, this can happen
                       #  due to the little numbers added by 
                       #   natural wrongs of computer
        X.append(x)
        for y in np.arange(a,b,dx):
            
            if j==L: break
            if x==a: Y.append(y) #Para que no haga L veces el append

            if np.sqrt(x*x+y*y)/abs(w) < 0.5: s.append(1.)

         #   elif np.sqrt(x*x+y*y)/abs(w) == 0.5: s.append(0.5)
            
            else: s.append(0.)
            j+=1
            
        i+=1
        j=0
    
    s = np.array(s).reshape(L,L)  #Vuelve el array s una matriz cuadrada
    
    if M == True: return X,Y,s
    if M == False: return s

    """
    #Another option to do it...  

    #s = np.zeros((L,L)) #n    
                                                                                                                                                     
    #i=0;j=0 #n                                                                                                                                                                     

    #x=np.arange(a,b,dx)                                                                                                                                                            
    #y=x                                                                                                                                                                            
    #for j in np.arange(0,L):                                                                                                                                                       
    for x in np.arange(a,b,dx):

        X.append(x[j])

        #for i in np.arange(0,L):                                               
                                                                                                    
        for y in np.arange(a,b,dx):

            if x==a: Y.append(y[i]) #Para que no haga L veces el append                                                                                                             

            if np.sqrt(x[j]*x[j]+y[i]*y[i])/abs(w) < 0.5: s.append(1.)#s[i,j]=1.                                                                                                    

         #   elif np.sqrt(x*x+y*y)/abs(w) == 0.5: s.append(0.5)                                                                                                                     

            else: s[i,j]=0.#s.append(0.)                                                                                                                                            
            #print i,j,a,b,dx,L,y                                                                                                                                                   
            #i+=1 #i para y: filas                                                                                                                                                 

        #j+=1  #j para x: columnas                                                                                                                                                  
        #i=0                                                                                                                                                                        

    #s = np.array(s).reshape(L,L)  #Vuelve el array s una matriz cuadrada                                                                                                           
    if M == True: return X,Y,s
    if M == False: return s
    
    """
