#!/usr/bin/env python
# coding: utf-8

# In[15]:


import math


# # GLOBAL Equations 

# In[16]:


def Rs_std(p,pb,t,so,sg):
    if p>pb:
        p=pb
    API = 141.5/so-131.5
    yg=0.00091*t-0.0125*API
    Rs=sg*(p/(18*10**yg))**1.204
    return Rs
def Rs_Vqz(p,pb,t,so,sg):
    P_sep=590
    T_sep=65
    if p>pb:
        p=pb
    API = 141.5/so-131.5
    sg_c = sg*(1+(5.912*(10**-5)*API*T_sep*math.log(P_sep/114.7,10)))
    if API <=30:
        a=1.0937
        c=0.0362*sg_c*math.exp(25.7240*API/(t+459.67))
    else:
        a=1.187
        c=0.0178*sg_c*math.exp(23.931*API/(t+459.67))
    Rs=c*p**a
    return Rs

def Rs_PF(p,pb,t,so,sg):
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        x=(7.916*10**(-4)*API**(1.5410))-(4.561*10**(-5)*t**(1.3911)) # t in F
        Rs=((p/112.727+12.340)*sg**(0.8439)*10**x)**1.73184
        return Rs
    
def Rs_Krt(p,pb,t,so,sg):
    P_sep=590
    T_sep=65
    if p>pb:
        p=pb
    API = 141.5/so-131.5
    sg_c=sg*(1+0.159*(API**0.4078)*(T_sep**(-0.2466))*math.log(P_sep/114.7,10))
    if API<=30:
        x=13.1405*API/(t+459.67) #T in F
        c=10**x
        Rs=0.05958*c*(sg_c**0.7972)*p**1.0014
    else:
        x=11.2895*API/(t+459.67)
        c=10**x
        Rs=0.0315*c*(sg_c**0.7587)*p**1.0937
    return Rs   

def co_Vzq(p,pb,t,sg,so):
    API=141.5/so-131.5
    P_sep=590
    T_sep=65
    Rs=Rs_Vqz(pb,pb,t,so,sg)
    if p<=pb:
        co=math.exp(-7.573-1.450*math.log(p)-0.383*math.log(pb)+1.402*math.log(t+459.67)+0.256*math.log(API)+0.449*math.log(Rs))
        return co
    else:
        co=(-1433+5*Rs+17.2*t-1180*sg+12.61*API)/(p*10**5)
        return co
def z_drk(p,t,sg):
    Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
    Tpc = 169.2 + 349.5*sg - 74*(sg**2)
    tpr = (t+459.67)/Tpc
    ppr = p/Ppc
    a = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
    rho = 0.27*ppr/tpr
    zfactor = (a[0] + a[1]/tpr + a[2]/(tpr**3) + a[3]/(tpr**4) + a[4]/(tpr**5))*rho + (a[5] + a[6]/tpr + a[7]/(tpr**2))*(rho**2) - a[8]*(a[6]/tpr + a[7]/(tpr**2))*(rho**5) + a[9]*(1+a[10]*(rho**2))*((rho**2)/(tpr**3))*math.exp(-1*a[10]*(rho**2)) + 1
    error = 999
    while (error > 0.0001):
        rho = 0.27*ppr/(zfactor*tpr)
        newzfactor = (a[0] + a[1]/tpr + a[2]/(tpr**3) + a[3]/(tpr**4) + a[4]/(tpr**5))*rho + (a[5] + a[6]/tpr + a[7]/(tpr**2))*(rho**2) - a[8]*(a[6]/tpr + a[7]/(tpr**2))*(rho**5) + a[9]*(1+a[10]*(rho**2))*((rho**2)/(tpr**3))*math.exp(-1*a[10]*(rho**2)) + 1
        error = abs((newzfactor - zfactor)/zfactor)
        zfactor = newzfactor
    return zfactor
def z_hall(p,t,sg):
    Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
    Tpc = 169.2 + 349.5*sg - 74*(sg**2)
    tpr= (t+459.67)/Tpc
    ppr=p/Ppc
    ttpr = 1/tpr
    aval = 0.06125*ttpr*math.exp(-1.2*((1 - ttpr)**2))
    bval = ttpr*(14.76 - 9.76*ttpr + 4.58*(ttpr**2))
    cval = ttpr*(90.7 - 242.2*ttpr + 42.4*(ttpr**2))
    dval = 2.18 + 2.82*ttpr
    yval = 0.061
    fval = -1*aval*ppr + (yval + yval**2 + yval**3 - yval**4)/((1 - yval)**3) - bval*(yval**2) + cval*(yval**dval)
    error = 999
    while (abs(fval) > 0.00000001 and error > 0.00000001):
        fval = -1*aval*ppr + (yval + yval**2 + yval**3 - yval**4)/((1 - yval)**3) - bval*(yval**2) + cval*(yval**dval)
        dfval = (1 + 4*yval + 4*(yval**2) - 4*(yval**3) + yval**4)/((1 - yval)**4) - 2*bval*yval + cval*dval*(yval**(dval - 1))
        nyval = yval - fval/dfval
        error = abs(nyval - yval)
        yval = nyval
    zfactor = aval*ppr/yval
    return zfactor
def Bw_mcn(p,pb,t):
    cw=3*10**-6
    if p<=pb:
        vt=-1.0001*10**-2+t*1.33391*10**-4+(t**2)*5.50654*10**-7
        vp=-1.95301*p*t*10**-9-1.72834*(p**2)*t*10**-13-3.58922*p*10**-7-2.25341*(p**2)*10**-10
        Bw=(1+vt)*(1+vp)
        return Bw
    else:
        vt=-1.0001*10**-2+t*1.33391*10**-4+(t**2)*5.50654*10**-7
        vp=-1.95301*pb*t*10**-9-1.72834*(pb**2)*t*10**-13-3.58922*pb*10**-7-2.25341*(pb**2)*10**-10
        Bwp=(1+vt)*(1+vp)
        Bw=Bwp*math.exp(-cw*(p-pb))
        return Bw


# # Density

# In[17]:


#Standing
def rho_Std(p,pb,t,sg,so):# t in f
    co_pb = co_Vzq(pb,pb,t,sg,so) 
    API = 141.5/so - 131.5
    Rs=Rs_std(p,pb,t,so,sg)
    if p<=pb:
        rho_a=38.52/(10**(0.00326*API))+(94.75-33.93*math.log(API,10))*math.log(sg,10)
        rho_po = (Rs*sg+4600*sg)/(73.71+Rs*sg/rho_a)
        drho_p = (0.167+16.181/(10**(0.0425*rho_po)))*(p/1000)-0.01*(0.299+263/(10**(0.0603*rho_po)))*((p/1000)**2)
        rho_p=rho_po+drho_p
        rho_t=(0.0133+152.4/(rho_p**2.45))*(t-60)-(0.0000081-0.0622/(10**(0.0764*rho_p)))*((t-60)**2)
        rho = rho_p-rho_t
        return rho
    else: 
        rho_a=38.52/(10**(0.00326*API))+(94.75-33.93*math.log(API,10))*math.log(sg,10)
        rho_po = (Rs*sg+4600*sg)/(73.71+Rs*sg/rho_a)
        drho_p = (0.167+16.181/(10**(0.0425*rho_po)))*(pb/1000)-0.01*(0.299+263/(10**(0.0603*rho_po)))*((pb/1000)**2)
        rho_p=rho_po+drho_p
        rho_t=(0.0133+152.4/(rho_p**2.45))*(t-60)-(0.0000081-0.0622/(10**(0.0764*rho_p)))*((t-60)**2)
        rho_pb = rho_p-rho_t
        rho=rho_pb*math.exp(co_pb*(p-pb)) #co=spivey,valko,mc cain
        return rho


# # Bo

# In[18]:


def BO(p,pb,t,so,sg,method):
    if method=='Standing':
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_std(p,pb,t,so,sg)
        Cn = Rs*((sg/so)**0.5)+1.25*t
        Bo = 0.9759+0.00012*(Cn**1.2)
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Vasquez&Begg':
        #Vasquez
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_Vqz(p,pb,t,so,sg)
        P_sep=590
        T_sep=65
        API=141.5/so - 131.5
        sg_c = sg*(1+(5.912*(10**-5)*API*T_sep*math.log(P_sep/114.7,10)))
        if API <= 30:
            Bo = 1+4.677*10**(-4)*Rs+1.751*10**(-5)*(API/sg_c)*(t+459.67-60)-1.8106*10**(-8)*(Rs*(API/sg_c)*(t-60))
        else:
            Bo = 1+4.677*10**(-4)*Rs+1.1*10**(-5)*(API/sg_c)*(t+459.67-60)+1.337*10**(-9)*(Rs*(API/sg_c)*(t-60))
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Petrosky': 
        #Petrosky
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_PF(p,pb,t,so,sg)
        Bo=1.0113+7.2046*10**(-5)*((Rs**0.3738)*(sg**0.2194)/(so**0.6265)+0.24626*t**(0.5371))**3.0936
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Farshad':
        #Farshad
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_PF(p,pb,t,so,sg)
        G=(Rs**0.5956)*(sg**0.2369)*(so**-1.3282)+0.0976*t
        Bo=1+10**(-2.6541+0.5576*math.log(G,10)+0.3331*math.log(G,10)**2)
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Marhoun':
        #Marhoun
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_std(p,pb,t,so,sg)
        Bo=1+0.177342*(10**-3)*Rs+0.220163*(10**-3)*Rs*(sg/so)+4.292580*(10**-6)*Rs*(1-so)*(t+459.67-60)+0.528707*(10**-3)*(t-60)
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Kartoatmodjo':
        #Kartoatmodjo
        
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_Krt(p,pb,t,so,sg)
        P_sep=590
        T_sep=65
        API = 141.5/so -131.5
        sg_c=sg*(1+0.159*(API**0.4078)*(T_sep**-0.2446)*math.log(P_sep/114.7,10))
        Bo=0.98496+0.0001*((Rs**0.755)*(sg_c**0.25)/(so**1.5)+0.45*t)**1.5
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Almedhaideb':
        #Almedhaideb
        
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_std(p,pb,t,so,sg)
        Bo=1.122018+1.41*(10**-6)*Rs*t/(so**2)
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Shammasi':
        #Al-Shammasi
        
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_std(p,pb,t,so,sg)
        Bo=1+5.53*(10**-7)*(Rs*(t-60))+0.000181*(Rs/so)+0.000449*((t-60)/so)+0.000206*(Rs*sg/so)
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo
    elif method=='Sharkawy':
        #El-Sharkawy
        
        co = co_Vzq(pb,pb,t,sg,so)
        Rs=Rs_std(p,pb,t,so,sg)
        Bo=1+40.428*(10**-5)*Rs+63.802*(10**-5)*(t-60)+0.78*(10**-6)*Rs*(t-60)*sg/so
        if p<=pb:
            return Bo
        else:
            Bo = Bo * math.exp(-co*(p-pb))
            return Bo


# # GOR

# In[19]:


def GOR(p,pb,t,so,sg,method):
    if method=='Standing':
        #Standing
        
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        yg=0.00091*t-0.0125*API
        Rs=sg*(p/(18*10**yg))**1.204
        return Rs
    elif method=='Vasquez&Begg':
        #Vasquez
        
        P_sep=590
        T_sep=65
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        sg_c = sg*(1+(5.912*(10**-5)*API*T_sep*math.log(P_sep/114.7,10)))
        if API <=30:
            a=1.0937
            c=0.0362*sg_c*math.exp(25.7240*API/(t+459.67))
        else:
            a=1.187
            c=0.0178*sg_c*math.exp(23.931*API/(t+459.67))
        Rs=c*p**a
        return Rs
    elif method=='Lasater':
        #Lasater
        
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        a=p*sg/(t+459.67) 
        if API<=40:
            M=630-10*API
        else:
            M=73110*(API)**(-1.562)
        if a<3.29:
            y=0.359*math.log((1.473*p*sg)/(t+459.67)+0.476)
        else:
            y=((0.121*p*sg/(t+459.67))-0.236)**(-0.281)
        Rs=132755*so*y/M/(1-y)
        return Rs
    elif method=='Petrosky&Farshad':
        #Petrosky&Farshad
        
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        x=(7.916*10**(-4)*API**(1.5410))-(4.561*10**(-5)*t**(1.3911)) # t in F
        Rs=((p/112.727+12.340)*sg**(0.8439)*10**x)**1.73184
        return Rs
    elif method=='Kartoatmodjo':
        #Kartoatmodjo-Schmidt
        
        P_sep=590
        T_sep=65
        if p>pb:
            p=pb
        API = 141.5/so-131.5
        sg_c=sg*(1+0.159*(API**0.4078)*(T_sep**(-0.2466))*math.log(P_sep/114.7,10))
        if API<=30:
            x=13.1405*API/(t+459.67) #T in F
            c=10**x
            Rs=0.05958*c*(sg_c**0.7972)*p**1.0014
        else:
            x=11.2895*API/(t+459.67)
            c=10**x
            Rs=0.0315*c*(sg_c**0.7587)*p**1.0937
        return Rs   


# # Oil Viscosity

# In[21]:


def MIUO(p,pb,t,so,sg,method):
    if method=='Khan':
        #Khan
        
        Rs=Rs_std(p,pb,t,so,sg)    
        a=0.09*(sg**0.5)/(Rs**(1/3)*((t+459.67)/459.67)**4.5)/(1-so)**3
        if p<=pb:
            visc=a*((p/pb)**-0.14)*math.exp((-2.5*10**-4)*(p-pb))
        else:
            visc=a*math.exp((9.6*10**-5)*(p-pb))
        return visc
    elif method=='Begss':
        #Begs
        
        Rs=Rs_std(p,pb,t,so,sg) 
        API = 141.5/so-131.5
        x=(t**-1.168)*math.exp(6.9824-0.04658*API)
        a=(10**x)-1
        A=10.715*(Rs+100)**-0.515
        B=5.44*(Rs+150)**-0.338
        visc1=A*a**B
        if p<=pb:
            return visc1
        else:
            m=2.6*(p**1.187)*math.exp(-11.513-8.98*(10**-5)*p)
            visc = visc1*(p/pb)**m
            return visc
    elif method=='Beal':
        #Beal
        Rs=Rs_std(p,pb,t,so,sg) 
        API = 141.5/so-131.5
        x=(t**-1.168)*math.exp(6.9824-0.04658*API)
        a=(10**x)-1
        A=10.715*(Rs+100)**-0.515
        B=5.44*(Rs+150)**-0.338
        visc1=A*a**B
        if p<=pb:
            return visc1
        else:
            visc=visc1+0.001*(0.024*visc1**1.6+0.038*visc1**0.56)*(p-pb)
            return visc    


# # Oil Compressibility

# In[22]:


def co_Vzq(p,pb,t,sg,so):
    API=141.5/so-131.5
    P_sep=590
    T_sep=65
    Rs=Rs_Vqz(pb,pb,t,so,sg)
    if p<=pb:
        co=math.exp(-7.573-1.450*math.log(p)-0.383*math.log(pb)+1.402*math.log(t+459.67)+0.256*math.log(API)+0.449*math.log(Rs))
        return co
    else:
        co=(-1433+5*Rs+17.2*t-1180*sg+12.61*API)/(p*10**5)
        return co


# # Z factor

# In[23]:


def ZF(p,t,sg,method):
    if method=='Papay':  
        #Papay
        Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
        Tpc = 169.2 + 349.5*sg - 74*(sg**2)
        Tpr = (t+459.67)/Tpc
        Ppr = p/Ppc
        a = 1.39*(Tpr-0.92)**0.5 - 0.36*Tpr - 0.1
        c = 0.132 - 0.32*math.log10(Tpr)
        e = 9*(Tpr-1)
        f = 0.3106 - 0.49*Tpr + 0.1824*Tpr**2
        d = 10**f
        b = (0.62-0.23*Tpr)*Ppr + ((0.066/(Tpr-0.86)) - 0.037)*Ppr**2 + (0.32*Ppr**6/10**e)
        z = a + ((1-a)/math.exp(b)) + c*Ppr**d
        return z
    elif method=='Dranchuk Abou Kassem':
        #Dranchuk
        Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
        Tpc = 169.2 + 349.5*sg - 74*(sg**2)
        tpr = (t+459.67)/Tpc
        ppr = p/Ppc
        a = [0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134, 0.7210]
        rho = 0.27*ppr/tpr
        zfactor = (a[0] + a[1]/tpr + a[2]/(tpr**3) + a[3]/(tpr**4) + a[4]/(tpr**5))*rho + (a[5] + a[6]/tpr + a[7]/(tpr**2))*(rho**2) - a[8]*(a[6]/tpr + a[7]/(tpr**2))*(rho**5) + a[9]*(1+a[10]*(rho**2))*((rho**2)/(tpr**3))*math.exp(-1*a[10]*(rho**2)) + 1
        error = 999
        while (error > 0.0001):
            rho = 0.27*ppr/(zfactor*tpr)
            newzfactor = (a[0] + a[1]/tpr + a[2]/(tpr**3) + a[3]/(tpr**4) + a[4]/(tpr**5))*rho + (a[5] + a[6]/tpr + a[7]/(tpr**2))*(rho**2) - a[8]*(a[6]/tpr + a[7]/(tpr**2))*(rho**5) + a[9]*(1+a[10]*(rho**2))*((rho**2)/(tpr**3))*math.exp(-1*a[10]*(rho**2)) + 1
            error = abs((newzfactor - zfactor)/zfactor)
            zfactor = newzfactor
        return zfactor
    elif method=='Hall Yaborough':
        #HallYaborough
        Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
        Tpc = 169.2 + 349.5*sg - 74*(sg**2)
        tpr= (t+459.67)/Tpc
        ppr=p/Ppc
        ttpr = 1/tpr
        aval = 0.06125*ttpr*math.exp(-1.2*((1 - ttpr)**2))
        bval = ttpr*(14.76 - 9.76*ttpr + 4.58*(ttpr**2))
        cval = ttpr*(90.7 - 242.2*ttpr + 42.4*(ttpr**2))
        dval = 2.18 + 2.82*ttpr
        yval = 0.061
        fval = -1*aval*ppr + (yval + yval**2 + yval**3 - yval**4)/((1 - yval)**3) - bval*(yval**2) + cval*(yval**dval)
        error = 999
        while (abs(fval) > 0.00000001 and error > 0.00000001):
            fval = -1*aval*ppr + (yval + yval**2 + yval**3 - yval**4)/((1 - yval)**3) - bval*(yval**2) + cval*(yval**dval)
            dfval = (1 + 4*yval + 4*(yval**2) - 4*(yval**3) + yval**4)/((1 - yval)**4) - 2*bval*yval + cval*dval*(yval**(dval - 1))
            nyval = yval - fval/dfval
            error = abs(nyval - yval)
            yval = nyval
        zfactor = aval*ppr/yval
        return zfactor


# # Gas Density

# In[24]:


def rho_g(p,t,sg):
    M=sg/28.9625
    R=10.732
    z=z_drk(p,t,sg)
    rho=p*M/z/R/(t+459.67)
    return rho


# # Gas Formation Volume Factor

# In[25]:


def Bg(p,t,sg):
    z=z_drk(p,t,sg)
    Bg=(14.65*(t+459.67)*z)/(519.67*p)
    return Bg


# # Gas Viscosity

# In[26]:


#Lee
def miug_lee(p,t,sg):
    rhog=rho_g(p,t,sg)
    M=sg*28.9625
    x=3.448+986.4/(t+459.67)+0.01009*M
    y=2.447-0.2224*x
    K=(9.379+0.01607*M)*((t+459.67)**1.5)/(209.2+19.26*M+(t+459.67))
    miug = K*math.exp(x*rhog**y)/10**4
    return miug


# # Gas Compressibility

# In[27]:


def CG(p,t,sg,method):
    if method=='Hall Yaborough':
        #Hall Yaborough
        z=z_hall(p,t,sg)
        Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
        Tpc = 169.2 + 349.5*sg - 74*(sg**2)
        Tpr = (t+459.67)/Tpc
        Ppr = p/Ppc
        t=1/Tpr
        a=0.06125*t*math.exp(-1.2*(1-t)**2)
        y=a*Ppr/z
        f_y=(1+4*y+4*y**2-4*y**3+y**4)/((1-y)**4)-(29.52*t-19.52*t**2+9.16*t**3)*y+(2.18+2.82*t)*(90.7*t-242.2*t**2+42.4*t**3)*y**(1.18+2.82*t)
        y_p=1*a/f_y
        z_p=a/y-a*Ppr/(y**2)*y_p
        cpr=1/Ppr-1/z*z_p
        cg=cpr/Ppc
        return cg
    elif method=='Dranchuk Abou Kassem':
        #Dranchuk
        z=z_drk(p,t,sg)
        Ppc = 756.8 - 131.07*sg - 3.6*(sg**2)
        Tpc = 169.2 + 349.5*sg - 74*(sg**2)
        Tpr = (t+459.67)/Tpc
        Ppr = p/Ppc
        a=0.3265
        b=-1.0700
        c=-0.5339
        d=0.01569
        e=-0.05165
        f=0.5475
        g=-0.7361
        h=0.1844
        i=0.1056
        j=0.6134
        k=0.7210
        rho=0.27*Ppr/z/Tpr
        z_p=(a+b/Tpr+c/Tpr**3+d/Tpr**4+e/Tpr**5)+2*(f+g/Tpr+h/Tpr**2)*rho-5*i*(g/Tpr+h/Tpr**2)*rho**4+2*j*(1+k*rho**2-(k**2)*(rho**4))*rho/(Tpr**3)*math.exp(-1*k*rho**2)
        x=rho/z*z_p
        cpr=(1-x/(1+x))/Ppr
        cg=cpr/Ppc
        return cg


# # FVF Water

# In[28]:


def BW(p,pb,t,nacl,method):
    if method=='Meehan':   
        #Meehan
        nacl=nacl/1000000
        if p<=pb:
            a=0.9947+t*5.8*10**-6+(t**2)*1.02*10**-6
            b=-4.228*10**-6+t*1.8376*10**-8-(t**2)*6.77*10**-11
            c=1.3*10**-10-t*1.3855*10**-12+(t**2)*4.285*10**-15
        else:
            a=0.9911+t*6.35*10**-6+(t**2)*8.5*10**-7
            b=-1.093*10**-6+t*3.497*10**-9-(t**2)*4.57*10**-12
            c=-5*10**-11-t*6.429*10**-13+(t**2)*1.43*10**-15
        S=1+nacl*((5.1*10**-8)*p+(5.47*10**-6-p*1.96*10**-10)*(t-60)+(-3.23*10**-8+p*8.5*10**-12)*(t-60)**2)
        Bw=(a+b*p+c*p**2)*S
        return Bw
    elif method=='Mc Cain':
        #McCain
        cw=3*10**-6
        if p<=pb:
            vt=-1.0001*10**-2+t*1.33391*10**-4+(t**2)*5.50654*10**-7
            vp=-1.95301*p*t*10**-9-1.72834*(p**2)*t*10**-13-3.58922*p*10**-7-2.25341*(p**2)*10**-10
            Bw=(1+vt)*(1+vp)
            return Bw
        else:
            vt=-1.0001*10**-2+t*1.33391*10**-4+(t**2)*5.50654*10**-7
            vp=-1.95301*pb*t*10**-9-1.72834*(pb**2)*t*10**-13-3.58922*pb*10**-7-2.25341*(pb**2)*10**-10
            Bwp=(1+vt)*(1+vp)
            Bw=Bwp*math.exp(-cw*(p-pb))
            return Bw


# # Water Density

# In[29]:


#McCain
def rhow_mcn(p,pb,t,nacl):
    nacl=nacl/1000000
    Bw=Bw_mcn(p,pb,t)
    rho=62.368+0.438603*nacl+1.60074*(nacl**2)*10**-3
    rhow=rho/Bw
    return rhow


# # Water Viscosity

# In[31]:


def MIUW(p,t,nacl,method):
    if method=='Meehan':
        #Meehan
        nacl=nacl/1000000
        f=1+3.5*(10**-12)*(t-40)*p**2
        a=-0.04518+0.009313*nacl-0.000393*nacl**2
        b=70.634+0.09576*nacl**2
        miu=a+b/t
        miuw=miu*f
        return miuw
    elif method=='Mc Cain':
        #McCain
        nacl=nacl/1000000
        a=1.09574*10**2-8.40564*nacl+(nacl**2)*3.13314*10**-1+(nacl**3)*8.72213*10**-3
        b=-1.12166+nacl*2.63951*10**-2-(nacl**2)*6.79461*10**-4-(nacl**3)*5.47119*10**-5+(nacl**4)*1.55586*10**-6
        miu=a*t**b
        miuw=miu*(0.9994+4.0295*p*10**-5+3.1062*(p**2)*10**-9)
        return miuw


# In[ ]:




