from math import sin, cos, tan, asin, acos, atan, radians, degrees, sqrt, pi
import numpy as np
import matplotlib.pyplot as plt

#Constantes
mb = 0.293 #Masse baton (kg)
mp = 0.20275579 #Masse poids (kg)
g = 9.81 #Accélération gravitationnelle (m/s^2)
hi = 1 #Hauteur initiale de la masse (m)
Ib = 0.042033 #Moment d'inertie du baton (kg*m^2)
r = 0.0508 #Rayon de la poulie (m)
eff = 0.5 #Efficacité du système
Lb = 0.63 #Longueur du baton (m)
Lc = 0.291 #Longueur du centre de masse depuis le point de pivot (m)
Ang = pi/4 #Angle du baton (rad)
Gr = 1/10 #Rapport de réduction

print("CALCULS AVANT BILANS")
#Torque nécessaire pour lever baton
T = mb*g*Lc*sin(Ang)
print("Torque nécessaire pour lever le baton: ", T, "N*m")

#Torque émise par le système
Te = eff*mp*g*r/Gr
print("Torque émise par le système: ", Te, "N*m")

#Hauteur finale du baton
nbtourpoulie = hi/(2*pi*r)
nbtourbaton = nbtourpoulie*Gr
hbf = nbtourbaton*Lc*cos(Ang)

print("Nombre de tours de la poulie: ", nbtourpoulie)
print("Nombre de tours du baton: ", nbtourbaton)
if nbtourbaton > 0.5:
    print("Le baton a trop tourné")

if hbf < 0:
    hbf = hbf*-1 + Lc*cos(Ang)

#hbf = hbf*sin(Ang)


#Bilan 1
print("")
print("BILAN 1")
Ppi = mp*g*hi #Énergie potentielle initiale du poids (J)
Vp = sqrt((2*g*(mp*hi*eff-mb*hbf))/(mp+(Ib*Gr**2)/r**2)) #Vitesse finale du poids (m/s)
print("Vitesse finale du poids: ", Vp, "m/s")
print("Hauteur finale du baton : ", hbf, "m")

#Bilan 2
print("")
print("BILAN 2")
wb = Vp*Gr/r #Vitesse angulaire du baton (rad/s)
print("Vitesse angulaire du baton: ", wb, "rad/s")
hbf2 = ((Ib*(wb**2))/2 + mb*g*hbf)/(mb*g*cos(Ang)) #Hauteur finale du baton (m)
print("Hauteur finale du baton: ", hbf2, "m")

#Bilan 3
print("")
print("BILAN 3")
wbf = sqrt(2*mb*g*hbf2/Ib) #Vitesse angulaire du baton (rad/s)
print("Vitesse angulaire finale du baton: ", wbf, "rad/s")

#Pour simulation
Gr = np.linspace(1/20, 1/2, 100)

nbtourpoulie = hi/(2*pi*r)
nbtourbaton = nbtourpoulie*Gr
hbf = nbtourbaton*Lc*cos(Ang)

index = np.where(hbf < 0)[0]
if len(index) > 0:
    hbf[index] = hbf[index]*-1 + Lc*cos(Ang)

Ppi = mp*g*hi #Énergie potentielle initiale du poids (J)
Vp = np.sqrt((2*g*(mp*hi*eff-mb*hbf))/(mp+(Ib*Gr**2)/r**2)) #Vitesse finale du poids (m/s)
#Ti = Ppi*eff*r/(Vp*Gr)
Ti = eff*mp*g*r/Gr

wb = Vp*Gr/r #Vitesse angulaire du baton (rad/s)
hbf2 = ((Ib*wb**2)/2 + mb*g*hbf)/(mb*g*cos(Ang)) #Hauteur finale du baton (m)

wbf = np.sqrt(2*mb*g*hbf2/Ib) #Vitesse angulaire du baton (rad/s)

#Graphique
plt.figure(1)
plt.suptitle("Bilan 1")

plt.subplot(131)
plt.title("hauteur b")
plt.plot(Gr, hbf)

plt.subplot(132)
plt.title("vitesse p")
plt.plot(Gr, Vp)

plt.subplot(133)
plt.title("torque")
plt.plot(Gr, Ti)
plt.plot(Gr, T*np.ones(len(Gr)))

plt.figure(2)
plt.suptitle("Bilan 2")
plt.subplot(121)
plt.title("hauteur b")
plt.plot(Gr, hbf2)

plt.subplot(122)
plt.title("vitesse b")
plt.plot(Gr, wb)

plt.figure(3)
plt.suptitle("Bilan 3")
plt.title("vitesse b")
plt.plot(Gr, wbf)

plt.show()