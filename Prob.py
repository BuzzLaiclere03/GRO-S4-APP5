from math import sin, cos, tan, asin, acos, atan, radians, degrees, sqrt, pi

#Constantes
mb = 0.293 #Masse baton (kg)
mp = 0.20275579 #Masse poids (kg)
g = 9.81 #Accélération gravitationnelle (m/s^2)
hi = 0.1 #Hauteur initiale (m)
Ib = 0.042033 #Moment d'inertie du baton (kg*m^2)
r = 0.0508 #Rayon de la poulie (m)
eff = 0.75 #Efficacité du système
Lb = 0.63 #Longueur du baton (m)
Lc = 0.291 #Longueur du centre de masse depuis le point de pivot (m)
Ang = pi/4 #Angle du baton (rad)
Gr = 1.5 #Rapport de réduction

print("CALCULS AVANT BILANS")
#Torque nécessaire pour lever baton
T = mb*g*Lc*sin(Ang)
print("Torque nécessaire pour lever le baton: ", T, "N*m")

#Torque émise par le système
Te = eff*mp*g*hi/Gr
print("Torque émise par le système: ", Te, "N*m")

#Hauteur finale du baton
nbtourpoulie = hi/(2*pi*r)
nbtourbaton = nbtourpoulie*Gr
hbf = nbtourbaton*Lb*cos(Ang)

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
Vp = sqrt((2*g*(mp*hi*eff+mb*hbf))/(mp*Ib*(1+(Gr*r)**2))) #Vitesse finale du poids (m/s)
print("Vitesse finale du poids: ", Vp, "m/s")
print("Hauteur finale du baton : ", hbf, "m")

#Bilan 2
print("")
print("BILAN 2")
wb = Vp*Gr*r #Vitesse angulaire du baton (rad/s)
print("Vitesse angulaire du baton: ", wb, "rad/s")
hbf2 = ((Ib*wb**2)/2 + mb*g*hbf)/(mb*g) #Hauteur finale du baton (m)
print("Hauteur finale du baton: ", hbf2, "m")

#Bilan 3
print("")
print("BILAN 3")
wbf = sqrt(2*mb*g*hbf2/Ib) #Vitesse angulaire du baton (rad/s)
print("Vitesse angulaire finale du baton: ", wbf, "rad/s")
