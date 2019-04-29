import numpy as np
import matplotlib.pyplot as plt

#intial variables
massList = np.array([])
numPosMass= 50 

numNegMass= 0
radius = 5000
numSim = 100
G = 1.0

#sets the masses
for i in range(0,numPosMass+numNegMass):
    massList = np.append(massList,0.01)

#Puts the objects in random places
xPos = np.array(np.random.uniform(-radius,radius,len(massList)))
yPos = np.array(np.random.uniform(-radius,radius,len(massList)))
zPos = np.array(np.random.uniform(-radius,radius,len(massList)))

#For each object, finds two different random angles for velocity direction
phi = np.random.uniform(0, 2*np.pi, numPosMass)
theta = np.random.uniform(0,(np.pi)/2, numPosMass)
#finds what the circlar velocity should be given position
cirVel = np.sqrt(G*massList[i]/np.sqrt(xPos**2+yPos**2+zPos**2))
#finds the cartesian coponents of velocties
xVel = cirVel*np.sin(theta)*np.cos(phi)
yVel = cirVel*np.sin(theta)*np.sin(phi)
zVel = cirVel*np.cos(theta)

def updateVelocities(xVel,yVel,zVel):
    a_x = np.array([])
    a_y = np.array([])
    a_z = np.array([])

    for i in range(0,len(massList)):       
        temp_ax_array = np.array([])
        temp_ay_array = np.array([])
        temp_az_array = np.array([])
        for j in range(0,len(massList)):
            if j == i:
                continue                
            temp_ax = G*massList[j]/((xPos[j]-xPos[i])**2)
            if xPos[j]-xPos[i]>0:
                temp_ax_array = np.append(temp_ax_array,temp_ax)
            else:
                temp_ax_array = np.append(temp_ax_array,-1*temp_ax)  
        
            temp_ay = G*massList[i]*massList[j]/((yPos[j]-yPos[i])**2)
            if yPos[j]-yPos[i]>0:
                temp_ay_array = np.append(temp_ay_array,temp_ay)
            else:
                temp_ay_array = np.append(temp_ay_array,-1*temp_ay)
            
            temp_az = G*massList[i]*massList[j]/((zPos[j]-zPos[i])**2)
            if zPos[j]-zPos[i]>0:
                temp_az_array = np.append(temp_az_array,temp_az)
            else:
                temp_az_array = np.append(temp_az_array,-1*temp_az)
        a_x = np.append(a_x,np.mean(temp_ax_array))
        a_y = np.append(a_y,np.mean(temp_ay_array))
        a_z = np.append(a_z,np.mean(temp_az_array))
       
    xVel = xVel + a_x
    yVel = yVel + a_y
    zVel = zVel + a_z
    return xVel,yVel,zVel
    
def applyBoundaryConditions(xPos,yPos,zPos,xVel,yVel,zVel):
    for i in range(0,len(massList)):
        if xPos[i] > radius:
            xVel[i] = -1*xVel[i]
            xPos = radius
        if xPos[i] < -1*radius:
            xVel[i] = -1*xVel[i]
            xPos = -1*radius

        if yPos[i] > radius:
            yVel[i] = -1*yVel[i]
            yPos = radius
        if yPos[i] < -1*radius: 
            yVel[i] = -1*yVel[i]
            yPos = -1*radius
            
        if zPos[i] > radius:
            zVel[i] = -1*zVel[i]
            zPos = radius
        if zPos[i] < -1*radius:
            zVel[i] = -1*zVel[i]
            zPos = -1*radius
    return xPos,yPos,zPos,xVel,yVel,zVel

for i in range(0,numSim):
    if i % 10 == 0:
        print('Running interation', i)
    xVel,yVel,zVel = updateVelocities(xVel,yVel,zVel)
    xPos = xPos + xVel
    yPos = yPos + yVel
    zPos = zPos + zVel
    xPos,yPos,zPos,xVel,yVel,zVel = applyBoundaryConditions(xPos,yPos,zPos,xVel,yVel,zVel)
    
circVelocity = np.sqrt(xVel**2+yVel**2+zVel**2)
rad = np.sqrt(xPos**2+yPos**2+zPos**2)
    
dotProd = xPos*xVel + yPos*yVel +zPos*zVel    

#sorts according to the radius
for iter_num in range(len(rad)-1,0,-1):
    for idx in range(iter_num):
        if rad[idx]>rad[idx+1]:
            
            temp = rad[idx]
            rad[idx] = rad[idx+1]
            rad[idx+1] = temp
            
            temp2 = circVelocity[idx]
            circVelocity[idx] = circVelocity[idx+1]
            circVelocity[idx+1] = temp2
    
#compute the moving averages
radius_movingAvg = np.array([])
velocity_movingAvg = np.array([])

print('Computing moving averages')
t = 10 #moving average variable     
for i in range(len(rad)-t):
    velocity_temp_array = circVelocity[i:i+t]
    radius_temp_array = rad[i:i+t]
    radius_movingAvg = np.append(radius_movingAvg,np.mean(radius_temp_array))
    velocity_movingAvg = np.append(velocity_movingAvg,np.mean(velocity_temp_array))

plt.plot(radius_movingAvg, velocity_movingAvg,'ro')
plt.xlabel('Radius')
plt.ylabel('Circular Velocity')