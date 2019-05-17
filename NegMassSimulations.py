import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

#intial variables
numPosMass= 50
numNegMass= 270

#the ratio of dark matter to normal matter is roughly 27:5
radius = 2500
G = 1.0
numSim = 1
timeStep = 0.01
guassVeloComp = 0.3
t = 5 #moving average variable     

totalParticles = numPosMass+numNegMass
massList = np.array([])
#sets the masses
for i in range(0,numPosMass):
    massList = np.append(massList,0.01)

for i in range(0,numNegMass):
    massList = np.append(massList,-1*0.01)

#Puts the objects in random places
xPos = np.array(np.random.uniform(-radius,radius,totalParticles))
yPos = np.array(np.random.uniform(-radius,radius,totalParticles))
zPos = np.array(np.random.uniform(-radius,radius,totalParticles))

#For each object, finds two different random angles for velocity direction
phi = np.random.uniform(0, 2*np.pi, numPosMass)
theta = np.random.uniform(0,(np.pi)/2, numPosMass)
#finds what the circlar velocity should be given position
cirVel = np.sqrt(G*massList[0:numPosMass]/np.sqrt(xPos[0:numPosMass]**2+yPos[0:numPosMass]**2+zPos[0:numPosMass]**2))
#finds the cartesian coponents of velocties
xVel = cirVel*np.sin(theta)*np.cos(phi) + np.random.normal(0.0, guassVeloComp, 1)
yVel = cirVel*np.sin(theta)*np.sin(phi) + np.random.normal(0.0, guassVeloComp, 1)
zVel = cirVel*np.cos(theta) + np.random.normal(0.0, guassVeloComp, 1)

#negative mass particles start at rest
for i in range(0,numNegMass):
    xVel = np.append(xVel,0)
    yVel = np.append(yVel,0)
    zVel = np.append(zVel,0)


a_x = np.zeros(totalParticles)
a_y = np.zeros(totalParticles)
a_z = np.zeros(totalParticles)

def updateVelocities(xVel,yVel,zVel,a_x,a_y,a_z):
    a_x_new = np.array([])
    a_y_new = np.array([])
    a_z_new = np.array([])
    #goes through all the object
    for i in range(0,len(massList)):       
        temp_ax_array = np.array([])
        temp_ay_array = np.array([])
        temp_az_array = np.array([])
        #for a given object, finds the accelation due to all other objects
        for j in range(0,len(massList)):
            #skips the iteration if it is the object itself
            if j == i:
                continue                
            temp_ax = G*massList[j]/((xPos[j]-xPos[i])**2)
            #if the object is on the right, then acceleration is negative
            if xPos[j]-xPos[i]>0:
                temp_ax_array = np.append(temp_ax_array,-1*temp_ax)
            #if the object is on the left, then acceleration is positive
            else:
                temp_ax_array = np.append(temp_ax_array,temp_ax)  
        
            temp_ay = G*massList[j]/((yPos[j]-yPos[i])**2)
            if yPos[j]-yPos[i]>0:
                temp_ay_array = np.append(temp_ay_array,-1*temp_ay)
            else:
                temp_ay_array = np.append(temp_ay_array,temp_ay)
            
            temp_az = G*massList[j]/((zPos[j]-zPos[i])**2)
            if zPos[j]-zPos[i]>0:
                temp_az_array = np.append(temp_az_array,-1*temp_az)
            else:
                temp_az_array = np.append(temp_az_array,temp_az)
                
        #calculates the net acceleration of a given object due to all other        
        a_x_new = np.append(a_x_new,np.mean(temp_ax_array))
        a_y_new = np.append(a_y_new,np.mean(temp_ay_array))
        a_z_new = np.append(a_z_new,np.mean(temp_az_array))
    
    #calculates the velocities using Velocity Verlet algorithym    
    xVel = xVel + ((a_x + a_x_new)*timeStep)/2
    yVel = yVel + ((a_y + a_y_new)*timeStep)/2
    zVel = zVel + ((a_z + a_z_new)*timeStep)/2
    
    return xVel,yVel,zVel,a_x_new,a_y_new,a_z_new

#applies boundary conditions and makes sure that the object is inside the box 
#if the object is outside the box, flips the velocity and brings it back to the edge
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

for i in tqdm(range(0,int(numSim/timeStep))):
    xVel,yVel,zVel,a_x,a_y,a_z = updateVelocities(xVel,yVel,zVel,a_x,a_y,a_z)
    #calculates the position using Velocity Verlet algorithym   
    xPos = xPos + xVel*timeStep + (a_x*timeStep**2)/2
    yPos = yPos + yVel*timeStep + (a_y*timeStep**2)/2
    zPos = zPos + zVel*timeStep + (a_z*timeStep**2)/2
    xPos,yPos,zPos,xVel,yVel,zVel = applyBoundaryConditions(xPos,yPos,zPos,xVel,yVel,zVel)
    

pos_x,pos_y,pos_z,pos_xVel, pos_yVel, pos_zVel = np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

#Makes a list of the position and velocities of objects with positive mass

pos_x = xPos[0:numPosMass]
pos_y = yPos[0:numPosMass]
pos_z = zPos[0:numPosMass]

pos_xVel = xVel[0:numPosMass]
pos_yVel = yVel[0:numPosMass]
pos_zVel = zVel[0:numPosMass]
        
rad = np.sqrt(pos_x**2+pos_y**2+pos_z**2)

#r_sqrd = pos_x**2+pos_y**2+pos_z**2

#xcomp = np.array((pos_x**2*pos_xVel + pos_x*pos_y*pos_yVel + pos_x*pos_z*pos_zVel)/r_sqrd - pos_xVel)
#ycomp = np.array((pos_x*pos_y*pos_xVel + pos_y**2*pos_yVel + pos_y*pos_z*pos_zVel)/r_sqrd - pos_yVel)
#zcomp = np.array((pos_x*pos_z*pos_xVel + pos_y*pos_z*pos_yVel + pos_z**2*pos_zVel)/r_sqrd - pos_zVel)

circVelocity = np.array((pos_xVel**2 + pos_yVel**2 + pos_zVel**2)**0.5)

#sorts according to the radius
print('Sorting by radius')

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

for i in range(len(rad)-t):
    velocity_temp_array = circVelocity[i:i+t]
    radius_temp_array = rad[i:i+t]
    radius_movingAvg = np.append(radius_movingAvg,np.mean(radius_temp_array))
    velocity_movingAvg = np.append(velocity_movingAvg,np.mean(velocity_temp_array))

test_vel = np.sqrt(G*0.01/radius_movingAvg)

plt.plot(radius_movingAvg, velocity_movingAvg,'ro',label='calculated')
#plt.plot(radius_movingAvg,test_vel,'bo',label='actual')
plt.legend(loc='upper right')
plt.xlabel('Radius')
plt.ylabel('Circular Velocity')

##find the density
#s = int(len(rad)/20)
#
#location = np.array([])
#density = np.array([])
#counter = 0
#t=1
#for i in range(0,len(radius_movingAvg)-s):
#     for i in range(s,s+20) :
#         if radius_movingAvg[i] <= radius*t/s:
#             counter = counter + 1
#     density = np.append(density, counter) 
#     location = np.append(location,np.mean(radius_movingAvg[i:i+s]))
#     counter = 0
#     i = i+s
#     t = t+1
#        
#plt.plot(location, density,'go')

