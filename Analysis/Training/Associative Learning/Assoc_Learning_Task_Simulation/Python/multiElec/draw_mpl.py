import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def draw_mpl(toDraw, **kwargs):
    if 'TimeStamp' in kwargs:
            this_TimeStamp = kwargs['TimeStamp']
    elif 'time' in kwargs:
        if kwargs['time'] in toDraw.TimeVector:
            this_TimeStamp = np.where(toDraw.TimeVector == kwargs['time'])[0][0]
        elif (kwargs['time'] < min(toDraw.TimeVector)) or (kwargs['time'] > max(toDraw.TimeVector)):
            print('Input time exceeds simulation period.')
            this_TimeStamp = np.argmin(abs(toDraw.TimeVector - kwargs['time']))
        else:
            this_TimeStamp = np.argmin(abs(toDraw.TimeVector - kwargs['time']))
    else:
        this_TimeStamp = 0
    this_TimeStamp = 901
    """
    Set up the figrue
    """
    plt.style.use('classic')
    fig, ax = plt.subplots()
    fig.set_size_inches(10,10)
    fig.set_facecolor((0.8,0.8,0.8))
    ax.set_facecolor((0.2,0.2,0.2))
    ax.set_title(f'time = {toDraw.TimeVector[this_TimeStamp]}')
    Lx, Ly = toDraw.connectivity.length_x, toDraw.connectivity.length_y
    ax.axis([-.2*Lx,1.2*Lx,-.2*Lx,1.2*Lx])
    
    """
    Plot junctions
    """
    xi, yi = toDraw.connectivity.xi, toDraw.connectivity.yi
    junctions = Line2D([xi],[yi], 
                    color='white', marker='o', 
                    ms=5, alpha=0.8,
                    ls='None')
    ax.add_line(junctions) 
    
    """
    Plot nanowires
    """
    xa,ya = toDraw.connectivity.xa, toDraw.connectivity.ya
    xb,yb = toDraw.connectivity.xb, toDraw.connectivity.yb
    for this_wire in range(xa.size):
        nanowires = Line2D([xa[this_wire],xb[this_wire]],[ya[this_wire],yb[this_wire]], 
                            color='white', alpha = 0.5)
        ax.add_line(nanowires)

    """
    Plot currents
    """
    sources = np.add(toDraw.sources, -1)
    drains = np.add(toDraw.drains, -1)
    sourcePoints = np.array([])
    drainPoints = np.array([])
    xc,yc = toDraw.connectivity.xc, toDraw.connectivity.yc
    adjMat = toDraw.adjMat
    edgeList = toDraw.connectivity.edge_list
    current = toDraw.junctionVoltage[this_TimeStamp,:]/toDraw.junctionResistance[this_TimeStamp,:]
    arrowCenterX = np.array([])
    arrowCenterY = np.array([])
    arrowX = np.array([])
    arrowY = np.array([])
    sectionCurrent = np.array([])

    wireAngles = np.arctan2(yb-ya, xb-xa) % np.pi
    wireAngles = np.where(wireAngles <= np.pi/2, wireAngles, wireAngles - np.pi)

    for this_wire in range(xc.size):
        contactWires = np.where(adjMat[this_wire,:]!=0)[0]
        contactJunctions = np.where((edgeList[:,0]==this_wire)+(edgeList[:,1]==this_wire))[0]
        if xa[this_wire] != xb[this_wire]:
            sortedInd = np.argsort(xi[contactJunctions])
        else:
            sortedInd = np.argsort(yi[contactJunctions])
        sortedContactWires = contactWires[sortedInd]
        sortedContactJunctions = contactJunctions[sortedInd]
        
        direction = ((sortedContactWires < this_wire) - 0.5) * 2
        wireCurrents = np.zeros(contactJunctions.size-1)
        wireCurrents = np.cumsum(current[sortedContactJunctions[0:-1]]*direction[0:-1])
        arrowCenterX = np.append(arrowCenterX, np.mean([2*xi[sortedContactJunctions[0:-1]], xi[sortedContactJunctions[1:]]], axis=0)/1.5)
        arrowCenterY = np.append(arrowCenterY, np.mean([2*yi[sortedContactJunctions[0:-1]], yi[sortedContactJunctions[1:]]], axis=0)/1.5)
        # arrowX = np.append(arrowX, np.cos(wireAngles[this_wire]) * wireCurrents)
        # arrowY = np.append(arrowY, np.sin(wireAngles[this_wire]) * wireCurrents)
        sectionCurrent =  np.append(sectionCurrent, wireCurrents)
        arrowX = np.append(arrowX, np.sign(wireCurrents)*(xi[sortedContactJunctions[1:]]-xi[sortedContactJunctions[0:-1]]))
        arrowY = np.append(arrowY, np.sign(wireCurrents)*(yi[sortedContactJunctions[1:]]-yi[sortedContactJunctions[0:-1]]))

        if this_wire in np.append(sources, drains):
            if xa[this_wire] != xb[this_wire]:
                if xa[this_wire] < xb[this_wire]:
                    contactEnd = [xb[this_wire], yb[this_wire]]
                else:
                    contactEnd = [xa[this_wire], ya[this_wire]]
            else:
                if ya[this_wire] < yb[this_wire]:
                    contactEnd = [xb[this_wire], yb[this_wire]]
                else:
                    contactEnd = [xa[this_wire], ya[this_wire]]

            if this_wire in sources:
                sourcePoints = np.append(sourcePoints, contactEnd)
            else:
                drainPoints = np.append(drainPoints, contactEnd)

            totalCurrent = np.sum(current[sortedContactJunctions]*direction)
            arrowCenterX = np.append(arrowCenterX, np.mean([2*contactEnd[0], xi[sortedContactJunctions[-1]]])/1.5)
            arrowCenterY = np.append(arrowCenterY, np.mean([2*contactEnd[1], yi[sortedContactJunctions[-1]]])/1.5)
            arrowX = np.append(arrowX, np.sign(totalCurrent)*(contactEnd[0] - xi[sortedContactJunctions[-1]]))
            arrowY = np.append(arrowY, np.sign(totalCurrent)*(contactEnd[1] - yi[sortedContactJunctions[-1]]))
            sectionCurrent = np.append(sectionCurrent, totalCurrent)

    # normCurrent = np.exp(abs(sectionCurrent)/np.sum(abs(sectionCurrent)))-1
    # coeff = 30
    # rgba = np.zeros((arrowCenterX.size,4))
    # rgba[:,0] = 1
    # rgba[:,3] = coeff*normCurrent
    # ax.quiver(arrowCenterX, arrowCenterY, 
    #         arrowX , arrowY, 
    #         units = 'xy', scale = 2,
    #         color = rgba)
    #         # width = 0.005, alpha = 0.5)
    """
    Think about arrow scale
    """
    rgba = np.zeros((arrowCenterX.size,4))
    rgba[:,0] = 1
    rgba[:,3] = (abs(sectionCurrent) > 1e-5)*0.9
    ax.quiver(arrowCenterX, arrowCenterY, 
            arrowX , arrowY, 
            color = rgba, 
            width = 0.005)

    """
    Plot Electrodes
    """
    sourcePoints = sourcePoints.reshape(len(sources),2)
    drainPoints = drainPoints.reshape(len(drains),2)
    S = Line2D(sourcePoints[:,0], sourcePoints[:,1], 
                color='g', marker='h', 
                ms=20, alpha=0.8,
                ls='None')
    ax.add_line(S)

    D = Line2D(drainPoints[:,0], drainPoints[:,1], 
                color='r', marker='h', 
                ms=20, alpha=0.8,
                ls='None')
    ax.add_line(D)
    
    return ax