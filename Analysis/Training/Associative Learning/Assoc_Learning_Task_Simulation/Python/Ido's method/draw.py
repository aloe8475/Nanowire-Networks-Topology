import numpy as np
import plotly.graph_objects as go
from plotly import figure_factory as ff

def draw_wires(fig, toDraw, TimeStamp=0, PathInd=[]):
    wires_toDraw = []
    sectionCenterX = np.array([])
    sectionCenterY = np.array([])
    sectionCurrentX = np.array([])
    sectionCurrentY = np.array([])

    for this_wire in toDraw.Nanowires:
        if this_wire.index == 0:
            continue
        if not this_wire in wires_toDraw:
            wires_toDraw.append(this_wire)
            sectionCenterX = np.append(sectionCenterX, this_wire.sectionCenterX)
            sectionCenterY = np.append(sectionCenterY, this_wire.sectionCenterY)
            sectionCurrentX = np.append(sectionCurrentX, this_wire.sectionCurrentX[TimeStamp,:])
            sectionCurrentY = np.append(sectionCurrentY, this_wire.sectionCurrentY[TimeStamp,:])

    for this_wire in wires_toDraw:
        xa, ya, xb, yb = this_wire.wireEnds
        xc, yc = this_wire.centerPosition
        this_opacity = 0.2 + 0.8*(this_wire.onPathChecker[TimeStamp] != 0)
        if this_wire.index in PathInd:
            this_color = 'yellow'
        else:
            this_color = 'white'

        line = go.Scatter(x=[xa,xb,xc],
                        y=[ya,yb,yc],
                        name = 'Nanowire #'+str(this_wire.index),
                        mode="lines",
                        hoverinfo='name',
                        opacity=this_opacity,
                        line=dict(color=this_color,
                                width=1.5))

        fig.add_trace(line)

        if this_wire.isElectrode != 0:
            if this_wire.isElectrode == 1:
                electrodeColor = 'green'
            else:
                electrodeColor = 'red'
            temp = go.Scatter(x = [this_wire.contactEnd[0]],
                            y = [this_wire.contactEnd[1]],
                            mode='markers',
                            name='Drain',
                            opacity=0.9,
                            marker=dict(symbol = 'hexagon',
                                        color=electrodeColor,
                                        size=30),
                            hoverinfo = 'name'
                    )
            fig.add_trace(temp)

    quivers = ff.create_quiver(sectionCenterX, sectionCenterY,
                                1e8*sectionCurrentX , 1e8*sectionCurrentY,
                                scale= 0.1,
                                arrow_scale=0.3,
                                name='Current',
                                opacity = 0.9,
                                hoverinfo = 'none',
                                line=dict(color='red',
                                        width=2))
    fig.add_trace(quivers.data[0])
    return fig
    
def draw_junctions(fig, toDraw, TimeStamp=0, JunctionsToObserve=[]):
    junctions_toDraw = []

    for this_junction in toDraw.Junctions:
        if this_junction.index == 0:
            continue
        if not this_junction in junctions_toDraw:
            junctions_toDraw.append(this_junction)
    
    for this_junction in junctions_toDraw:
        xi, yi = this_junction.position
        if this_junction.onCurrentPath[TimeStamp]:
            this_color = 'red'
        elif this_junction.switch[TimeStamp]:
            this_color = 'green'
        else:
            this_color = 'white'

        temp = go.Scatter(x = [xi],
                            y = [yi],
                            mode='markers',
                            name='Junction #' + str(this_junction.index),
                            opacity=0.9,
                            marker=dict(color=this_color,
                                        size=5),
                            hoverinfo = 'name'
                )
        fig.add_trace(temp)

        if this_junction.index in JunctionsToObserve:
            temp = go.Scatter(x = [xi],
                            y = [yi],
                            mode='markers',
                            name='Junction #' + str(this_junction.index),
                            opacity=0.7,
                            marker=dict(symbol = 'star',
                                        color='orange',
                                        size=15),
                            hoverinfo = 'name'
                )
            fig.add_trace(temp)
        
    return fig

def draw(toDraw, PathHighlight=[], JunctionsToObserve=[], **kwargs):
    Lx, Ly = toDraw.gridSize

    if len(PathHighlight) != 0:
        PathInd = PathHighlight
    elif hasattr(toDraw, 'shortestPaths'):
        PathInd = toDraw.shortestPaths[0]
    else:
        PathInd = []

    for i in PathInd:
        if i not in toDraw.WiresInd:
            print(f'Nanowire #{i} is not in this region!')

    for i in JunctionsToObserve:
        if i not in toDraw.JunctionsInd:
            print(f'Junction #{i} is not in this region.')
    
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

    fig = go.Figure()
    draw_wires(fig, toDraw, TimeStamp = this_TimeStamp, PathInd = PathInd)
    draw_junctions(fig, toDraw, TimeStamp = this_TimeStamp, JunctionsToObserve = JunctionsToObserve)  
    fig.update_layout({'height':800, 'width':800},
                    showlegend=False,
                    plot_bgcolor='rgb(51,51,51)',
                    paper_bgcolor='rgb(211,211,211)')
    fig.update_xaxes(title_text=r'x $\mu$ m',
                    range=[-.2*Lx,1.2*Lx],
                    showgrid=False,
                    zeroline=False)
    fig.update_yaxes(title_text=r'x $\mu$ m',
                    range=[-.2*Lx,1.2*Lx],
                    showgrid=False,
                    zeroline=False)

    return fig

if __name__ == "__main__":
    from utils import *

    SimulationOptions = simulation_options__(dt = 1e-3, T = 3,
                                        interfaceContactWires = [73, 30])

    Connectivity = connectivity__(
        filename = '2016-09-08-155153_asn_nw_00100_nj_00261_seed_042_avl_100.00_disp_10.00.mat')

    JunctionState = junctionState__(Connectivity.numOfJunctions)
    Stimulus = stimulus__(biasType = 'DC', 
                        TimeVector = SimulationOptions.TimeVector, 
                        onTime = 0, offTime = 1,
                        onAmp = 1.1, offAmp = 0.005)

    Equations = equations__(Connectivity, SimulationOptions.contactWires)
    this_realization = simulateNetwork(SimulationOptions, Connectivity, JunctionState, 
                                        Stimulus, Equations, 
                                        simpleOutput=False,
                                        useSparse=False)
    this_realization.allocateData()
    draw(this_realization, time=0.9)