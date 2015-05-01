#encoding=utf8
import math
import random
import odsreader

class Curve(object):
    def __init__(self,points): #points are tuple  (pressure,flow)
        self.points=sorted(points,key=lambda x:x[0])
    def get_flow(self,pressure):
        for p1,p2 in zip(self.points,self.points[1:]):
            pressure1,flow1=p1
            pressure2,flow2=p2
            
            if pressure>=pressure1 and pressure<=pressure2:
                k=(pressure-pressure1)/(pressure2-pressure1)
                return flow1 + k*(flow2-flow1)
        minpressure,minpressureflow=self.points[0]
        maxpressure,maxpressureflow=self.points[-1]
        if pressure>=maxpressure:
            return maxpressureflow
        if pressure<=minpressure:
            return minpressureflow        
        raise Exception("Can't get flow for pressure %f"%(pressure,)) 


class VirtualPoint(object):
    def __init__(self,name):
        self.name=name
        self.pipes=[]
        self.pipes_A=[]
        self.pipes_B=[]
        self.flowsum=0.0
        self.pressure=0.0
        self.zero=False
    def __repr__(self):
        return "Point(%s - %s)"%(self.name,id(self))
        
    def connect(self,pipe_end):
        self.pipes.append(pipe_end.pipe)
        if pipe_end.idx==0:
            self.pipes_A.append(pipe_end.pipe)
        else:
            self.pipes_B.append(pipe_end.pipe)
        pipe_end.connected_point=self
        
        
class PipeEnd(object):
    def __init__(self,pipe,idx):
        self.pipe=pipe
        self.idx=idx
        self.connected_point=None

fluid_kinematic_viscosity=9.7937E-7
fluid_density=(789*0.3 + 1000.0*0.7)
        

class FanCoil(object):
    def __init__(self,name,flow_lph,pressure_kpa,maxpower,stryp=None,valve_kvs=None):
        self.name=name
        self.A=PipeEnd(self,0)
        self.B=PipeEnd(self,1)
        self.equivdists=0.0
        self.flow=0.0
        self.stryp=stryp
        self.maxpower=maxpower
        self.valve_kvs=valve_kvs
        
        #flow=math.sqrt(abs(pdiff)*self.k) #P in bar
        #flow**2 = pdiff*self.k
        #k = flow**2/pdiff
        flow=flow_lph/1000.0/3600.0
        pressure=pressure_kpa*1000.0
        self.k=flow**2.0/pressure
        self.maxflow=flow
        self.valveloss=0.0
    def printstate(self):
        ratio=self.flow/self.maxflow
        pdiff=(self.A.connected_point.pressure-self.B.connected_point.pressure)        
        print "Power of %s: %fW (%.0f%%) (dP: %f kPa, dP(valve): %f kPa, flow: %f l/h)"%(self.name,ratio*self.maxpower,100.0*ratio,
            (pdiff-self.valveloss)/1e3,self.valveloss/1e3,self.flow*1000*3600.0)
    def __repr__(self):
        return "Fancoil(%s)"%(self.name,)
    def calcflow(self,pdiff):
        
        if abs(pdiff)<1e-4:
            return 0.0
        else:
            valveloss=0
            if self.valve_kvs and abs(self.flow)>1e-20:
                # Kvs = F*sqrt(SG/dP)
                # Kvs/F = sqrt(SG/dP)
                # (Kvs/F)**2.0 = SG/dP
                # dP*(Kvs/F)**2.0 = SG
                # dP = SG/((Kvs/F)**2.0)
                SG = fluid_density/1000.0
                Kvs = self.valve_kvs
                F = self.flow*3600
                valveloss=SG/((Kvs/F)**2.0) * 100.0e3
                self.valveloss=valveloss
            flow=math.sqrt((abs(pdiff-valveloss))*self.k)
            #if flow>self.maxflow:
            #    flow=self.maxflow
            if self.stryp:
                flow*=(100.0-self.stryp)/100.0
                
            if pdiff<0:
                return -flow
            else:
                return flow
                
        
class Pipe(object):
    def __init__(self,name,length,diameter,roughness=0.03e-3):
        self.name=name
        self.A=PipeEnd(self,0)
        self.B=PipeEnd(self,1)
        self.equivdists=0.0
        self.flow=0.0
        self.length=length
        self.diameter=diameter
        self.roughness=roughness

    def __repr__(self):
        return "Point(%s)"%(self.name,)
    def calcflow(self,pdiff):
        r=self.diameter/2.0
        area=math.pi*r**2.0        
        flowspeed = abs(self.flow/area)
        diameter=self.diameter
        
        flowspeed=solve_flow(pdiff,diameter,self.length+self.equivdists,flowspeed,self.roughness)
        
        if pdiff<0:
            return -flowspeed*area
        else:
            return flowspeed*area

def solve_flow(pdiff,diameter,length,flowspeed,roughness=0.03e-3):
    r=diameter/2.0
    area=math.pi*r**2.0        
    reynolds=  flowspeed * diameter /(fluid_kinematic_viscosity)
    
    if reynolds<200:
        reynolds=200
    
    A=(-2.457*math.log( ( ((7.0/reynolds)**0.9 + 0.27*(roughness/diameter))  )) )**16.0

    B=(37530.0/reynolds)**16.0

    f = 8 * (((8.0/reynolds)**12.0+(1.0/((A+B)**1.5))) ** (1/12.0))
    
    
    flowspeed = math.sqrt(abs(pdiff) / ((f * length / diameter  *  fluid_density)/2.0) )
    
    
    return flowspeed
    
    
class FittingEndpoint(object):
    def __init__(self,kind):
        self.kind=kind
    
class Fitting(object):
    pass

class Joint(Fitting):
    def __init__(self,name,*pipes):
        self.point=VirtualPoint(name)
        for pipe in pipes:
            if not isinstance(pipe,PipeEnd):
                raise Exception("TPiece must connect to objects of type PipeEnd")
            self.point.connect(pipe)



class Pump(object):
    def __repr__(self):
        return "Pump(%s)"%(self.name,)
    def setzero(self):
        self.A.connected_point.zero=True
        self.A.connected_point.pressure=0.0

    def __init__(self,name,curve,zero=False):
        self.name=name
        self.flow=0.0

        self.zero=zero
        self.equivdists=0.0
           
        self.A=PipeEnd(self,0)
        self.B=PipeEnd(self,1)
        
        self.curve=curve
    def calcflow(self,pdiff):

        flow=self.curve.get_flow(-pdiff)
        if self.equivdists>0: raise Exception("This program does not support connecting pumps directly to T-junctions or elbows. You must connect it to a short straight pipe first (may use 0 length)")
        return flow
        
def find_all_pipes(point,seen_points,seen_pipes):
    if point in seen_points:
        return
    seen_points.add(point)

    for pipe in point.pipes:
        seen_pipes.add(pipe)
        if pipe.A.connected_point==None:
            print "Pipe:",pipe.name
            raise Exception("A-end of pipe %s is unconnected"%(pipe.name,))
        if pipe.B.connected_point==None:
            print "Pipe:",pipe.name
            raise Exception("A-end of pipe %s is unconnected"%(pipe.name,))
        find_all_pipes(pipe.A.connected_point,seen_points,seen_pipes)
        find_all_pipes(pipe.B.connected_point,seen_points,seen_pipes)
        

class Simulator(object):
    def simulate(self,pipes):
        if len(pipes)==0: raise Exception("Must have at least one pipe")
        firstpipe=pipes[0]
        
        expectedpipeset=set(pipes)
        actualpipeset=set()
        find_all_pipes(firstpipe.A.connected_point,set(),actualpipeset)
        if actualpipeset!=expectedpipeset:
            print "Actual set:",actualpipeset
            print "Expected set:",expectedpipeset
            print "Extra: ",expectedpipeset-actualpipeset
            print "Missing: ",actualpipeset-expectedpipeset
            raise Exception("All pipes must be connected together, and all connected pipes must be given to this method")
            
        allpipes=list(actualpipeset)
        
        pointset=set()
        for pipe in allpipes:
            pointset.add(pipe.A.connected_point)
            pointset.add(pipe.B.connected_point)

        allpoints=list(pointset)
        
        
        zeropoint=[x for x in allpoints if x.zero]
        if len(zeropoint)!=1:
            pumps=[x for x in allpipes if isinstance(x,Pump) and x.zero]
            if len(pumps)==0:       
                pumps=[x for x in allpipes if isinstance(x,Pump)]
            if len(pumps)==1:
                pumps[0].setzero()
            else:
                raise Exception("There must be exactly one point in the graph with pressure defined to zero. Set zero-flag on one of the pumps")        
        pressure_adjust=1e6
        orig_factor=pressure_adjust
        cowardice=0.1
        last_residual=1e10
        kfactor=0.25
        for iter in xrange(50000):
            
            for point in allpoints:
                point.flowsum=0.0

            for pipe in allpipes:
                pdiff=(pipe.A.connected_point.pressure-pipe.B.connected_point.pressure)        
                flow=pipe.calcflow(pdiff)
                epsilon=10.0
                flowr=pipe.calcflow(pdiff + epsilon)
                flowl=pipe.calcflow(pdiff - epsilon)
                
                pipe.derivative=(flowr-flowl)/(2.0*epsilon)
                pipe.derivative_r=(flowr-flow)/epsilon
                pipe.derivative_l=(flow-flowl)/epsilon                                
                
                pipe.flow=flow
                pipe.A.connected_point.flowsum-=pipe.flow
                pipe.B.connected_point.flowsum+=pipe.flow

            #for point in allpoints:
            if True:
                kfun=lambda k:abs(k.flowsum) if not k.zero else -1e60
                point=max(allpoints,key=kfun)
                biggest_flowsum=abs(point.flowsum)
                biggest_flownode=point
                #print "Chose node:",point.name,point.flowsum
                if point.zero:
                    point.pressure=0.0
                else:
                    
                    draining_pipes=point.pipes_A
                    sourcing_pipes=point.pipes_B
                    
                    #P = k * flow**2 => flow = ((P-Pn)/K)**0.5
                    
                    # flowsum = sum(n,flow[n])
                    # flowsum = sum(n,(P/K-P[n]/K[n])**0.5)
                    # flowsum = sum(n,(P/K-P[n]/K[n])**0.5)
                        
                    
                    epsilon_flow=0.0
                    doprint=False
                    #if point.name=="K6" or 1:

                    if point.flowsum>0:
                        #pressure is too low                        
                        for drain_pipe in draining_pipes:
                            epsilon_flow+=drain_pipe.derivative
                            if doprint:print "incpresdrain Pipe:",drain_pipe,"w/ flow",drain_pipe.flow,"Contributes",drain_pipe.derivative_r,"to epsflow"
                        for source_pipe in sourcing_pipes:
                            epsilon_flow+=source_pipe.derivative
                            if doprint:print "incpressource Pipe:",source_pipe,"w/ flow",source_pipe.flow,"Contributes",source_pipe.derivative_l,"to epsflow"
                    if point.flowsum<0:
                        #pressure is too high
                        
                        for drain_pipe in draining_pipes:
                            epsilon_flow+=drain_pipe.derivative
                            if doprint:print "lowerpresdrain Pipe:",drain_pipe,"w/ flow",drain_pipe.flow,"Contributes",drain_pipe.derivative_l,"to epsflow"
                        for source_pipe in sourcing_pipes:
                            epsilon_flow+=source_pipe.derivative
                            if doprint:print "lowerpressource Pipe:",source_pipe,"w/ flow",source_pipe.flow,"Contributes",source_pipe.derivative_r,"to epsflow"
                            
                    if abs(epsilon_flow)>1e-15:
                        one_pascal_change=epsilon_flow
                        pressure_change=.5*point.flowsum/epsilon_flow
                        if pressure_change>1e3:
                            pressure_change=1e3
                        if pressure_change<-1e3:
                            pressure_change=-1e3
                        if doprint:print point.name,"Flowsum:",point.flowsum*3600*1000.0,"epsilon flow:",epsilon_flow,"pressure:",point.pressure
                        point.pressure+=pressure_change
                        if doprint or True:print "Decided to change pressure by",pressure_change,"to",point.pressure
                        
                        for x in xrange(1):
                            for point1 in allpoints:
                                point1.flowsum=0.0

                            for pipe in allpipes:
                                pdiff=(pipe.A.connected_point.pressure-pipe.B.connected_point.pressure)        
                                flow=pipe.calcflow(pdiff)
                                epsilon=10.0
                                flowr=pipe.calcflow(pdiff + epsilon)
                                flowl=pipe.calcflow(pdiff - epsilon)
                                
                                pipe.derivative=(flowr-flowl)/(2.0*epsilon)
                                pipe.derivative_r=(flowr-flow)/epsilon
                                pipe.derivative_l=(flow-flowl)/epsilon                                
                                
                                pipe.flow=flow
                                pipe.A.connected_point.flowsum-=pipe.flow
                                pipe.B.connected_point.flowsum+=pipe.flow
                                if pipe.A.connected_point.name==point.name:
                                    print "flowsum:",pipe.A.connected_point.flowsum*3600.0*1000.0 ," abs flow:",pipe.flow*3600.0*1000
                        print "Result:",point.flowsum*3600.0*1000.0
                    else:
                        pass#print "Epsilon flow of point was zero, flowsum was:",point.flowsum
                        
                    
                    #if abs(point.flowsum)>biggest_flowsum:
                    #    biggest_flownode=point
                    #    biggest_flowsum=abs(point.flowsum)
                    
                    #adjust=(0.5+0.5*random.random())*pressure_adjust*point.flowsum
                    #point.pressure+=adjust
                
            cowardice+=0.1
            if cowardice>1.0:
                cowardice=1.0
            print "Flowsum:",biggest_flowsum*1000.0*3600.0,"l/h",biggest_flownode.name
            print "Residual:",biggest_flowsum,"Last residual:",last_residual
            if biggest_flowsum>last_residual:                
                kfactor*=0.8
            else:
                kfactor+=0.1
            print "Kfacotr:",kfactor
            last_residual=(biggest_flowsum)
            #if biggest_flowsum*3600.0*1000.0<0.1:
            #    break
            pressure_adjust/=1.0001
        print "End adjust:",pressure_adjust/orig_factor
        overspeed=[]
        for pipe in allpipes:
            pressure=pipe.B.connected_point.pressure-pipe.A.connected_point.pressure
            pressure_kpa=pressure/1000.0
            cooling_power=pipe.flow*1000.0 * 4185 * (12.0-7.0)

            if hasattr(pipe,'diameter'):
                r=pipe.diameter/2.0
                area=math.pi*r**2.0        
                flowspeed=pipe.flow/area
                if flowspeed>1.4:
                    overspeed.append((pipe.name,flowspeed))
                reynolds=  flowspeed * pipe.diameter /(fluid_kinematic_viscosity)
                flowspeed = "%f m/s (re: %f)"%(flowspeed,reynolds)
            else:
                flowspeed=""
            print "Flow in pipe %s: %f l/h (dP: %f kPa, %f kW, %s)"%(pipe.name,pipe.flow*1000.0*3600.0,pressure_kpa,cooling_power,flowspeed)
        for point in allpoints:
            print "Pressure at point %s: %f kPa (residual: %f l/h)"%(point.name,point.pressure/1000.0,point.flowsum*1000.0*3600.0)


        print
        for name,speed in overspeed:
            print "OVERSPEED: Pipe %s, speed = %f m/s"%(name,speed)

        print                        
    
def kpa_lph_to_si_units(points):
    for pressure_kpa,flow_lph in points:
        yield (pressure_kpa*1000.0,flow_lph/1000.0 / 3600.0)
def kpa_lpm_to_si_units(points):
    for pressure_kpa,flow_lpm in points:
        yield (pressure_kpa*1000.0,flow_lpm/1000.0 / 60.0)
        
      
