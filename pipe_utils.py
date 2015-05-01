#encoding=utf8
import math
import random
import odsreader
import scipy.optimize

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

class Elbow90deg(object): #curved, curve_radius/pipe_diameter = 1
    def __init__(self,count):
        self.LD=30.0
        self.count=count

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
            
            flow=math.sqrt((abs(pdiff))*self.k)
            #if flow>self.maxflow:
            #    flow=self.maxflow
            if self.stryp:
                flow*=(100.0-self.stryp)/100.0
                
            if pdiff<0:
                return -flow
            else:
                return flow
                
        
class Pipe(object):
    def __init__(self,name,length,diameter,restrictions=None,roughness=0.03e-3):
        self.name=name
        self.A=PipeEnd(self,0)
        self.B=PipeEnd(self,1)
        self.equivdists=0.0
        self.flow=0.0
        self.length=length
        self.diameter=diameter
        self.roughness=roughness
        self.restrictions=[]
        if restrictions!=None:
            self.restrictions=restrictions
        
    def __repr__(self):
        return "Point(%s)"%(self.name,)
    def calcflow(self,pdiff):
        r=self.diameter/2.0
        area=math.pi*r**2.0        
        flowspeed = abs(self.flow/area)
        diameter=self.diameter
        
        flowspeed=solve_flow(pdiff,diameter,self.length+self.equivdists,flowspeed,self.restrictions,self.roughness)
        
        if pdiff<0:
            return -flowspeed*area
        else:
            return flowspeed*area

def solve_flow(pdiff,diameter,length,flowspeed,restrictions,roughness=0.03e-3):
    r=diameter/2.0
    area=math.pi*r**2.0        
    oldflow=None
    
        
    
    while True:
        reynolds=  flowspeed * diameter /(fluid_kinematic_viscosity)
        
        reynoldst=reynolds
        if reynoldst<1:
            reynoldst=1
        #print "Reynolds:",reynolds
        A=(-2.457*math.log( ( ((7.0/reynoldst)**0.9 + 0.27*(roughness/diameter))  )) )**16.0

        B=(37530.0/reynoldst)**16.0

        f = 8 * (((8.0/reynoldst)**12.0+(1.0/((A+B)**1.5))) ** (1/12.0))
        
        equiv_len=0.0
        for restr in restrictions:
            if hasattr(restr,'K'):
                #K=f​*(L/D)
                #LD = K/f            
                LD = restr.K/f                            
            elif hasattr(restr,'LD'):
                LD = restr.LD
            else:
                raise Exception("Restriction must have K-value or L/D value")  
            equiv_len += LD*diameter*restr.count
                
        flowspeed = math.sqrt(abs(pdiff) / ((f * (length+equiv_len) / diameter  *  fluid_density)/2.0) )
        newreynolds=  flowspeed * diameter /(fluid_kinematic_viscosity)
        if oldflow!=None:
            if abs(oldflow)<1e-12 and abs(flowspeed)<1e-12:
                #print "Zero flow detected"
                break
            change=abs(oldflow-flowspeed)/(abs(oldflow)+abs(flowspeed))
            #print "Oldflow:",oldflow,"newflow:",flowspeed,"change ratio:",change,"old rey:",reynolds,"new rey:",newreynolds
            if change<0.00001:
                break
        oldflow=flowspeed
            
    #print "Volumetric flow:",flowspeed*area*3600
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
        
        
        zeropoints=[x for x in allpoints if x.zero]
        if len(zeropoints)!=1:
            pumps=[x for x in allpipes if isinstance(x,Pump) and x.zero]
            if len(pumps)==0:       
                pumps=[x for x in allpipes if isinstance(x,Pump)]
            if len(pumps)==1:
                pumps[0].setzero()
            else:
                raise Exception("There must be exactly one point in the graph with pressure defined to zero. Set zero-flag on one of the pumps")        

        zeropoint,=[x for x in allpoints if x.zero]
        start_flowsum=None    
        for iterionnr in xrange(50000):
            
            for point in allpoints:
                point.flowsum=0.0

            for pipe in allpipes:
                pdiff=(pipe.A.connected_point.pressure-pipe.B.connected_point.pressure)        
                flow=pipe.calcflow(pdiff)
                
                pipe.flow=flow
                pipe.A.connected_point.flowsum-=pipe.flow
                pipe.B.connected_point.flowsum+=pipe.flow

            #for point in allpoints:
            biggest_flowsum,biggest_flownode=max([(abs(p2.flowsum),p2) for p2 in allpoints],key=lambda x:x[0])
            if start_flowsum==None:
                start_flowsum=biggest_flowsum
            ln_start=-math.log(start_flowsum)
            
            max_residual=0.25/3600.0/1000.0
            
            ln_end=-math.log(max_residual)
            perc=100*((-math.log(biggest_flowsum))-ln_start)/(ln_end-ln_start)
            #print "Residual:",biggest_flowsum*3600*1000.0,biggest_flownode.name,
            print perc,"%"
            
            
            
            
            if biggest_flowsum<max_residual:
                break
            
            zeroadj=zeropoint.pressure
            if abs(zeroadj)>10000:
                for point in allpoints:
                    point.pressure-=zeroadj
            
            for point in allpoints:
                #kfun=lambda k:abs(k.flowsum) if not k.zero else -1e60
                
                #point=max(allpoints,key=kfun)
                

                #print "Chose node:",point.name,point.flowsum
                #for p in allpoints:
                #    if abs(p.flowsum)*3600*1000.0>2000:
                #        print "Big flow:",p.name,p.flowsum
                if point.zero and False:
                    point.pressure=0.0
                else:

                    manual=[False]
                    def evaluate(pressure):
                        point.pressure=pressure
                        point.flowsum=0.0
                        for pipe in point.pipes:
                            pdiff=(pipe.A.connected_point.pressure-pipe.B.connected_point.pressure)        
                            flow=pipe.calcflow(pdiff)
                            if manual[0]:
                                print "Pipe:",pipe.name,"pdiff:",pdiff,"flow:",flow*3600.0*1000
                            pipe.flow=flow
                            if pipe.A.connected_point==point:
                                pipe.A.connected_point.flowsum-=pipe.flow
                            if pipe.B.connected_point==point:
                                pipe.B.connected_point.flowsum+=pipe.flow
                        #print "F(%f) = %f"%(pressure,point.flowsum)
                        return point.flowsum*1e6
                    point.pressure=scipy.optimize.bisect(evaluate,-100e3,100e3,disp=True,xtol=1e-7)
                    #y=evaluate(point.pressure)*1e-6
                    #print "Solution quality:",y  * 1000* 3600.0,"F(%f)=%f"%(point.pressure,y*3600.0*1000.0)
                    #if abs(y*1000*3600.0)>0.01:
                    #    manual[0]=True
                    #    while True:
                    #        inp=raw_input()
                    #        if not inp: break
                    #        p=evaluate(float(inp))*1e-6
                    #        print "Result:",p*1000*3600.0
                            
                    
                    #if abs(point.flowsum)>biggest_flowsum:
                    #    biggest_flownode=point
                    #    biggest_flowsum=abs(point.flowsum)
                    
                    #adjust=(0.5+0.5*random.random())*pressure_adjust*point.flowsum
                    #point.pressure+=adjust
                
            #print "Flowsum:",biggest_flowsum*1000.0*3600.0,"l/h",biggest_flownode.name
            
            
        print "Needed %d iterations"%(iterionnr,)
        zeroadj=zeropoint.pressure            
        for point in allpoints:
            point.pressure-=zeroadj

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
            print "Flow in pipe %s: %f l/h (dP: %f kPa, %f kW, %s)"%(pipe.name,pipe.flow*1000.0*3600.0,pressure_kpa,cooling_power/1e3,flowspeed)
        for point in allpoints:
            print "Pressure at point %s: %f kPa (residual: %f l/h (zero:%s))"%(point.name,point.pressure/1000.0,point.flowsum*1000.0*3600.0,point.zero)


        print
        for name,speed in overspeed:
            print "OVERSPEED: Pipe %s, speed = %f m/s"%(name,speed)
        f=open("graph.dot","w")
        f.write("digraph {\n")
        f.write("\tcharset=\"UTF-8\";\n")
        for point in allpoints:
            f.write((u"%s [label=\"%s %.0f kPa\"];"%(point.name.replace("-",'_'),point.name.replace("-",'_'),point.pressure/1e3)).encode('utf8'))
        for pipe in sorted(allpipes,key=lambda x:-1 if isinstance(x,Pump) else 0):
            f.write((u"\t%s -> %s [label=\"%s\"];\n"%(pipe.A.connected_point.name.replace("-","_"),pipe.B.connected_point.name.replace("-","_"),
                u"%s %.0f l/h"%(pipe.name,pipe.flow*3600.0*1000.0))).encode('utf8'))
        f.write("}\n")
        f.close()
        print                        
    
def kpa_lph_to_si_units(points):
    for pressure_kpa,flow_lph in points:
        yield (pressure_kpa*1000.0,flow_lph/1000.0 / 3600.0)
def kpa_lpm_to_si_units(points):
    for pressure_kpa,flow_lpm in points:
        yield (pressure_kpa*1000.0,flow_lpm/1000.0 / 60.0)
        
      
