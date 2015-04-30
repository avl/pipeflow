#encoding=utf8
import math
import random
import odsreader

class Curve(object):
    def __init__(self,points): #points are tuple  (pressure,flow)
        self.points=sorted(points,key=lambda x:x[0])
        #print "Curvp:",self.points
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
        self.flowsum=0.0
        self.pressure=0.0
        self.zero=False
    def __repr__(self):
        return "Point(%s - %s)"%(self.name,id(self))
        
    def connect(self,pipe_end):
        self.pipes.append(pipe_end.pipe)
        pipe_end.connected_point=self
        
        
class PipeEnd(object):
    def __init__(self,pipe,idx):
        self.pipe=pipe
        self.idx=idx
        self.connected_point=None

fluid_kinematic_viscosity=9.7937E-7
#fluid_kinematic_viscosity=1e-6
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
        
        #flow=math.sqrt(abs(pdiff)*self.k)
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
        print "Power of %s: %f (dP: %f kPa, dP(valve): %f kPa, flow: %f l/h)"%(self.name,ratio*self.maxpower,
            (pdiff-self.valveloss)/1e3,self.valveloss/1e3,self.flow*1000*3600.0)
    def __repr__(self):
        return "Fancoil(%s)"%(self.name,)
    def calcflow(self):
        #print "Fancoil k:",self.k
        pdiff=(self.A.connected_point.pressure-self.B.connected_point.pressure)        
        
        #pdiff=0.02e3
        if abs(pdiff)<1e-2:
            self.flow=0
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
                #print "FloW: l/h",self.flow*3600*1000.0
                #print "Valve loss:",valveloss*1e3,"kPa"
                self.valveloss=valveloss
            flow=math.sqrt((abs(pdiff-valveloss))*self.k)
            if flow>self.maxflow:
                flow=self.maxflow
            if self.stryp:
                flow*=(100.0-self.stryp)/100.0
                
            if pdiff<0:
                self.flow=-flow
            else:
                self.flow=flow
        #print "Fancoil flow:",self.flow*3600.0*1000.0,"l/h at",pdiff/1000.0,"kPa"
                
        
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
    def calcflow(self):
        r=self.diameter/2.0
        area=math.pi*r**2.0        
        flowspeed = abs(self.flow/area)
        pdiff=(self.A.connected_point.pressure-self.B.connected_point.pressure)        
        diameter=self.diameter
        
        flowspeed=solve_flow(pdiff,diameter,self.length+self.equivdists,flowspeed,self.roughness)
        
        if pdiff<0:
            self.flow=-flowspeed*area
        else:
            self.flow=flowspeed*area
        #print "Pipe flow at pressure diff %f kPa: %f l/h"%(pdiff/1000.0,self.flow*1000.0*3600.0)
        
        #pdiff=self.A.connected_point.pressure-self.B.connected_point.pressure 
        #self.flow=pdiff*30.0 #temp         

def solve_flow(pdiff,diameter,length,flowspeed,roughness=0.03e-3):
    r=diameter/2.0
    area=math.pi*r**2.0        
    reynolds=  flowspeed * diameter /(fluid_kinematic_viscosity)
    
    if reynolds<200:
        reynolds=200
    
    #reynolds=27000
    
    #reynolds=1000.0
    #print "Reynolds:",reynolds
    #print "Roughness/diameter:",roughness/diameter
    #print "Thelog:",math.log((7/reynolds)**0.9+0.27*(roughness/diameter))
    A=(-2.457*math.log( ( ((7.0/reynolds)**0.9 + 0.27*(roughness/diameter))  )) )**16.0
    #print "A:",A
    B=(37530.0/reynolds)**16.0
    #print "A+B:",A+B
    f = 8 * (((8.0/reynolds)**12.0+(1.0/((A+B)**1.5))) ** (1/12.0))
    
    
    #print "Wrong f:",f
    #f=0.0055*(1 + (2e4 * (roughness/diameter) + 1e6/reynolds)**(1/3.0))
    
    #f=0.04
    #f=0.2095
    #print "Got f:",f
    #f=0.026
    #print "Right f:",f
    #print "Resistance coeff:",f

    # pdiff = f * self.length / self.diameter  *  fluid_density * flowspeed**2.0 / 2.0
    # pdiff = ((f * self.length / self.diameter  *  fluid_density)/2.0) * flowspeed**2.0
    # pdiff / ((f * self.length / self.diameter  *  fluid_density)/2.0) = flowspeed**2.0
    # sqrt( pdiff / ((f * self.length / self.diameter  *  fluid_density)/2.0) ) = flowspeed
    
    #print "Pdiff, bar:",pdiff/100e3
    #print "Pdiff, mbar:",1000.0*pdiff/100e3
    
    flowspeed = math.sqrt(abs(pdiff) / ((f * length / diameter  *  fluid_density)/2.0) )
    
    
    #print "Volume flow:",3600*math.pi*flowspeed* (diameter/2.0)**2.0,"m3/h","flowspeed:",flowspeed,"m/s","m^3/s:",math.pi*flowspeed* (diameter/2.0)**2.0
    return flowspeed
    
    
class FittingEndpoint(object):
    def __init__(self,kind):
        self.kind=kind
    
class Fitting(object):
    pass

"""
class Appliance(object):
    def __init__(self,name,pipes,RCK):
        self.name=name
        self.point=VirtualPoint(name)
        for pipe in pipes:
            self.point.connect(pipe)
            pipe.equivdists.append(RCK)
"""            

class TPiece(Fitting):
    def __init__(self,name,pipe1,pipe2,tee,equiv_dists=[0.0,0.0,0.0]):
        self.point=VirtualPoint(name)
        if not isinstance(pipe1,PipeEnd):
            raise Exception("TPiece must connect to objects of type PipeEnd")
        if not isinstance(pipe2,PipeEnd):
            raise Exception("TPiece must connect to objects of type PipeEnd")
        if not isinstance(tee,PipeEnd):
            raise Exception("TPiece must connect to objects of type PipeEnd")
        pipe1.pipe.equivdists+=(equiv_dists[0])#RCK
        pipe2.pipe.equivdists+=(equiv_dists[1])
        tee.pipe.equivdists+=(equiv_dists[2])
        self.point.connect(pipe1)
        self.point.connect(pipe2)
        self.point.connect(tee)

        
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

    def __init__(self,name,curve,zero=False):
        self.name=name
        self.flow=0.0

        self.equivdists=0.0
           
        self.A=PipeEnd(self,0)
        self.B=PipeEnd(self,1)
        
        self.curve=curve
    def calcflow(self):
        pdiff=self.B.connected_point.pressure-self.A.connected_point.pressure 
        self.flow=self.curve.get_flow(pdiff)
        if self.equivdists>0: raise Exception("This program does not support connecting pumps directly to T-junctions or elbows. You must connect it to a short straight pipe first (may use 0 length)")
       
        
        #print "Pump flow at pressure diff %f kPa: %f l/h"%(pdiff/1000.0,self.flow*1000.0*3600.0)
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
            pumps=[x for x in allpipes if isinstance(x,Pump)]
            if len(pumps)==1:
                pumps[0].setzero()
            else:
                raise Exception("There must be exactly one point in the graph with pressure defined to zero. Set zero-flag on one of the pumps")        
        pressure_adjust=1e6
        orig_factor=pressure_adjust
        for iter in xrange(50000):
            
            for point in allpoints:
                point.flowsum=0.0

            for pipe in allpipes:
                #print "Pipe:",pipe.name,"flow was:",pipe.flow,
                pipe.calcflow()
                #print "new flow:",pipe.flow
                pipe.A.connected_point.flowsum-=pipe.flow
                pipe.B.connected_point.flowsum+=pipe.flow
                #print "points:",pipe.A.connected_point,pipe.B.connected_point
                #print "Flowsums:",pipe.A.connected_point.flowsum,pipe.B.connected_point.flowsum
            
            for point in allpoints:
                if point.zero:
                    point.pressure=0.0
                else:
                    adjust=(0.5+0.5*random.random())*pressure_adjust*point.flowsum
                    #print "Point:",point.name,"pressure was:",point.pressure,"flow sum:",point.flowsum,"adjust:",adjust,"adjust factor:",pressure_adjust
                    point.pressure+=adjust
                    #print "new pressure:",point.pressure
                
            
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
        
        
        
def run():
    odsdoc=odsreader.ODSReader("pipeplan_smaller.ods")

    pipes=[]
    def getpipe(nr):
        return pipes[nr-1]
    for idx,row in enumerate(odsdoc.getSheet("Pipes")):
        if (not any(row)): continue

        if idx==0: continue
        number,name,descr,diam,length=row[:5]
        number=int(number)
        diam=float(diam)
        length=float(length)
        
        pipe=Pipe(name,length,diam*1e-3)
        pipes.append(pipe)

   
    x=1
    FKV=crc84=FanCoil("FK-V",1188.0,20.0,6979.0,valve_kvs=5.5,stryp=100.0*1)
    FKA=crc24a=FanCoil("FK-A",324.0,15.9,1954.0,valve_kvs=1.7,stryp=100*1)
    FKB=crc24b=FanCoil("FK-B",324.0,15.9,1954.0,valve_kvs=1.7,stryp=100*1)
    FKG=crc54a=FanCoil("FK-G",756.0,35.5,4395.0,valve_kvs=2.8,stryp=100*(0))
    FKM=crc54b=FanCoil("FK-M",756.0,35.5,4395.0,valve_kvs=2.8,stryp=100*(1))
    fancoils=[FKV,FKA,FKB,FKG,FKM]

    
    
    #PCS44 pump: pump=Pump("pump",Curve(kpa_lph_to_si_units([(0.0,2800),(53.0,520),(60.0,0)])))
    pump=Pump("Grundfos Alpha2",Curve(kpa_lph_to_si_units([(0.0,4200),(10.0,4000.0),(20.0,2900.0),(30.0,2300.0),(40.0,1600.0),(50.0,1100.0),(58.0,600.0),(60.0,0)])))
    
    K1=Joint("K1",pump.B,getpipe(1).A,getpipe(2).A,getpipe(3).A)
    K2=Joint("K2",getpipe(2).B,getpipe(7).A,getpipe(8).A)
    K3=Joint("K3",getpipe(5).A,getpipe(9).B,getpipe(10).B)
    K4=Joint("K4",getpipe(3).B,getpipe(11).A,getpipe(12).A)
    K5=Joint("K5",getpipe(13).B,getpipe(14).B,getpipe(6).A)
    K6=Joint("K6",pump.A,getpipe(4).B,getpipe(5).B,getpipe(6).B)
    
    Joint("ansl1", getpipe(1).B,FKG.A)
    Joint("ansl2", getpipe(4).A,FKG.B)
    Joint("ansl3", getpipe(7).B,FKA.A)
    Joint("ansl4", getpipe(9).A,FKA.B)
    Joint("ansl5", getpipe(8).B,FKB.A)
    Joint("ansl6", getpipe(10).A,FKB.B)
    Joint("ansl7", getpipe(11).B,FKV.A)
    Joint("ansl8", getpipe(13).A,FKV.B)
    Joint("ansl9", getpipe(12).B,FKM.A)
    Joint("ansl10",getpipe(14).A,FKM.B)
    
    K1.point.pressure=31e3
    
    

    """
    ror1=Pipe("ror1",

    till_kok=Pipe("till_kök",7.0,26.0e-3)

    till_sovrum=Pipe("till_sovrum",4.0,20.0e-3)
    
    fran_sovrum=Pipe("från_sovrum",4.0,20.0e-3)

    fran_kok=Pipe("från_kök",7.0,26.0e-3)
    
    
    #crc84=FanCoil("Crc84",1188.0,20.0)
    #crc54=FanCoil("Crc54",756.0,35.5)
    crc34=FanCoil("Crc34",468.0,13.1)
    crc34b=FanCoil("Crc34B",468.0,13.1)
    """
    
    
        
    
    #TPiece("T1",pipe1.A,pipe2.A,pump.A)
    #TPiece("T1",pipe1.B,pipe2.B,pump.B)
    
    """
    Joint("pumpansl.1",pump.A,till_kok.A)
    Joint("pumpansl.2",pump.B,fran_kok.B)
    
    TPiece("takskarv1",till_kok.B,crc84.A,till_sovrum.A)
    
    
    TPiece("takskarv2",fran_kok.A,crc84.B,fran_sovrum.B)
    
    Joint("sovansl.1",till_sovrum.B,crc34.A,crc34b.A)
    Joint("sovansl.2",crc34.B,fran_sovrum.A,crc34b.B)
    """
    
    
    
    
    simulator=Simulator()
    simulator.simulate(pipes+fancoils+[pump])

    drill_hole=200.0
    cooling_power=pump.flow*1000.0 * 4185 * (12.0-7.0)
    print "Pump flow:",pump.flow*1000.0*3600,"l/h"
    print "Cooling power:",cooling_power/1e3,"kW"
    print "Cooling power per m hole:",cooling_power/drill_hole,"W/m"
        
    for fk in fancoils:
        fk.printstate()        
if __name__=='__main__':
    run()
        