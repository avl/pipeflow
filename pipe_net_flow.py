#encoding=utf8
import math
import random
import odsreader
from pipe_utils import *
                
        
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
    FKV=crc84=FanCoil("FK-V",1188.0,20.0,6979.0,valve_kvs=5.5,stryp=100.0*0)
    FKA=crc24a=FanCoil("FK-A",324.0,15.9,1954.0,valve_kvs=1.7,stryp=100*1)
    FKB=crc24b=FanCoil("FK-B",324.0,15.9,1954.0,valve_kvs=1.7,stryp=100*1)
    FKG=crc54a=FanCoil("FK-G",756.0,35.5,4395.0,valve_kvs=2.8,stryp=100*(0))
    FKM=crc54b=FanCoil("FK-M",756.0,35.5,4395.0,valve_kvs=2.8,stryp=100*(1))
    fancoils=[FKV,FKA,FKB,FKG,FKM]

    
    
    #PCS44 pump: pump=Pump("pump",Curve(kpa_lph_to_si_units([(0.0,2800),(53.0,520),(60.0,0)])))
    #pump=Pump("Grundfos Alpha2",Curve(kpa_lph_to_si_units([(0.0,4200),(10.0,4000.0),(20.0,2900.0),(30.0,2300.0),(40.0,1600.0),(50.0,1100.0),(58.0,600.0),(60.0,0)])))
    pump=Pump("Linear pump",Curve(kpa_lph_to_si_units([(0.0,4200),(60.0,0)])))
    
    K1=Joint("K1",pump.B,getpipe(1).A,getpipe(2).A,getpipe(3).A)
    K2=Joint("K2",getpipe(2).B,getpipe(7).A,getpipe(8).A)
    K3=Joint("K3",getpipe(5).A,getpipe(9).B,getpipe(10).B)
    K4=Joint("K4",getpipe(3).B,getpipe(11).A,getpipe(12).A)
    K5=Joint("K5",getpipe(13).B,getpipe(14).B,getpipe(6).A)
    K6=Joint("K6",pump.A,getpipe(4).B,getpipe(5).B,getpipe(6).B)
    
    Joint("FKG-in", getpipe(1).B,FKG.A)
    Joint("FKG-ut", getpipe(4).A,FKG.B)
    Joint("FKA-in", getpipe(7).B,FKA.A)
    Joint("FKA-ut", getpipe(9).A,FKA.B)
    Joint("FKB-in", getpipe(8).B,FKB.A)
    Joint("FKB-ut", getpipe(10).A,FKB.B)
    Joint("FKV-in", getpipe(11).B,FKV.A)
    Joint("FKV-ut", getpipe(13).A,FKV.B)
    Joint("FKM-in", getpipe(12).B,FKM.A)
    Joint("FKM-ut",getpipe(14).A,FKM.B)
    
    K1.point.pressure=31e3
    
        
    
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
        
