#include <iostream>

#include "Circuit.h"
#include "resistor.h"
#include "Capacitor.h"
#include "Inductor.h"
#include "Source.h"
#include "SourceFunc.h"
#include "Analysis.h"
#include "Probes.h"

int main(int argc, char **argv) {
    
    Circuit C1;
    
    C1<< new resistor("R1" , "n1" , "n2", 1);
    
    C1<< new nonlin_resistor("R1" , "n1" , "n2", "(10^(-2))*(exp((v(n1)-v(n2))/25)-1)");
    
    C1<< new Inductor("L1" , "n2" , "n3", 1);
    
    C1<< new Capacitor ("C1" , "n3" , "0", 1);
    
    //C1<< new CurrentSource ("I1", "0", "n1" , new DCSource(1) ) ;
    //C1<< new CurrentSource ("I1", "0", "n1" , new PWLSource(10 , 0.0 , 0.0 , 1.0 , 1.0 ,2.0 , 1.0, 3.0, 0.0, 10.0, 0.0) ) ;
    C1<< new VoltageSource ("V1", "0", "n1" , new PWLSource(10 , 0.0 , 0.0 , 1.0 , 1.0 ,2.0 , 1.0, 3.0, 0.0, 10.0, 0.0) ) ;
    
    C1<< new DC;
    
    C1<< new  transient(0, 10 , 0.1);
    
    C1 << new VoltageProbe("V(C1)" , "n3" , "0");
    
    C1.start_analysis();
    
    C1.plot_probes();
    
    return 0;
}
