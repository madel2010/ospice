#include <iostream>

#include "Circuit.h"
#include "subcircuit.h"
#include "resistor.h"
#include "Capacitor.h"
#include "Inductor.h"
#include "Source.h"
#include "SourceFunc.h"
#include "Analysis.h"
#include "Probes.h"

int main(int argc, char **argv) {
    
    Circuit C1;
    
    std::vector<std::string> term;
    term.push_back("s1");
    term.push_back("s2");
    SubCircuit S1("Xtest" , term);

    //C1<< new resistor("R1" , "n1" , "n2", 1);
    S1 << new resistor("R1" , "s1" , "n2", 1);
    
    //C1<< new Inductor("L1" , "n2" , "n3", 1);
    
    //C1<< new Capacitor ("C1" , "n3" , "0", 1);
    S1<< new Capacitor ("C1" , "n2" , "s2", 1);
    
    //C1<< new CurrentSource ("I1", "0", "n1" , new DCSource(1) ) ;
    //C1<< new CurrentSource ("I1", "0", "n1" , new PWLSource(10 , 0.0 , 0.0 , 1.0 , 1.0 ,2.0 , 1.0, 3.0, 0.0, 10.0, 0.0) ) ;
    C1<< new VoltageSource ("V1", "0", "n1" , new PWLSource(10 , 0.0 , 0.0 , 1.0 , 1.0 ,2.0 , 1.0, 3.0, 0.0, 10.0, 0.0) ) ;
    
    //C1<< new nonlin_resistor("R2" , "n2" , "0", "(10^(-12))*(exp(40*v(n2))-1)");
    
    std::vector<std::string> inst_term;
    inst_term.push_back("n1");
    inst_term.push_back("0");
    C1<< S1.create_instance("SI1" , inst_term);
    
    
    C1<< new DC;
    
    C1<< new  transient(0, 10 , 0.01);

    //C1 << new VoltageProbe("V(R2)" , "n1" , "0");
    
    C1.start_analysis();
    
    C1.plot_probes();
    
    return 0;
}
