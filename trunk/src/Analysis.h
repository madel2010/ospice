/*
    Copyright (c) 2012, <copyright holder> <email>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY <copyright holder> <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Circuit.h"
#include "BMatrix.h"

class Analysis
{
protected:
    bool simulation_done; // This is a variable to check if we already did the analysis
    
public:
    Analysis(){};
    virtual ~Analysis(){};
    
    virtual void simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Dense<double> &B, Circuit* circ)=0;
};


///This class is for the DC Analysis
class DC: public Analysis
{
 
private:  
    BMatrix::Dense<double> dc_solution; //the solution to the DC

public:
    DC();

    ~DC(){};
    void simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Dense<double> &B, Circuit* circ);
     
     const BMatrix::Dense<double>& get_solution() { return dc_solution ; }
};

///This class is for the DC Analysis
class transient: public Analysis
{
 
private:  
    BMatrix::Dense<double> tr_solution; //the solution to the DC
    double h; //the step size
    
    double start_time; //the starting time of the simulation_done
    double end_time; //the end time of the simulation
    
public:
    transient(double _start_time, double _end_time , double _h):start_time(_start_time), end_time(_end_time), h(_h){};

    void simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Dense<double> &B, Circuit* circ);
    
};

#endif // ANALYSIS_H
