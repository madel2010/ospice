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

class envelope_following;

class Analysis
{
protected:
    bool simulation_done; // This is a variable to check if we already did the analysis
    
public:
    Analysis(){};
    virtual ~Analysis(){};
    
    virtual void simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ)=0;
    virtual std::ostream& print(std::ostream &out) const =0;
    friend std::ostream& operator<< (std::ostream &out, const Analysis &B);
};


///This class is for the DC Analysis
class DC: public Analysis
{
 
private:  
    BMatrix::Dense<double> dc_solution; //the solution to the DC
 
    //newton Iteration with scaling B. Return NULL if not converged
    //B_scale -> for source stepping
    //Gmin_scale ->for stepping the conductances connected to the ground
    bool Newton_iter(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &Gmin, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, 
		     Circuit* circ, BMatrix::Dense<double> &solution, int B_scale=1, long int Gmin_scale=0);
public:
    DC();

    ~DC(){};
    void simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ);
     
    const BMatrix::Dense<double>& get_solution() { return dc_solution ; }
    
    //Print the output
    std::ostream& print(std::ostream &out)const;
};

///This class is for the transient Analysis
class transient: public Analysis
{
 friend envelope_following;
private:  
    BMatrix::Dense<double> tr_solution; //the solution to the DC
    
    bool set_initial_condition; //if we want to set an initial condition not to use the DC analysis
    BMatrix::Dense<double> my_initial_condition; //an initial condition for the circuit
        
    double start_time; //the starting time of the simulation_done
    double end_time; //the end time of the simulation
    double h; //the step size
    
    //time points to save solution for
    std::vector<double> save_solution_at;
    std::vector<BMatrix::Dense<double> > saved_solution;
    
    //to do the transient analysis. this is called from simulate(). it can also compute the Sensitivty Matrix
    bool perform_BE(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J,const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, 
			   BMatrix::Dense<double>& solution, double h, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix);
public:
    transient(double _start_time, double _end_time , double _h):start_time(_start_time), end_time(_end_time), h(_h){
      set_initial_condition = false;
    };

    void simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ);
    
    //this function if we want to return the Sensitivty Matrix
    void simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix);
    
    void set_save_solution_time(std::vector<double> _save_solution_at){save_solution_at = _save_solution_at;}

    std::vector<BMatrix::Dense<double> > get_saved_solution(){return saved_solution;};
    
    void use_initial_condition(BMatrix::Dense<double>& initial_condition){
	my_initial_condition = initial_condition;
	set_initial_condition = true;
    }
    
    //print the output
    std::ostream& print(std::ostream &out)const;
};


///This class is for Enveope Following Analaysis
class envelope_following: public Analysis
{
    
private:
      double start_time;
      double end_time;
      double H; //The jump
      double h; //the small step size for normal transient anlysis
      double T; //The period of high frequency component
    
      transient * my_transient; //the transient analysis that will be used

      
public:
     envelope_following(double _start_time, double _end_time , double _H, double _h, double _T):start_time(_start_time), end_time(_end_time), H(_H), h(_h), T(_T){
        my_transient = new  transient(start_time, start_time+T , h);

    };
     
     envelope_following(double _start_time, double _end_time , double _H, double _T):start_time(_start_time), end_time(_end_time), H(_H), T(_T){
       h = T/100;
       my_transient = new  transient(start_time, start_time+T , h);

    };

     void simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ);

     std::ostream& print(std::ostream &out) const;
};

#endif // ANALYSIS_H
