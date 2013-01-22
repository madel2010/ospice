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


#include "Analysis.h"


///THE DC ANALYSIS
DC::DC(){
    simulation_done = false;
    dc_solution = NULL;
}



void DC::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Dense<double> &B, Circuit* circ){

    if(circ->is_circ_linear()){
      
	circ->update_sources(0);
	
	dc_solution = B;
	G.solve(dc_solution);  //Note: solve function rewrites the dc_solution
	
	//Now update the probes of the circuit with the solution
	circ->update_probes( 0 , (*dc_solution) ); //this function takes the time and the solution vector
	
    }
     
}

/*--------------------The transient analysis------------*/
void transient::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Dense<double> &B, Circuit* circ){

    BMatrix::Sparse<double> scaled_C = C/h; // let us save the C/h because we have constant step size
    BMatrix::Sparse<double> G_p_C = G+C/h; // let us save the C/h because we have constant step size
    
    if(circ->is_circ_linear()){
	double time = start_time+h;
	
	tr_solution = circ->get_dc_solution();
	
	while(time <= end_time){

	    circ->update_sources(time);
	    
	    tr_solution = ( scaled_C*tr_solution ) + B;
	    std::cout<<(G_p_C);
	    G_p_C.solve(tr_solution);
	    
	    circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value
	    time+= h;
	}
    }
     
}
