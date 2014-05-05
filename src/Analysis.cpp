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
}



void DC::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Sparse<double> &J, BMatrix::Dense<double> &B, BMatrix::Dense<double> &fx, Circuit* circ){

    if(circ->is_circ_linear()){
      
	circ->update_sources(0);
	
	dc_solution = B;
	G.solve(dc_solution);  //Note: solve function rewrites the dc_solution
	
	//Now update the probes of the circuit with the solution
	circ->update_probes( 0 , (*dc_solution) ); //this function takes the time and the solution vector
	
    }else{
	dc_solution.create(circ->size_of_mna(),1);
	dc_solution = 0;
	    
	BMatrix::Dense<double> Phi = BMatrix::Dense<double>(circ->size_of_mna(),1);
	
	bool convergence=false;
	
	while(!convergence){
	    
	    circ->update_fx((*dc_solution));
	    circ->update_J((*dc_solution));
	    
	    //G here is actually G+J (check  circ->update_G_P_J())
	    Phi = G*dc_solution + fx - B;

	    if(Phi.norm()<=1e-12){
		convergence = true;
	    }else{
	   
	    	//G here is actually G+J (check  circ->update_G_P_J())
		(G+J).solve(Phi);  //Note: solve function rewrites the Phi
		dc_solution -= Phi;
	    }
	}
	circ->update_probes(0 , (*dc_solution) ); //This function should send to all the probes in the circuit the solution to add its value

    }
     
}

/*--------------------The transient analysis------------*/
void transient::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Sparse<double> &J, BMatrix::Dense<double> &B, BMatrix::Dense<double> &fx, Circuit* circ){

    BMatrix::Sparse<double> scaled_C = C/h; // let us save the C/h because we have constant step size
    
    if(circ->is_circ_linear()){
	double time = start_time+h;
	
	BMatrix::Sparse<double> G_p_C = G+C/h; // let us save the C/h because we have constant step size
    
	tr_solution = circ->get_dc_solution();
	
	while(time <= end_time){

	    circ->update_sources(time);
	    
	    tr_solution = ( scaled_C*tr_solution ) + B;
	    
	    G_p_C.solve(tr_solution);
	    
	    circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value
	    time+= h;
	}
    }else{
	  double time = start_time+h;
	
	  tr_solution = circ->get_dc_solution();
	
	  BMatrix::Dense<double> Phi(circ->size_of_mna(),1), pre_tr_solution(circ->size_of_mna(),1);
	  
	  while(time <= end_time){
  
	      std::cout<<"Current time = "<<time<<std::endl;
	      
	      circ->update_sources(time);
	    
	      BMatrix::Sparse<double> G_p_C = G+C/h; // let us save the C/h because we have constant step size
    
	      Phi = BMatrix::Dense<double>(circ->size_of_mna(),1);

	      pre_tr_solution = tr_solution;
	      
	      bool convergence=false;
	      int number_of_iterations = 0;;
	      while(!convergence){
	    
		  circ->update_fx((*tr_solution));
		  circ->update_J((*tr_solution));
		  
		  //G here is actually G+J (check  circ->update_G_P_J())
		  Phi = G_p_C*tr_solution + fx - B - ( scaled_C*pre_tr_solution );

		  if(Phi.norm()<=1e-12 || (number_of_iterations>15) ){
		      convergence = true;
		  }else{
	   
		      //G here is actually G+J (check  circ->update_G_P_J())
		      (G_p_C+J).solve(Phi);  //Note: solve function rewrites the Phi

		      tr_solution -= Phi;
     
		}
		number_of_iterations++;
	      }

	      circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value
	      time+= h;
	}
	
    }
     
}
