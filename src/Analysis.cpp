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
#include <string>


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
      perform_BE(G, C, J, B, fx, circ, NULL);

}

void transient::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Sparse<double> &J, BMatrix::Dense<double> &B, BMatrix::Dense<double> &fx, Circuit* circ, BMatrix::Sparse< double >& _Sensitivty_Matrix){
      BMatrix::Sparse< double >* Sensitivty_Matrix = &_Sensitivty_Matrix;
      perform_BE(G, C, J, B, fx, circ, Sensitivty_Matrix);
}

void transient::perform_BE(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Sparse<double> &J, BMatrix::Dense<double> &B, BMatrix::Dense<double> &fx, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix){
     BMatrix::Sparse<double> scaled_C = C/h; // let us save the C/h because we have constant step size
     BMatrix::Sparse<double> temp, G_p_C_p_J;
    
    if(circ->is_circ_linear()){
	double time = start_time+h;
	
	BMatrix::Sparse<double> G_p_C = G+scaled_C; // G_p_C = G+C/h
    
        if(!set_initial_condition){
	  tr_solution = circ->get_dc_solution();
	}else{
	  tr_solution = my_initial_condition;
	}
	
	while(time <= end_time){

	    circ->update_sources(time);
	    
	    tr_solution = ( scaled_C*tr_solution ) + B;
	    
	    G_p_C.solve(tr_solution);
	    
	    //calculate sensitivity matrix
	    if(Sensitivty_Matrix){
	        //sensetivity_matrix = (J\(C/h)) * sensetivity_matrix;
	        temp  = scaled_C;
	        G_p_C.solve(temp);
		(*Sensitivty_Matrix) = temp*(*Sensitivty_Matrix);
	    }
	    
	    circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value
	    time+= h;
	}
    }else{
	  double time = start_time+h;
	
	  if(!set_initial_condition){
	    tr_solution = circ->get_dc_solution();
	  }else{
	    tr_solution = my_initial_condition;
	  }
	
	  BMatrix::Dense<double> Phi(circ->size_of_mna(),1), pre_tr_solution(circ->size_of_mna(),1);
	  
	  while(time <= end_time){
  
	      std::cout<<"Current time = "<<time<<std::endl;
	      
	      circ->update_sources(time);
	    
	      BMatrix::Sparse<double> G_p_C = G+scaled_C; // let us save the C/h because we have constant step size
    
	      pre_tr_solution = tr_solution;
	      
	      bool convergence=false;
	      int number_of_iterations = 0;;
	      while(!convergence){
	    
		  circ->update_fx((*tr_solution));
		  circ->update_J((*tr_solution));
		  
		  Phi = G_p_C*tr_solution + fx - B - ( scaled_C*pre_tr_solution );

		  if(Phi.norm()<=1e-12 || (number_of_iterations>15) ){
		      convergence = true;
		      
		  //calculate sensitivity matrix
		  if(Sensitivty_Matrix){
		      //sensetivity_matrix = (J\(C/h)) * sensetivity_matrix;
		      temp  = scaled_C;
		      G_p_C_p_J.solve(temp);
		      (*Sensitivty_Matrix) = temp*(*Sensitivty_Matrix);
		     
		  }
	    
		  }else{
		      G_p_C_p_J = G_p_C+J;
		      G_p_C_p_J.solve(Phi);  //Note: solve function rewrites the Phi

		      tr_solution -= Phi;
     
		}
		number_of_iterations++;
	      }

	      circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value

	      if(std::find(save_solution_at.begin(),save_solution_at.end(),time)!=save_solution_at.end()){
		  saved_solution.push_back(tr_solution);
	      }
	      time+= h;
	}
	
    }
}

/*************Envelope Following Class -------------------------*/
void envelope_following::simulate(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, BMatrix::Sparse<double> &J, BMatrix::Dense<double> &B, BMatrix::Dense<double> &fx, Circuit* circ){
  int MNA_size = circ->size_of_mna();
  
  my_transient->save_solution_at.push_back(start_time);
  my_transient->save_solution_at.push_back(start_time+T);
  my_transient->simulate(G, C, J, B, fx,circ);
  
  circ->update_probes(start_time , (*my_transient->saved_solution[0]) ); //This function should send to all the probes in the circuit the solution to add its value
  circ->update_probes(start_time+T , (*my_transient->saved_solution[1]) ); //This function should send to all the probes in the circuit the solution to add its value

 double my_time = start_time+H;
 BMatrix::Dense<double> dym(MNA_size,1), Yn(MNA_size,1), Yn1(MNA_size,1); 
 BMatrix::Sparse<double> dyn1dyn(MNA_size,MNA_size);
 
 BMatrix::Sparse<double> I = EYE(MNA_size);

 BMatrix::Dense<double> Phi(MNA_size,1);

 while(my_time <= end_time){
	  dym = (my_transient->saved_solution[1] - my_transient->saved_solution[0])/T; 
	  
	  //Initial Guess
	  Yn = my_transient->saved_solution[0];
	  
	  std::cout<<"Current time = "<<my_time<<std::endl;
	      	    	      
	  bool convergence=false;
	  int number_of_iterations = 0;;
	  while(!convergence){

	    
	          //Get Yn1
	          my_transient->use_initial_condition(Yn);
		  
		  my_transient->start_time = my_time;
		  my_transient->end_time = my_time+T;
		  
		  my_transient->save_solution_at.clear();
		  my_transient->save_solution_at.push_back(my_time+T);
		  
		  my_transient->simulate(G, C, J, B, fx,circ,dyn1dyn);
		  
		  
		  Yn1 = my_transient->saved_solution[0];
 		  Phi = Yn - my_transient->saved_solution[0] - (Yn1-Yn + my_transient->saved_solution[1] - my_transient->saved_solution[0])*H/(2*T);

		  if(Phi.norm()<=1e-6 || (number_of_iterations>10) ){
		      convergence = true;
		  }else{
	   
		      (I - (dyn1dyn-I)*H/(2*T)).solve(Phi);  //Note: solve function rewrites the Phi

		      Yn -= Phi;
     
		  }
		number_of_iterations++;
	  }
	  
	  circ->update_probes(my_time , (*Yn) );
	  circ->update_probes(my_time+T , (*Yn1) );

	  
 }
  
}
