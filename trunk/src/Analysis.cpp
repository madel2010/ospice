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
#include <math.h>
#include "debug.h"
#include "util.h"

BMatrix::Sparse<double> DC::Jac = BMatrix::Sparse<double>();
BMatrix::Sparse<double> transient::G_p_C_p_J = BMatrix::Sparse<double>();

///THE DC ANALYSIS
DC::DC(){
    simulation_done = false;
    Jac.Freeze_structure = true; //not do sparse ordering every time
}

bool DC::Newton_iter(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &Gmin, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, 
		     Circuit* circ, BMatrix::Dense<double> &solution, int B_scale, long int Gmin_scale){

    
    BMatrix::Dense<double> Phi;
    
    int Iter_Number = 0;
    double Phi_norm=0 , Phi_norm_old=0;
    bool convergence=false;	
    while(!convergence && !std::isnan(Phi_norm) && (Iter_Number < 400 || Phi_norm<Phi_norm_old)){
	    
	circ->update_fx((*solution));
	circ->update_J((*solution));
	    
	if(Gmin_scale!=0.0){ //do the Gmin algorithm
	  BMatrix::Sparse<double> A = G+Gmin/Gmin_scale;
	  Phi = A*solution + fx - B;
	  Jac = A+J;
	}else if(B_scale!=1){
	  Phi = G*solution + fx - B/B_scale;
	  Jac = G+J;
	}else{
	  Phi = G*solution + fx - B;
	  Jac = G+J;
	} 
		
	Phi_norm_old = Phi_norm;
	Phi_norm = Phi.norm();
	
	if(Phi_norm<=1e-6){
		convergence = true;
	}else{
		Jac.solve(Phi);  //Note: solve function rewrites the Phi
		solution -= Phi;
                Iter_Number++;
	}

	_DD(2){
	   if(convergence) std::cout<<"DC Analysis: Phi = " << Phi_norm <<std::endl; 
	}

     }

     return convergence;
}

void DC::simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ){

     circ->update_sources(0);
    
    
     dc_solution.create(circ->size_of_mna(),1);
     dc_solution.reset();

     BMatrix::Sparse<double> Gmin(circ->size_of_mna() , circ->size_of_mna());
     
     //First try with B_scale=1 Gmin_scale=0
     bool Newton_iter_result =  Newton_iter(G, Gmin, C, J, B,fx, circ, dc_solution);
     
     
     if(!Newton_iter_result){ //No convergence 
	//Try G scaling
        _DD(1){ std::cout<<"Normal DC did not converge, trying Gmin steping"<<std::endl;}

        for (auto n: circ->get_node_indeces()){
	    Gmin.put(n,n, 1e14);
	}
 	dc_solution.reset();
	for (long int i=1; i<=1e23; i=i*10){
    		 Newton_iter_result = Newton_iter(G, Gmin, C, J, B, fx, circ, dc_solution, 1 , i);

		 if(!Newton_iter_result){ //One step did not converge 
			//exit 
                        break;
		 }
        }
        
        if(Newton_iter_result){
	    //we need to do a final iteration without Gmin
	    Newton_iter_result = Newton_iter(G, Gmin, C, J, B, fx, circ, dc_solution);
	}
    
     }
     
     if(!Newton_iter_result){ //No convergence
        //Try source stepping
 	_DD(1){ std::cout<<"Normal DC did not converge, trying DC source steping"<<std::endl;}
 	
 	dc_solution.reset();
	for (int i=1000; i>=1; i--){
    		 Newton_iter_result = Newton_iter(G, Gmin, C, J, B, fx, circ, dc_solution, i , 0);

		 if(!Newton_iter_result){ //One step did not converge 
			//exit with Error
			throw std::runtime_error(std::string("DC did not converge using source steping"));
                        break;
		 }
        }
     }
     circ->update_probes(0 , (*dc_solution) ); //This function should send to all the probes in the circuit the solution to add its value

     _DD(1){
	 std::cout<<"DC Analysis Report:"<<std::endl;
	 Jac.report_timing();
      }
}

std::ostream& DC::print(std::ostream &out , const Circuit* circ)const{
    std::map<std::string , int>::const_iterator MNA_iter = circ->get_mna_variable_iterator();
    std::string name; // whether it V, I or Q
    for(int i=0; i<circ->size_of_mna(); i++){
	
	if(MNA_iter->first.substr(MNA_iter->first.size() - 1) == "I"){
		name = "I(" + MNA_iter->first.substr(0,MNA_iter->first.size() - 2) + ")";
	}else if(MNA_iter->first.substr(MNA_iter->first.size() - 1) == "Q"){
		name = "Q(" + MNA_iter->first.substr(0,MNA_iter->first.size() - 2) + ")";
	}else{
		name = "V(" + MNA_iter->first + ")";
	}
	out<< name << "=" << dc_solution.get(MNA_iter->second,0)<<std::endl;
	
	MNA_iter++;
    }
    
    return out;
}


/*--------------------The transient analysis------------*/
void transient::simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ){
      simulate(G,C, J,B, fx,  circ, nullptr);
}

void transient::simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix){
     
      G_p_C_p_J.Freeze_structure = true; //do not sparse ordering every time

      BMatrix::Dense<double> TR_B;

      if(!set_initial_condition){
	   tr_solution = circ->get_dc_solution();
      }else{
	   tr_solution = my_initial_condition;
      }
	  
      double time = start_time+h;
      int s = 1; //to find the first time point and use BE
      bool convergence;
      while(time <= end_time){	  
	std::cout<<"Current time = "<<time;
	
	TR_B = B; //B at previous time point	
	
	if(circ->config.integration_order==1 || s==1){
		circ->update_sources(time);
		convergence = perform_BE(G, C, J, B, fx, tr_solution ,h, circ, Sensitivty_Matrix);
	}else if(circ->config.integration_order==2){
		circ->update_sources(time);
		TR_B += B;
		convergence = perform_TR(G, C, J, TR_B, fx, tr_solution ,h, circ, Sensitivty_Matrix);
	}else{
		throw std::runtime_error("Transient analysis: Integration order not known");
	}

	if(!convergence){
		throw std::runtime_error("Transient analysis: No convergence");
	}

	circ->update_probes(time , (*tr_solution) ); //This function should send to all the probes in the circuit the solution to add its value

	if(std::find(save_solution_at.begin(),save_solution_at.end(),time)!=save_solution_at.end()){
	    saved_solution.push_back(tr_solution);
	}
	 
	time+= h;
	s++;
      }

      _DD(1){
	 std::cout<<"Transient Analysis Report:"<<std::endl;
	 G_p_C_p_J.report_timing();
      }
}

//solution should be the previous point and is overwritten by the new solution
bool transient::perform_BE(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, 
			   BMatrix::Dense<double>& solution, double h, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix){
  
     int MNA_size = circ->size_of_mna();
     
     BMatrix::Sparse<double> scaled_C = C/h; // let us save the C/h because we have constant step size
    
     BMatrix::Dense<double> Phi(MNA_size,1);
     BMatrix::Dense<double> scaled_C_times_pre_solution = scaled_C*solution;

     BMatrix::Sparse<double> temp;	  

     BMatrix::Sparse<double> G_p_C = G+scaled_C; // let us save the C/h because we have constant step size
    
     
      
     bool convergence=false;
     int number_of_iterations = 0;

     double max_phi;

     while(!convergence){
	 
        //if the number of iterations exeeded a certain value, then stop
        if(number_of_iterations>30){
	    convergence=false;
	    break;
	}
	
	circ->update_fx((*solution));
	circ->update_J((*solution));
		  
	Phi = G_p_C*solution + fx - B - ( scaled_C_times_pre_solution );
	
	G_p_C_p_J = G_p_C+J;
	G_p_C_p_J.solve(Phi);  //Note: solve function rewrites the Phi

	solution -= Phi;

	number_of_iterations++;
	
	max_phi = max_abs(Phi);
	if( max_phi < max_abs(solution*circ->config.abstol + circ->config.abstol ) ){
	      convergence = true;
		      
	      //calculate sensitivity matrix
	      if(Sensitivty_Matrix){
		      //sensetivity_matrix = (J\(C/h)) * sensetivity_matrix;
		      temp  = scaled_C;
		      G_p_C_p_J.solve(temp);
		      (*Sensitivty_Matrix) = temp*(*Sensitivty_Matrix);
	      }

	      break;
	    
        }
	
    }

    std::cout<<" PHI="<<max_phi<<std::endl;
    return convergence;
    
}

//solution should be the previous point and is overwritten by the new solution
bool transient::perform_TR(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, 
			   BMatrix::Dense<double>& solution, double h, Circuit* circ, BMatrix::Sparse< double >* Sensitivty_Matrix){
  
     int MNA_size = circ->size_of_mna();
     
     BMatrix::Sparse<double> scaled_C = C*(2/h); // let us save the C/h because we have constant step size
    
     BMatrix::Dense<double> Phi(MNA_size,1);
     BMatrix::Dense<double> scaled_C_times_pre_solution = (G-scaled_C)*solution;

     BMatrix::Sparse<double> temp;	  

     BMatrix::Sparse<double> G_p_C = G+scaled_C; // let us save the C/h because we have constant step size
    
     
      
     bool convergence=false;
     int number_of_iterations = 0;

     //get fx_pre at previous time point
     circ->update_fx((*solution));
     BMatrix::Dense<double> fx_pre = fx;

     double max_phi;

     while(!convergence){
	 
        //if the number of iterations exeeded a certain value, then stop
        if(number_of_iterations>30){
	    convergence=false;
	    break;
	}
	
	circ->update_fx((*solution));
	circ->update_J((*solution));
		  
	//NOTE: B = B - B_pre ... check simulate()
	Phi = G_p_C*solution + fx + fx_pre - B + ( scaled_C_times_pre_solution );

	G_p_C_p_J = G_p_C+J;
        G_p_C_p_J.solve(Phi);  //Note: solve function rewrites the Phi

        solution -= Phi;

	number_of_iterations++;

	max_phi = max_abs(Phi);
	if(max_phi < max_abs(solution*circ->config.abstol + circ->config.abstol) ){
	      convergence = true;
		      
	      //calculate sensitivity matrix
	      if(Sensitivty_Matrix){
		      //sensetivity_matrix = (J\(C/h)) * sensetivity_matrix;
		      temp  = scaled_C;
		      G_p_C_p_J.solve(temp);
		      (*Sensitivty_Matrix) = temp*(*Sensitivty_Matrix);
	      }
	    
        }

	
    }

    std::cout<<" PHI="<<max_phi<<std::endl;
    return convergence;
    
}

std::ostream& transient::print(std::ostream &out , const Circuit* circ)const{
    
    out<<"Transient is done"<<std::endl;
    
    return out;
}

/*************Envelope Following Class -------------------------*/
void envelope_following::simulate(const BMatrix::Sparse<double> &G, const BMatrix::Sparse<double> &C, const BMatrix::Sparse<double> &J, const BMatrix::Dense<double> &B, const BMatrix::Dense<double> &fx, Circuit* circ){
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
		  
		  my_transient->simulate(G, C, J, B, fx,circ,&dyn1dyn);
		  
		  
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

std::ostream& envelope_following::print(std::ostream &out , const Circuit* circ)const{
    
    out<<"EF is done"<<std::endl;
    
    return out;
}


