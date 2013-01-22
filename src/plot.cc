#include "plot.h"
#include <string.h>

#include <stdarg.h> //for the valist, va_start, va_end
#include <stdio.h> //for vasprintf
#include <stdlib.h>

#include <iostream>

#include <unistd.h>

int Plot::plot_fd;
bool Plot::ready;

Plot::Plot(){
	active = false;

	xlogscale = false;
	ylogscale =  false;
	x_label = NULL;
	y_label = NULL;

	if(!Plot::ready) init();
}

void Plot::init(){
	if(!Plot::ready){
		pid_t pid;
		int fd[2];
		pipe(fd);

		if( (pid=fork()) == -1){
			fprintf(stderr,"Fork error. Exiting.\n");  /* something went wrong */
		}

		if(pid){ //Parent
			Plot::plot_fd = fd[1];  //save the write end
			close(fd[0]); // close read end 

			//write commands to gnuplot
			write(plot_fd, "set terminal wxt persist enhanced font \"arial-12\"\n",4+9+4+8+9+5+11);
			
			Plot::ready = true;

		}else{ //Child
      			dup2(fd[0], 0); //duplicate read end to stdin
					//This is used that the external program will not write its output to the terminal 

      			//close them as we already duplicated them
			close(fd[0]); 
      			close(fd[1]);

      			execlp("gnuplot", "gnuplot", "-persist", NULL);

		}

	}
}

void Plot::send_command(const char * format , ...){
	if(!Plot::ready) init();
	
	char* command;

	//add the command to a string
	va_list args;
	va_start(args, format); //initialize args
	vasprintf(&command , format , args); /*coppy the data to a string
					     Note: using vasprintf takes care of allocating memory to the command pointer*/
	va_end(args);
	
	write(Plot::plot_fd , command, strlen(command));
		
	free(command); 
}

void Plot::update(){
	if (x_label){
    		send_command("set xlabel \"%s\"\n", x_label);
	}
  	if (y_label){
		send_command("set ylabel \"%s\"\n", y_label);
	}
	if (xlogscale){
		//send_command("%s \n","unset logscale x ");
		send_command("%s \n", "set logscale x");
	}
}

void Plot::draw_curves(){
	
	//send the plot command
	send_command("set autoscale y\n");
	send_command("set autoscale x\n");
	send_command("set xrange [*:*] noreverse\n");
	send_command("set yrange [*:*] noreverse\n");
	send_command("plot ");

	//enter the title for each curve
	std::vector<Plotlist*>::iterator it;
	for(int i=0; i<curves.size(); i++){
 		send_command("\'-\' title \"%s\" with lines lw 2", curves[i]->title.c_str());
		if(i<curves.size()-1){
			send_command(", ");
		}
	}

	send_command("\n");

	//start adding the data
	for(it=curves.begin(); it!=curves.end(); it++){
		int n = (*it)->n;
		for(int i=0; i<n; i++){
 			send_command("%e %e\n", (*it)->x[i], (*it)->y[i]);
		}
		send_command("e\n");
	} 	
}

void Plot::xlabel(const char* label){
	x_label = label;
}

void Plot::ylabel(const char* label){
	y_label = label;
}

void Plot::set_xlogscale(bool x){
	xlogscale = x;
}

void Plot::set_ylogscale(bool y){
	ylogscale = y;
}

void Plot::plot(std::string function){
	//update the plot for the data
	update();

	std::string cmd = "plot ";
	cmd+= function;
	cmd += "\n";
	send_command("%s",cmd.c_str());
}

void Plot::plot(std::vector<double> x, std::vector<double> y, std::string name){
	//update the plot for the data
	update();

	//add a curve in the Plotlist structure
	Plotlist* p = new Plotlist;

	int n = x.size(); //number of points

	p->title = name;
	p->n = n;
	p->x = new double[n];
 	p->y = new double[n];

	for(int i=0; i<n; i++){
		p->x[i] = x[i];
		p->y[i] = y[i];
	}
	
	//add the curve to the curves list
	curves.push_back(p);

	//draw the curves
	draw_curves();
}

void Plot::plot(int n_points, double* x, double* y, std::string name){
	//update the plot for the data
	update();

	//add a curve in the Plotlist structure
	Plotlist* p = new Plotlist;

	p->title = name;
	p->n = n_points;
	p->x = new double[n_points];
 	p->y = new double[n_points];

	for(int i=0; i<n_points; i++){
		p->x[i] = x[i];
		p->y[i] = y[i];
	}
	
	//add the curve to the curves list
	curves.push_back(p);

	//draw the curves
	draw_curves();
}

Plot::~Plot(){

	std::vector<Plotlist*>::iterator it;
	for(it=curves.begin(); it!=curves.end(); it++){
		delete (*it);	
	}
}