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


#include "SourceFunc.h"
#include <stdexcept>
#include <stdarg.h>
#include <math.h>
#include "config.h"

/*---------------------THE PWL FUNCTIONS-------------*/
PWLSource::PWLSource(int arg_number, ...){
  
  //check that we have even number of arguments
  if(arg_number%2 != 0){
    throw std::runtime_error("PWL linear function has wrong number of arguments");
  }
  number_of_points = arg_number/2;
  
  values = new double [number_of_points];
  time_points = new double[number_of_points];
  
  va_list Args; //declasre the iterator through the arguments
  va_start(Args, arg_number); 

  for(int i = 0; i < number_of_points; i++ ){
        time_points[i] = va_arg(Args, double);
	values[i] = va_arg(Args, double);
  }

  va_end(Args); //end iterating through the arguments list
  
  
}

PWLSource::~PWLSource(){
  
  if(values) delete[] values;
  if(time_points) delete[] time_points;
  
}
    
double PWLSource::get_value(double _time){
    bool found_value = false;
    double result = 0;
    
    for (int i=1; i<number_of_points; i++){
	if(_time >= time_points[i-1] && _time < time_points[i]){
    
	    //Linear interpolate the value
	    result = values[i-1] + (values[i]-values[i-1])*(_time-time_points[i-1])/(time_points[i]-time_points[i-1]);
	
	    found_value = true;
	    break;
	}
    }
    
    //if we can not find a values, then just use the last value (we assume that it will be constant for time=infinty)
    if(!found_value){
	result = values[number_of_points];
    }
    
    return result;
}

/*---------------------THE PWL FUNCTIONS-------------*/
SinSource::SinSource(double _Vo, double _Va, double _Freq, double _Td , double _Df , double _Phase){
    Vo = _Vo;
    Va = _Va;
    Freq = _Freq;
    Td = _Td;
    Df = _Df;
    Phase = _Phase;
}
 
    
double SinSource::get_value(double time){
  return Vo+Va*sin(2*PI*(Freq*(time-Td)+Phase/360))*exp(-(time-Td)*Df);
}

