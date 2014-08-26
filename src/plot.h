#ifndef H_PLOT
#define H_PLOT

#include <vector>
#include <string>

struct Plotlist{
	int n;
	double* x;
	double* y;
	std::string title;

	Plotlist(){
		x = NULL;
		y = NULL;
		title = "";
	}
	
	~Plotlist(){
		if(x) delete[] x;
		if(y) delete[] y;
	}
};


class Plot{
	
private:
	bool xlogscale;
	bool ylogscale;
	const char* x_label;
	const char* y_label;
	bool active;
	

	std::vector<Plotlist*> curves;
	void draw_curves();

	void update();

public:
	static bool ready;
	static int plot_fd;

	Plot();
	
	void plot(int n_points, double* x, double* y, std::string name);
	void plot(std::vector<double> x, std::vector<double> y, std::string name);
	void plot(std::string function);  //if you want to just plot a function

	void init();
	void end();
	void send_command(const char * format , ...);
	void xlabel(const char* label);
	void ylabel(const char* label);
	void set_xlogscale(bool x);
	void set_ylogscale(bool y);
	~Plot();

};

#endif