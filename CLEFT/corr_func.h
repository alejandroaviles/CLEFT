#ifndef CORR_FUNC_H
#define CORR_FUNC_H

#include "q_depend_funcs.h"
//~ #include "lu_decomp.h"
#include "integral.h"
#include "prog_bar.h"

#include <vector>
#include <string>

struct corr_func_init
{
	double r_max, r_min;
	int r_bin_num;
	std::string file_name;
};

struct vec3
{
	double x, y, z;
};


class corr_func
{
	////////// Con-/destructor & initializer //////////
public:
	corr_func(  );
	~corr_func(  );
	void set_par( const corr_func_init & c_arg,
				  const q_func & q );
	
	////////// Integration kernel //////////
private:						// Data
	double b1, b2;
	q_func * qf;
private:						// Functions
	void M( const double & r, const vec3 & q );
	void M( const double & r, const double & q,
			  const double & mu );
	int delta_k( const int & i, const int & j );
	
	////////// Corr func xi //////////
private:						// Data
	double r_max, r_min;
	int r_bin_num;
	static const int num_bias_comp = 7;
	double bias_comp_inner[ num_bias_comp ];
	double bias_comp_outer[ num_bias_comp ];
	std::vector<double> r_buf;
	std::vector<double> xi_L_buf;
	std::vector<double> b10b20, b11b20, b10b21;
	std::vector<double> b12b20, b11b21, b10b22;
	std::vector<double> newb;

private:						// Functions
	void xi( const double & r );
	double xi_L( const double & r );
public:
	void get_xi(  );

	////////// Output //////////
private:
	std::string xi_file_name;
public:
	void output(  );
	
	////////// Progress bar //////////
private:
	prog_bar pg;
	
	////////// Mathematical const //////////
private:
	static const double pi = 3.14159265358979323846;
	static const double nearly_0 = 1e-5;
	static const double min_q_for_integration = 2e-1;
	static const double max_q_for_integration = 3e3;
};

#endif

