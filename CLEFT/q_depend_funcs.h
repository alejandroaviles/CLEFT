#ifndef Q_DEPEND_FUNCS_H
#define Q_DEPEND_FUNCS_H

#include "k_depend_funcs.h"
#include "integral.h"
#include "prog_bar.h"
#include <vector>
#include <string>
#include <map>

// All q-dependent functions in Jordan Carlson's notes

struct q_func_init
{
	std::string pow_spec_name;
	std::string k_file_name;
	std::string q_file_name;
};

struct q_func_vals
{
	double xi_L;
	
   double Xi0Lin, Xi0loop, Xi0eft; 
   double Xi1loop, Xi1eft; 
   double Xi2Lin, Xi2loop, Xi2eft; 
   double Xi3loop, Xi3eft;
   double Xi0LinZeroLag, Xi0loopZeroLag;
            
};


class q_func
{
	////////// Datatype def and declaration //////////
private:
	typedef std::vector<double> dvec;
	typedef std::map<double, double>::iterator itr_PL;
	
	////////// Con-/destructor & initializer //////////
public:
	q_func(  );
	~q_func(  );
	void set_par( const q_func_init & arg );

	////////// Do all preparations //////////
private:
	void cal_all( std::string pow_spec_name );
	void load_k( std::string pow_spec_name,
				 std::string k_file_name );
	void load_all( std::string pow_spec_name,
				   std::string k_file_name,
				   std::string q_file_name );

	////////// Functions of k //////////
public:
	const k_func & kfunc(  );
private:
	k_func kf;
	
	////////// Spherical Bessel functions //////////
private:
	double sph_bessel_j( int n, const double & x );
	// n for order.

	////////// Various functions //////////
public:
	const std::vector<double> & qvec(  );
private:						// Data
	std::map<double, int> q_index_buf;
	std::vector<double> q_buf;
	std::vector<double> k_intg_buf;
	std::vector<double> xi_L_buf;

	
	std::vector<double> Xi0Lin_buf, Xi0loop_buf, Xi0eft_buf; 
    std::vector<double> Xi1loop_buf, Xi1eft_buf;
    std::vector<double> Xi2Lin_buf, Xi2loop_buf, Xi2eft_buf;
    std::vector<double> Xi3loop_buf, Xi3eft_buf;
    std::vector<double> Xi0LinZeroLag_buf, Xi0loopZeroLag_buf;
	
	
private:						// Function
	void get_var_func(  );
	double xi_L_intg( const double & q );
	
	double Xi0Lin_intg( const double & q ), Xi0loop_intg( const double & q ), Xi0eft_intg( const double & q ); 
    double Xi1loop_intg( const double & q ), Xi1eft_intg( const double & q ); 
    double Xi2Lin_intg( const double & q ), Xi2loop_intg( const double & q ), Xi2eft_intg( const double & q ); 
    double Xi3loop_intg( const double & q ), Xi3eft_intg( const double & q );	
    double Xi0LinZeroLag_intg( const double & q ), Xi0loopZeroLag_intg( const double & q );
	

	
	// Linear interpolation 
	double interp_val( const double & q, const int & i,
					   const std::vector<double> & vec );
public:
	void var_func( const double & q, q_func_vals & res );

	////////// Save and load //////////
private:
	void save_q_func( std::string file_name );
	void load_q_func( std::string file_name );

	////////// Progress bar //////////
private:
	prog_bar pg;
	
	////////// Mathematical constants/func //////////
private:
	static const double nearly_0 = 1.e-3;
	static const double nearly_inf = 2.e2;
	static const double pi = 3.14159265358979323846;
	static const double one_over_pi2 = 0.05066059182116889;  //1/(2 Pi^2)
	static const int k_intg_points_multip = 3;
	integral intg;
};

#endif

