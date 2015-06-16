#include "q_depend_funcs.h"
#include "prog_bar.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

////////////////////////////////////////////////////////////
// Constructor, desctructor and initializer

q_func::q_func(  )
{

}

q_func::~q_func(  )
{

}

void q_func::set_par( const q_func_init & arg )
{
    if( arg.k_file_name == "none" )
	cal_all( arg.pow_spec_name );
    else if( arg.q_file_name == "none" )
	load_k( arg.pow_spec_name, arg.k_file_name );
    else
	load_all( arg.pow_spec_name,
		  arg.k_file_name,
		  arg.q_file_name );
    return;
}


////////////////////////////////////////////////////////////
// To cal(culate) or not to cal, it is a problem.

void q_func::cal_all( std::string pow_spec_name )
{
    kf.load_PL( pow_spec_name );
    kf.get_Q_func(  );
    kf.get_R_func(  );
    get_var_func(  );
    kf.save_k_func( "../data/k_func.dat" );
    save_q_func( "../data/q_func.dat" );
    return;
}

void q_func::load_k( std::string pow_spec_name,
		     std::string k_file_name )
{
    kf.load_PL( pow_spec_name );
    kf.load_k_func( k_file_name );
    get_var_func(  );
    save_q_func( "../data/q_func.dat" );
    return;
}

void q_func::load_all( std::string pow_spec_name,
		       std::string k_file_name,
		       std::string q_file_name )
{
    kf.load_PL( pow_spec_name );
    kf.load_k_func( k_file_name );
    load_q_func( q_file_name );
    return;
}

////////////////////////////////////////////////////////////
// Functions of k

const k_func & q_func::kfunc(  )
{
    return kf;
}

////////////////////////////////////////////////////////////
// Spherical Bessel functions

double q_func::sph_bessel_j( int n, const double & x )
{
    if( fabs( x ) < nearly_0 )
	switch( n )
	{
	case 0:
	    return 1. - x*x / 6.;
	case 1:
	    return x / 3. - pow( x, 3 ) / 30.;
	case 2:
	    return x*x / 15.;
	case 3:
	    return pow( x, 3 ) / 105.;
	default:
	    throw "Bessel function error.";
	}
    switch( n )
    {
    case 0:
	return sin( x ) / x;
    case 1:
	return sin( x ) / pow( x, 2 )
	    - cos( x ) / x;
    case 2:
	return ( 3. / pow( x, 2 ) - 1 ) * sin( x ) / x
	    - 3 * cos( x ) / pow( x, 2 );
    case 3:
	return ( 15. / pow( x, 3 ) - 6. / x ) * sin( x ) / x
	    - ( 15. / pow( x, 2 ) - 1. ) * cos( x ) / x;
    default:
	throw "Bessel function error.";
    }
    return 0.;
}

////////////////////////////////////////////////////////////
// Various functions

const std::vector<double> & q_func::qvec(  )
{
    return q_buf;
}

void q_func::get_var_func(  )
{
    double q( 0. ), k( 0. );
    std::pair<double, int> q_index_pair;

    const std::vector<double> & kv = kf.kvec(  );

    const double k_max = kv[ kv.size(  ) - 1 ]/2.;
    const double k_min = kv[ 0 ];
    const double q_max = 1. / kv[ 0 ];
    const double q_min = 1. / kv[ kv.size(  ) - 1 ];

    for( unsigned i = 0; i < kv.size(  ); ++ i )
    {
	q = log( q_min )
	    + ( log( q_max ) - log( q_min ) ) * i
	    / double( kv.size(  ) - 1 );
	q_buf.push_back( exp( q ) );
    }

    const unsigned n_k_intg
	= k_intg_points_multip * kv.size(  );
    const double d_logk
	= ( log( k_max ) - log( k_min ) )
	/ double( n_k_intg - 1 );
    for( unsigned i = 0; i < n_k_intg; ++ i )
    {
	k = exp( log( k_min ) + i * d_logk );
	k_intg_buf.push_back( k );
    }
    std::cout << "Generating q-dependent functions: ";
    pg.init( q_buf.size(  ) );

    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	pg.show( i );
	q = q_buf[ i ];
	q_index_pair.first = q;
	q_index_pair.second = i;
	q_index_buf.insert( q_index_pair );

	xi_L_buf.push_back( xi_L_intg( q ) );


	
   Xi0Lin_buf.push_back( Xi0Lin_intg( q ) ); 
   Xi0loop_buf.push_back( Xi0loop_intg( q ) ); 
   Xi0eft_buf.push_back( Xi0eft_intg( q ) ); 
   
   Xi0LinZeroLag_buf.push_back( Xi0LinZeroLag_intg( q ) );
   Xi0loopZeroLag_buf.push_back( Xi0loopZeroLag_intg( q ) );   
   

   Xi1loop_buf.push_back( Xi1loop_intg( q ) ); 
   Xi1eft_buf.push_back( Xi1eft_intg( q ) ); 
   
   Xi2Lin_buf.push_back( Xi2Lin_intg( q ) ); 
   Xi2loop_buf.push_back( Xi2loop_intg( q ) ); 
   Xi2eft_buf.push_back( Xi2eft_intg( q ) ); 

   Xi3loop_buf.push_back( Xi3loop_intg( q ) ); 
   Xi3eft_buf.push_back( Xi3eft_intg( q ) ); 
   
   

   

    }

    std::cout << "Various f(q) are obtained." << std::endl;
    return;
}

double q_func::xi_L_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = pow( k, 2 ) * exp( - pow ( k, 2) ) * PL * sph_bessel_j( 0, jx );//Anti-aliasing kernel: exp(-k^2)
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}





double q_func::Xi0Lin_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = PL * sph_bessel_j( 0, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}


double q_func::Xi0LinZeroLag_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	kernel_val = PL ;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}







double q_func::Xi0loop_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val = (9./98. * kf.Q_val( 1, k )  + 10.0 / 21 * kf.R_val( 1, k )) * sph_bessel_j( 0, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}



double q_func::Xi0loopZeroLag_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	kernel_val = 9./98. * kf.Q_val( 1, k )  + 10.0 / 21 * kf.R_val( 1, k );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}




double q_func::Xi0eft_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = k * k * PL * sph_bessel_j( 0, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}






double q_func::Xi1loop_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val =  (-3./7.) / k * (kf.Q_val( 1, k ) - 3.0 *  kf.Q_val( 2, k ) + 2. * kf.R_val( 1, k )
	- 6. * kf.R_val( 2, k ) ) * sph_bessel_j( 1, jx )  ;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}


double q_func::Xi1eft_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = -3./7. * k * PL * sph_bessel_j( 1, jx )  ;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}




double q_func::Xi2Lin_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = PL * sph_bessel_j( 2, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}

double q_func::Xi2loop_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val =  (9./98. * kf.Q_val( 1, k )  + 10.0 / 21. * kf.R_val( 1, k )) * sph_bessel_j( 2, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}


double q_func::Xi2eft_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = k * k * PL * sph_bessel_j( 2, jx );
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}






double q_func::Xi3loop_intg( const double & q )
{
    intg.clear(  );
    double k( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	jx = k * q;
	kernel_val =  (-3./7.) / k * (kf.Q_val( 1, k ) + 2.0 *  kf.Q_val( 2, k ) + 2. * kf.R_val( 1, k )
	+ 4. * kf.R_val( 2, k ) ) * sph_bessel_j( 3, jx )  ;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}


double q_func::Xi3eft_intg( const double & q )
{
    intg.clear(  );
    double k( 0. ), PL( 0. );
    double jx( 0. );			// Argument of spherical bessel
    double kernel_val( 0. );
    for( unsigned i = 0; i < k_intg_buf.size(  ); ++ i )
    {
	k = k_intg_buf[ i ];
	PL = kf.PL_val( k );
	jx = k * q;
	kernel_val = -3./7. * k * PL * sph_bessel_j( 3, jx )  ;
	intg.read( k, kernel_val );
    }

    return intg.result(  ) * one_over_pi2;
}








double q_func::interp_val( const double & q, const int & i,
			   const std::vector<double> & vec )
{
    if( i < 0 || i > int( q_buf.size(  ) - 2 ) )
	return 0.;
    else if( ( q_buf[ i+1 ] - q_buf[ i ] ) / q_buf[ i ]
	     < nearly_0 )
	return vec[ i ];
    else
	return vec[ i ] + ( vec[ i+1 ] - vec[ i ] )
	    / ( q_buf[ i+1 ] - q_buf[ i ] )
	    * ( q - q_buf[ i ] );
}







void q_func::var_func( const double & q, q_func_vals & res )
{
	

	
    std::map<double, int>::iterator p
	= q_index_buf.lower_bound( q );
    if( p != q_index_buf.begin(  ) )
	-- p;

    const int i = p->second;

    res.xi_L = interp_val( q, i, xi_L_buf );

   res.Xi0Lin = interp_val( q, i, Xi0Lin_buf ); 
   res.Xi0loop = interp_val( q, i, Xi0loop_buf ); 
   res.Xi0eft = interp_val( q, i, Xi0eft_buf ); 

   res.Xi1loop = interp_val( q, i, Xi1loop_buf ); 
   res.Xi1eft = interp_val( q, i, Xi1eft_buf );
   
   res.Xi2Lin = interp_val( q, i, Xi2Lin_buf ); 
   res.Xi2loop = interp_val( q, i, Xi2loop_buf ); 
   res.Xi2eft = interp_val( q, i, Xi2eft_buf ); 
   
   res.Xi3loop = interp_val( q, i, Xi3loop_buf ); 
   res.Xi3eft = interp_val( q, i, Xi3eft_buf ); 
   
   res.Xi0LinZeroLag = interp_val( q, i, Xi0LinZeroLag_buf ); 
   res.Xi0loopZeroLag = interp_val( q, i, Xi0loopZeroLag_buf );


    return;
}

////////////////////////////////////////////////////////////
// Test output

void q_func::save_q_func( std::string file_name )
{
    std::ofstream fout( file_name.c_str(  ) );

    dvec * p[ 14 ]  
	= { &q_buf, &xi_L_buf,  
	    &Xi0Lin_buf, &Xi0loop_buf, &Xi0eft_buf,
	    &Xi1loop_buf, &Xi1eft_buf,
	    &Xi2Lin_buf, &Xi2loop_buf, &Xi2eft_buf,
	    &Xi3loop_buf, &Xi3eft_buf,
	    &Xi0LinZeroLag_buf, &Xi0loopZeroLag_buf }; 


    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	for( int j = 0; j < 14; ++ j )  
	    fout << p[ j ]->at( i ) << ' ';
	fout << '\n';
    }
    fout << std::endl;

    return;
}

void q_func::load_q_func( std::string file_name )
{
    std::ifstream fin( file_name.c_str(  ) );
    if( !fin )
	throw "Unable to open q function file.";

    std::pair<double, int> q_index_pair;

    dvec * p[ 14 ]   
	= { &q_buf, &xi_L_buf,
	     &Xi0Lin_buf, &Xi0loop_buf, &Xi0eft_buf,
	    &Xi1loop_buf, &Xi1eft_buf,
	    &Xi2Lin_buf, &Xi2loop_buf, &Xi2eft_buf,
	    &Xi3loop_buf, &Xi3eft_buf, 
	    &Xi0LinZeroLag_buf, &Xi0loopZeroLag_buf };

	    
	    
    for( int i = 0; i < 14; ++ i )  
	p[ i ]->clear(  );

    double temp( 0. );
    while( !fin.eof(  ) )
    {

	for( int i = 0; i < 14; ++ i )  
	{
	    fin >> temp;
	    p[ i ]->push_back( temp );
	}
    }
    for( int i = 0; i < 14; ++ i ) 
	p[ i ]->pop_back(  );
	
	
	
    for( unsigned i = 0; i < q_buf.size(  ); ++ i )
    {
	q_index_pair.first = q_buf[ i ];
	q_index_pair.second = i;
	q_index_buf.insert( q_index_pair );

    }

    std::cout << "q functions loaded from: "
	      << file_name << std::endl;
    return;
}
