#include "corr_func.h"
//~ #include "lu_decomp.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

////////////////////////////////////////////////////////////
// Constructor, destructor and initializer

	
corr_func::corr_func(  )
{

}

corr_func::~corr_func(  )
{
	
}
	
void corr_func::set_par( const corr_func_init & c_arg,
						 const q_func & qf )
{
	this->r_max = c_arg.r_max;
	this->r_min = c_arg.r_min;
	this->r_bin_num = c_arg.r_bin_num;
	this->xi_file_name = c_arg.file_name;
	this->qf = ( q_func * )( &qf );

	std::cout << std::endl;
 	return;
}

////////////////////////////////////////////////////////////
// Integration kernel

void corr_func::M( const double & r, const vec3 & q )
{
	double q_norm = sqrt( q.x*q.x + q.y*q.y + q.z*q.z );
	const double q_vec[ 3 ] = { q.x, q.y, q.z };
	const double qh[ 3 ]
		= { q.x / q_norm, q.y / q_norm, q.z / q_norm };
	const double r_vec[ 3 ] = { 0, 0, r };
	
	q_func_vals qfv;
	qf->var_func( q_norm, qfv );
	


	////////////////////////////////////////////////////////////
/////////////////Choose the decomposition/////////////////

	double B_inv[ 3 ][ 3 ];
	double Cloop[ 3 ][ 3 ];
	double Ceft0[ 3 ][ 3 ];
	double Ceftq[ 3 ][ 3 ];
	double Bdetinvh( 0. );
	
	double Xup, Yup, Xloop, Yloop, f_X, h_Y;
	
	double Xlin0lag, Xlinq, Ylin;
	double Xlpt0lag, Xlptq, Ylpt;
	double Xeft0lag, Xeftq, Yeft;
	double Vloop, Tloop;
	double Veftalpha1, Teftalpha1;
	double Veftalpha3, Teftalpha3;
	
	Xlin0lag = 0.66666667 *qfv.Xi0LinZeroLag ;
	Xlinq = -0.66666666*(qfv.Xi0Lin + qfv.Xi2Lin);
	Ylin =   2. * qfv.Xi2Lin;
	
	Xlpt0lag = 0.66666666 * qfv.Xi0loopZeroLag;
	Xlptq = -0.66666666*(qfv.Xi0loop + qfv.Xi2loop);
	Ylpt =   2. * qfv.Xi2loop;	
	
	Xeft0lag = 2.0 / 3.0;
	Xeftq = -0.66666666*(qfv.Xi0eft + qfv.Xi2eft);
	Yeft  = 2.* qfv.Xi2eft;
	
	Vloop = 2./5. * qfv.Xi1loop - 3./5. * qfv.Xi3loop;
	Tloop = 3. * qfv.Xi3loop;
	
	Veftalpha1 = 2./5. * qfv.Xi1eft;
	Teftalpha1 = 0.0;
	
	Veftalpha3 = - 3./5. * qfv.Xi3eft;
	Teftalpha3 = 3. * qfv.Xieft;	
	
	
	Xup = Xlin0lag + Xlinq;
	Yup = Ylin;
	
	Xloop = Xlpt0lag + Xlptq;
	Yloop = Ylpt;
	
	
	
	

		
	
	f_X = 1/Xup;
	h_Y = -Yup / ( Xup*Xup + Xup*Yup );
	
	Bdetinvh = f_X*sqrt( f_X + h_Y );
	
	for( int i = 0; i < 3; ++ i ){
		for( int j = 0; j < 3; ++ j ){
	 	B_inv[ i ][ j ] = f_X * delta_k( i, j ) + h_Y * qh[ i ] * qh[ j ];
	 	Cloop[ i ][ j ] = Xloop * delta_k( i, j ) + Yloop * qh[ i ] * qh[ j ]; 
	 	Ceft0[ i ][ j ] =  Xeft0lag * delta_k( i, j );
	 	Ceftq[ i ][ j ] = Xeftq * delta_k( i, j ) + Yeft * qh[ i ] * qh[ j ]; 
		};
	};          
			
	
		
	
	
	
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////	
	
	



	////////// Vectors and tensors //////////



	double W[ 3 ][ 3 ][ 3 ];
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )
				W[ i ][ j ][ k ]
					= Vloop * ( qh[ i ] * delta_k( j, k )
					+ qh[ j ] * delta_k( k, i )+  qh[ k ] * delta_k( i, j ) )
					+ Tloop  * qh[ i ] * qh[ j ] * qh[ k ];
				 	

    double Weftalpha1[ 3 ][ 3 ][ 3 ];
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )			    
    Weftalpha1[ i ][ j ][ k ] =  Veftalpha1 * ( qh[ i ] * delta_k( j, k )
					+ qh[ j ] * delta_k( k, i )+  qh[ k ] * delta_k( i, j ) )
					+ Teftalpha1  * qh[ i ] * qh[ j ] * qh[ k ];
    
    double Weftalpha3[ 3 ][ 3 ][ 3 ];
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )			    
    Weftalpha3[ i ][ j ][ k ] =  Veftalpha3 * ( qh[ i ] * delta_k( j, k )
					+ qh[ j ] * delta_k( k, i )+  qh[ k ] * delta_k( i, j ) )
					+ Teftalpha3  * qh[ i ] * qh[ j ] * qh[ k ];   
    
	
	double g[ 3 ];
	for( int i = 0; i < 3; ++ i )
	{
		g[ i ] = 0.;
		for( int j = 0; j < 3; ++ j )
			g[ i ] += B_inv[ i ][ j ]
				* ( q_vec[ j ] - r_vec[ j ] );
	}

	double G[ 3 ][ 3 ];
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			G[ i ][ j ] = B_inv[ i ][ j ]
				- g[ i ] * g[ j ];

	double Gamma[ 3 ][ 3 ][ 3 ];
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )
				Gamma[ i ][ j ][ k ]
					= B_inv[ j ][ k ] * g[ i ]
					+ B_inv[ k ][ i ] * g[ j ]
					+ B_inv[ i ][ j ] * g[ k ]
					- g[ i ] * g[ j ] * g[ k ];
					
					
	    
	    
				

	////////// Sum them up! //////////

	for( int i = 0; i < num_bias_comp; ++ i )
		bias_comp_inner[ i ] = 0.;
	
	double temp( 0. );
	bias_comp_inner[ 0 ] += 1.;
	
	for( int i = 0; i < 3; ++ i ){
		for( int j = 0; j < 3; ++ j ){
		bias_comp_inner[ 1 ] += - Cloop[ i ][ j ] * G[ i ][ j ] / 2.;
		bias_comp_inner[ 2 ] += - Ceft0[ i ][ j ] * G[ i ][ j ] / 2.;
		bias_comp_inner[ 3 ] += - Ceftq[ i ][ j ] * G[ i ][ j ] / 2.;
		};
	};


	temp = 0.;
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )
				temp += W[ i ][ j ][ k ]
					* Gamma[ i ][ j ][ k ];
		bias_comp_inner[ 4 ] +=  temp / 6.;
		
		
	temp = 0.;
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )
				temp += Weftalpha1[ i ][ j ][ k ]
					* Gamma[ i ][ j ][ k ];
		bias_comp_inner[ 5 ] +=  temp / 6.;		
		
	temp = 0.;
	for( int i = 0; i < 3; ++ i )
		for( int j = 0; j < 3; ++ j )
			for( int k = 0; k < 3; ++ k )
				temp += Weftalpha3[ i ][ j ][ k ]
					* Gamma[ i ][ j ][ k ];
		bias_comp_inner[ 6 ] +=  temp / 6.;			


	temp = 0.;
	for( int i = 0; i < 3; ++ i )
		temp += ( q_vec[ i ] - r_vec[ i ] ) * g[ i ];

	static const double two_pi_cube = 248.05021344239853;
	// ( 2 \pi )^3
	const double gauss
		= exp( -0.5 * temp ) * Bdetinvh
		/ sqrt( two_pi_cube );
	for( int i = 0; i < num_bias_comp; ++ i )
		bias_comp_inner[ i ] *= gauss;
	
	return;
}


void corr_func::M( const double & r, const double & q,
				   const double & mu )
{
	vec3 qv;
	qv.x = q * sqrt( 1 - mu*mu );
	qv.y = 0;
	qv.z = q * mu;
	M( r, qv );
	return;
}

int corr_func::delta_k( const int & i, const int & j )
{
	return ( i == j ? 1 : 0 );
}

////////////////////////////////////////////////////////////
// Correlation function xi
//~ 


void corr_func::xi( const double & r )
{


	double q( 0. );
	double y( 0. );
	double mu( 0. ), beta( 0. );

	integral intg[ num_bias_comp ];
	for( int i = 0; i < num_bias_comp; ++ i )
		intg[ i ].clear(  );
	

	for( int i = 0; i < intg[ 0 ].q_intnum; ++ i )   
	{

		y = intg[ 0 ].midpoint_qi( i );   
		
		for( int k = 0; k < num_bias_comp; ++ k )
			intg[ k ].gl_clear(  );
		for( int j = 0; j < intg[ 0 ].gl_num; ++ j )
		{
			beta = intg[ 0 ].gl_xi( j );
			q = sqrt( y*y + r*r + 2 * r * y * beta );
			if( q < min_q_for_integration
				|| q > max_q_for_integration )
				continue;

			
			mu = ( r + y * beta ) / q;
			M( r, q, mu );
			
			
			for( int k = 0; k < num_bias_comp; ++ k )
			{
				q = bias_comp_inner[ k ];
				intg[ k ].gl_read( j, q );
			}
		}
		for( int k = 0; k < num_bias_comp; ++ k )
		{
			q = pow( y, 2 ) * intg[ k ].gl_result(  );
			intg[ k ].read( y, q );
		}
	}

	for( int k = 0; k < num_bias_comp; ++ k )
	{
		bias_comp_outer[ k ]
			= 2 * pi * intg[ k ].result(  );
		if( k == 0 )
			bias_comp_outer[ k ] -= 1.;
	}
	return;
}

double corr_func::xi_L( const double & r )
{
	static q_func_vals qfv;
	qf->var_func( r, qfv );
	
	return qfv.xi_L;
}

void corr_func::get_xi(  )
{
	double r( 0. );
	
	std::cout << "Generate xi: ";
	

	
	std::cout.flush(  );

	std::vector<double> * p[ num_bias_comp ]
		={& b10b20, & b11b20, & b10b21,
		  & b12b20, & b11b21, & b10b22,
		  & newb};

	pg.init( r_bin_num );
	for( int i = 0; i < r_bin_num; ++ i )
	{
		pg.show( i );
		r = ( r_max - r_min ) / ( r_bin_num - 1. )
			* float( i ) + r_min;
		xi( r );
		r_buf.push_back( r );
		for( int j = 0; j < num_bias_comp; ++ j )
			p[ j ]->push_back( bias_comp_outer[ j ] );
		xi_L_buf.push_back( xi_L( r ) );
	}



	std::cout << "Correlation function obtained. "
			  << std::endl;
	return;
}

////////////////////////////////////////////////////////////
// Output

void corr_func::output(  )
{
	std::ofstream fout( xi_file_name.c_str(  ) );
	std::vector<double> * p[ num_bias_comp + 2 ]
		={& r_buf, &xi_L_buf,
		  & b10b20, & b11b20, & b10b21,
		  & b12b20, & b11b21, & b10b22,
		  & newb};
	
	for( unsigned i = 0; i < r_buf.size(  ); ++ i )
	{
		for( int j = 0; j < num_bias_comp + 2 ; ++ j )
			fout << p[ j ]->at( i ) << ' ';
		fout << '\n';
	}
	fout.flush(  );

	return;
}


