#ifndef INTEGRAL_H
#define INTEGRAL_H

class integral
{
	////////// Con-/destructor & initializer//////////
public:
	integral(  );
	~integral(  );
	void clear(  );

	////////// Data buffer //////////
private:
	int counter;
	double x_buf[ 2 ];
	double y_buf[ 2 ];
	double intg_res;
	


    static const double qbuf2[ ]; //Modificacion	

	////////// Read data //////////
public:
	void read( const double & x, const double & y );
	double midpoint_qi( const int & i );  //Modificacion	
	static const int q_intnum; // Modificacion	


	////////// Feed back //////////
public:
	double result(  );

	////////// Gauss-Legendre //////////
public:
	void gl_clear(  );
	double gl_xi( const int & i );
	void gl_read( const int & i, const double & kernel );
	double gl_result(  );
	static const int gl_num;

private:
	double gl_intg_res;
	static const double x[  ];
    static const double w[  ];


};

#endif

