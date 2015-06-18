
############################################################
                        CLEFT-code 
############################################################



Alejandro Aviles
University of California, Berkeley.
aviles@berkeley.edu, avilescervantes@gmail.com

Zvonimir Vlah
Stanford University.
zvlah@stanford.edu




CLEFT-code is an extension of Lile Wang's integration c++ code on CLPT (https://github.com/wll745881210/CLPT_GSRSD) that includes EFT contributions as explained in arXiv:1506.XXXX. It also calculates the \Xi_ell(q) functions defined in eqs.(4.14) of arXiv:1506.05264.

The input is the linear power spectrum at the desire redshift of the output files. It should be included in the data/ directory  (or otherwise specified in the CLEFT/par.ini file). It should be a raw ascii file with columns:
 
#1 = k 
#2 = P(k)

(We define the power spectrum through <\delta(k) \delta(k')> = (2 \pi)^3 \delta_D(k+k') P(k).)

To obtain good precision the input PS should have values up to k_max ~ 100 h/Mpc. 


Run: Compile with the Makefile and run it as ./cleft par.ini


############################################
############################################
Output files:
############################################
############################################

1) xi.dat give the different terms in the CLEFT correlation function (\xi(r)) and also the linear correlation function. The columns are 

#1  = r
#2  = \xi_linear
#3  = \xi_ZA
#4  = \xi_Aloop
#5  = \xi_Aeft_alpha_0
#6  = \xi_Aeft_alpha_1
#7  = \xi_Wloop
#8  = \xi_Weft_alpha_2
#9  = \xi_Weft_alpha_3

Thus, to plot the full theory one should add:

\xi_LEFT = #3 + #4 + #7 + \alpha_0 #5 + \alpha_1 #6 +\alpha_2 #8 +\alpha_3 #9

2) data/k_func.dat output file gives the Q_n(k) and R_n(k) functions as written in appendix A of arXiv:1506.05264. The columns are

#1 = q
#2 = R_1
#3 = R_2
#4 = Q_1
#5 = Q_2
#6 = Q_3

You can use your own k_func.dat input file by writing its relative route in CLEFT/par.ini

3) data/q_func.dat output file gives the \Xi_ell(q) functions as defined in eqs.(4.14) of 1506.05264. Each \Xi_ell is split in linear (if it has), loop and eft pieces. The last two columns gives the linear and loop 1D displacement field dispersions. The columns are:

#1  =  q
#2  =  linear correlation function 
#3  =  \Xi_0_linear(q)
#4  =  \Xi_0_loop(q)
#5  =  \Xi_0_eft(q)
#6  =  \Xi_1_loop(q)
#7  =  \Xi_1_eft(q)
#8  =  \Xi_2_linear(q)
#9  =  \Xi_2_loop(q)
#10 =  \Xi_2_eft(q)
#11 =  \Xi_3_loop(q)
#11 =  \Xi_3_eft(q)
#13 =  \Xi_0_lin(q=0)  = 3*sigma^2_lin 
#14 =  \Xi_0_loop(q=0) = 3*sigma^2_loop 

You can use your own q_func.dat input file by writing its relative route in CLEFT/par.ini

#################################################
#################################################


You can use this code as it is or modify it as you want. Just cite the papers:

1) Vlah, White & Aviles [arXiv:1506.05264]
2) Wang, Reid & White [MNRAS 437 (2014) 588, arXiv:1306:1804] 

You may also want to cite:

3) Carlson, Reid & White [MNRAS 429 (2013) 1674, arXiv:1209:0780]  
4) Vlah, Seljak & Baldauf [Phys.Rev. D91, 023508 (2015), arXiv:1410.1617]








