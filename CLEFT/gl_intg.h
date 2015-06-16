#ifndef GL_INTG_H
#define GL_INTG_H

////////////////////////////////////////////////////////////
// Gauss-Legendre quadrature

// const int integral::gl_num = 32;
// const double integral::x[ 32 ]
// = { -0.99726386, -0.98561151, -0.96476226, 
// 	-0.93490608, -0.89632116, -0.84936761, 
// 	-0.7944838 , -0.73218212, -0.66304427, 
// 	-0.58771576, -0.50689991, -0.42135128, 
// 	-0.3318686 , -0.23928736, -0.14447196, 
// 	-0.04830767,  0.04830767,  0.14447196,  
// 	0.23928736,  0.3318686 , 0.42135128,  
// 	0.50689991,  0.58771576,  0.66304427,  
// 	0.73218212, 0.7944838 ,  0.84936761,  
// 	0.89632116,  0.93490608,  0.96476226, 
// 	0.98561151,  0.99726386 };
// const double integral::w[ 32 ]
// = { 0.00701861,  0.01627439,  0.02539207,  
// 	0.03427386,  0.0428359 , 0.05099806,  
// 	0.05868409,  0.06582222,  0.07234579, 
// 	0.0781939 , 0.08331192,  0.08765209,  
// 	0.09117388,  0.0938444 ,  0.09563872, 
// 	0.09654009,  0.09654009,  0.09563872, 
// 	0.0938444 ,  0.09117388, 0.08765209,  
// 	0.08331192,  0.0781939 ,  0.07234579,  
// 	0.06582222, 0.05868409,  0.05099806,  
// 	0.0428359 ,  0.03427386,  0.02539207, 
// 	0.01627439,  0.00701861 };


const int integral::gl_num = 128;
const double integral::x[ 128 ] =
{-0.9998248879471319,-0.9990774599773759,-0.997733248625514,
 -0.9957927585349812,-0.9932571129002129,-0.9901278184917344,
 -0.9864067427245862,-0.9820961084357185,-0.9771984914639074,
 -0.9717168187471366,-0.9656543664319653,-0.9590147578536999,
 -0.9518019613412644,-0.9440202878302202,-0.9356743882779164,
 -0.9267692508789478,-0.9173101980809605,-0.9073028834017568,
 -0.8967532880491582,-0.8856677173453972,-0.8740527969580318,
 -0.8619154689395485,-0.849262987577969,-0.8361029150609068,
 -0.8224431169556438,-0.8082917575079137,-0.7936572947621933,
 -0.778548475506412,-0.7629743300440947,-0.746944166797062,
 -0.7304675667419088,-0.7135543776835874,-0.6962147083695143,
 -0.6784589224477193,-0.6602976322726461,-0.6417416925623076,
 -0.6228021939105849,-0.6034904561585486,-0.5838180216287631,
 -0.5637966482266181,-0.5434383024128104,-0.5227551520511755,
 -0.5017595591361445,-0.480464072404172,-0.4588814198335522,
 -0.4370245010371042,-0.414906379552275,-0.3925402750332674,
 -0.369939555349859,-0.3471177285976355,-0.3240884350244134,
 -0.3008654388776772,-0.2774626201779044,-0.2538939664226943,
 -0.23017356422666,-0.2063155909020792,-0.1823343059853372,
 -0.1582440427142249,-0.1340591994611878,-0.1097942311276437,
 -0.0854636405045155,-0.06108196960413957,-0.03666379096873349,
 -0.01222369896061576,0.01222369896061576,0.03666379096873349,
 0.06108196960413957,0.0854636405045155,0.1097942311276437,
 0.1340591994611878,0.1582440427142249,0.1823343059853372,
 0.2063155909020792,0.23017356422666,0.2538939664226943,
 0.2774626201779044,0.3008654388776772,0.3240884350244134,
 0.3471177285976355,0.369939555349859,0.3925402750332674,
 0.414906379552275,0.4370245010371042,0.4588814198335522,
 0.480464072404172,0.5017595591361445,0.5227551520511755,
 0.5434383024128104,0.5637966482266181,0.5838180216287631,
 0.6034904561585486,0.6228021939105849,0.6417416925623076,
 0.6602976322726461,0.6784589224477193,0.6962147083695143,
 0.7135543776835874,0.7304675667419088,0.746944166797062,
 0.7629743300440947,0.778548475506412,0.7936572947621933,
 0.8082917575079137,0.8224431169556438,0.8361029150609068,
 0.849262987577969,0.8619154689395485,0.8740527969580318,
 0.8856677173453972,0.8967532880491582,0.9073028834017568,
 0.9173101980809605,0.9267692508789478,0.9356743882779164,
 0.9440202878302202,0.9518019613412644,0.9590147578536999,
 0.9656543664319653,0.9717168187471366,0.9771984914639074,
 0.9820961084357185,0.9864067427245862,0.9901278184917344,
 0.9932571129002129,0.9957927585349812,0.997733248625514,
 0.9990774599773759,0.9998248879471319};
const double integral::w[ 128 ]=
{0.0004493809602920904,0.001045812679340349,
 0.00164250301866903,0.002238288430962619,0.002832751471457991,
 0.003425526040910216,0.004016254983738642,0.004604584256702955,
 0.00519016183267633,0.005772637542865699,0.006351663161707189,
 0.006926892566898814,0.007497981925634729,0.008064589890486058,
 0.00862637779861675,0.009183009871660874,0.009734153415006806,
 0.01027947901583216,0.01081866073950308,0.01135137632408042,
 0.01187730737274028,0.01239613954395092,0.01290756273926735,
 0.01341127128861633,0.01390696413295199,0.01439434500416685,
 0.01487312260214731,0.01534301076886514,0.01580372865939935,
 0.01625500090978519,0.0166965578015892,0.01712813542311138,
 0.0175494758271177,0.01796032718500869,0.01836044393733134,
 0.01874958694054471,0.01912752360995095,0.0194940280587066,
 0.01984888123283086,0.02019187104213004,0.02052279248696007,
 0.02084144778075115,0.02114764646822135,0.02144120553920846,
 0.02172194953805208,0.02198971066846049,0.02224432889379977,
 0.02248565203274497,0.02271353585023646,0.02292784414368685,
 0.02312844882438703,0.02331522999406276,0.02348807601653591,
 0.02364688358444762,0.0237915577810034,0.02392201213670346,
 0.02403816868102405,0.02413995798901928,0.02422731922281525,
 0.02430020016797187,0.02435855726469063,0.02440235563384958,
 0.02443156909785005,0.02444618019626252,0.02444618019626252,
 0.02443156909785005,0.02440235563384958,0.02435855726469063,
 0.02430020016797187,0.02422731922281525,0.02413995798901928,
 0.02403816868102405,0.02392201213670346,0.0237915577810034,
 0.02364688358444762,0.02348807601653591,0.02331522999406276,
 0.02312844882438703,0.02292784414368685,0.02271353585023646,
 0.02248565203274497,0.02224432889379977,0.02198971066846049,
 0.02172194953805208,0.02144120553920846,0.02114764646822135,
 0.02084144778075115,0.02052279248696007,0.02019187104213004,
 0.01984888123283086,0.0194940280587066,0.01912752360995095,
 0.01874958694054471,0.01836044393733134,0.01796032718500869,
 0.0175494758271177,0.01712813542311138,0.0166965578015892,
 0.01625500090978519,0.01580372865939935,0.01534301076886514,
 0.01487312260214731,0.01439434500416685,0.01390696413295199,
 0.01341127128861633,0.01290756273926735,0.01239613954395092,
 0.01187730737274028,0.01135137632408042,0.01081866073950308,
 0.01027947901583216,0.009734153415006806,0.009183009871660874,
 0.00862637779861675,0.008064589890486058,0.007497981925634729,
 0.006926892566898814,0.006351663161707189,0.005772637542865699,
 0.00519016183267633,0.004604584256702955,0.004016254983738642,
 0.003425526040910216,0.002832751471457991,0.002238288430962619,
 0.00164250301866903,0.001045812679340349,0.0004493809602920904};




//Modificacion //To make q integration in corr_func::xi  function
const int integral::q_intnum = 1000;
const double integral::qbuf2[ 1000 ] =
 {0.01,0.41039,0.810781,1.21117,1.61156,2.01195,2.41234,2.81273,3.21312,3.61351,4.0139,4.41429,4.81468,5.21508,5.61547,6.01586,
	 6.41625,6.81664,7.21703,7.61742,8.01781,8.4182,8.81859,9.21898,9.61937,10.0198,10.4202,10.8205,11.2209,11.6213,12.0217,
	 12.4221,12.8225,13.2229,13.6233,14.0237,14.4241,14.8244,15.2248,15.6252,16.0256,16.426,16.8264,17.2268,17.6272,18.0276,
	 18.428,18.8283,19.2287,19.6291,20.0295,20.4299,20.8303,21.2307,21.6311,22.0315,22.4319,22.8323,23.2326,23.633,24.0334,
	 24.4338,24.8342,25.2346,25.635,26.0354,26.4358,26.8362,27.2365,27.6369,28.0373,28.4377,28.8381,29.2385,29.6389,30.0393,
	 30.4397,30.8401,31.2405,31.6408,32.0412,32.4416,32.842,33.2424,33.6428,34.0432,34.4436,34.844,35.2444,35.6447,36.0451,
	 36.4455,36.8459,37.2463,37.6467,38.0471,38.4475,38.8479,39.2483,39.6486,40.049,40.4494,40.8498,41.2502,41.6506,42.051,
	 42.4514,42.8518,43.2522,43.6526,44.0529,44.4533,44.8537,45.2541,45.6545,46.0549,46.4553,46.8557,47.2561,47.6565,48.0568,
	 48.4572,48.8576,49.258,49.6584,50.0588,50.4592,50.8596,51.26,51.6604,52.0608,52.4611,52.8615,53.2619,53.6623,54.0627,
	 54.4631,54.8635,55.2639,55.6643,56.0647,56.465,56.8654,57.2658,57.6662,58.0666,58.467,58.8674,59.2678,59.6682,60.0686,
	 60.4689,60.8693,61.2697,61.6701,62.0705,62.4709,62.8713,63.2717,63.6721,64.0725,64.4729,64.8732,65.2736,65.674,66.0744,
	 66.4748,66.8752,67.2756,67.676,68.0764,68.4768,68.8771,69.2775,69.6779,70.0783,70.4787,70.8791,71.2795,71.6799,72.0803,
	 72.4807,72.8811,73.2814,73.6818,74.0822,74.4826,74.883,75.2834,75.6838,76.0842,76.4846,76.885,77.2853,77.6857,78.0861,
	 78.4865,78.8869,79.2873,79.6877,80.0881,80.4885,80.8889,81.2892,81.6896,82.09,82.4904,82.8908,83.2912,83.6916,84.092,
	 84.4924,84.8928,85.2932,85.6935,86.0939,86.4943,86.8947,87.2951,87.6955,88.0959,88.4963,88.8967,89.2971,89.6974,90.0978,90.4982,90.8986,
	 91.299,91.6994,92.0998,92.5002,92.9006,93.301,93.7014,94.1017,94.5021,94.9025,95.3029,95.7033,96.1037,96.5041,96.9045,97.3049,97.7053,
	 98.1056,98.506,98.9064,99.3068,99.7072,100.108,100.508,100.908,101.309,101.709,102.11,102.51,102.91,103.311,103.711,104.112,104.512,
	 104.912,105.313,105.713,106.113,106.514,106.914,107.315,107.715,108.115,108.516,108.916,109.317,109.717,110.117,110.518,110.918,
	 111.319,111.719,112.119,112.52,112.92,113.32,113.721,114.121,114.522,114.922,115.322,115.723,116.123,116.524,116.924,117.324,
	 117.725,118.125,118.526,118.926,119.326,119.727,120.127,120.528,120.928,121.328,121.729,122.129,122.529,122.93,123.33,
	 123.731,124.131,124.531,124.932,125.332,125.733,126.133,126.533,126.934,127.334,127.735,128.135,128.535,128.936,129.336,
	 129.736,130.137,130.537,130.938,131.338,131.738,132.139,132.539,132.94,133.34,133.74,134.141,134.541,134.942,135.342,
	 135.742,136.143,136.543,136.944,137.344,137.744,138.145,138.545,138.945,139.346,139.746,140.147,140.547,140.947,141.348,
	 141.748,142.149,142.549,142.949,143.35,143.75,144.151,144.551,144.951,145.352,145.752,146.152,146.553,146.953,147.354,
	 147.754,148.154,148.555,148.955,149.356,149.756,150.156,150.557,150.957,151.358,151.758,152.158,152.559,152.959,153.36,
	 153.76,154.16,154.561,154.961,155.361,155.762,156.162,156.563,156.963,157.363,157.764,158.164,158.565,158.965,159.365,
	 159.766,160.166,160.567,160.967,161.367,161.768,162.168,162.568,162.969,163.369,163.77,164.17,164.57,164.971,165.371,165.772,
	 166.172,166.572,166.973,167.373,167.774,168.174,168.574,168.975,169.375,169.776,170.176,170.576,170.977,171.377,171.777,172.178,
	 172.578,172.979,173.379,173.779,174.18,174.58,174.981,175.381,175.781,176.182,176.582,176.983,177.383,177.783,178.184,178.584,
	 178.985,179.385,179.785,180.186,180.586,180.986,181.387,181.787,182.188,182.588,182.988,183.389,183.789,184.19,184.59,184.99,
	 185.391,185.791,186.192,186.592,186.992,187.393,187.793,188.193,188.594,188.994,189.395,189.795,190.195,190.596,190.996,191.397,
	 191.797,192.197,192.598,192.998,193.399,193.799,194.199,194.6,195.,195.401,195.801,196.201,196.602,197.002,197.402,197.803,198.203,
	 198.604,199.004,199.404,199.805,200.205,200.606,201.006,201.406,201.807,202.207,202.608,203.008,203.408,203.809,204.209,204.609,205.01,
	 205.41,205.811,206.211,206.611,207.012,207.412,207.813,208.213,208.613,209.014,209.414,209.815,210.215,210.615,211.016,211.416,211.817,
	 212.217,212.617,213.018,213.418,213.818,214.219,214.619,215.02,215.42,215.82,216.221,216.621,217.022,217.422,217.822,218.223,218.623,
	 219.024,219.424,219.824,220.225,220.625,221.025,221.426,221.826,222.227,222.627,223.027,223.428,223.828,224.229,224.629,225.029,
	 225.43,225.83,226.231,226.631,227.031,227.432,227.832,228.233,228.633,229.033,229.434,229.834,230.234,230.635,231.035,231.436,231.836,
	 232.236,232.637,233.037,233.438,233.838,234.238,234.639,235.039,235.44,235.84,236.24,236.641,237.041,237.442,237.842,238.242,238.643,239.043,
239.443,239.844,240.244,240.645,241.045,241.445,241.846,242.246,242.647,243.047,243.447,243.848,244.248,244.649,245.049,245.449,245.85,246.25,
	 246.65,247.051,247.451,247.852,248.252,248.652,249.053,249.453,249.854,250.254,250.654,251.055,251.455,251.856,252.256,252.656,253.057,
	 253.457,253.858,254.258,254.658,255.059,255.459,255.859,256.26,256.66,257.061,257.461,257.861,258.262,258.662,259.063,259.463,259.863,
	 260.264,260.664,261.065,261.465,261.865,262.266,262.666,263.066,263.467,263.867,264.268,264.668,265.068,265.469,265.869,266.27,266.67,
	 267.07,267.471,267.871,268.272,268.672,269.072,269.473,269.873,270.274,270.674,271.074,271.475,271.875,272.275,272.676,273.076,
	 273.477,273.877,274.277,274.678,275.078,275.479,275.879,276.279,276.68,277.08,277.481,277.881,278.281,278.682,279.082,279.482,279.883,
	 280.283,280.684,281.084,281.484,281.885,282.285,282.686,283.086,283.486,283.887,284.287,284.688,285.088,285.488,285.889,286.289,
	 286.69,287.09,287.49,287.891,288.291,288.691,289.092,289.492,289.893,290.293,290.693,291.094,291.494,291.895,292.295,292.695,
	 293.096,293.496,293.897,294.297,294.697,295.098,295.498,295.898,296.299,296.699,297.1,297.5,297.9,298.301,298.701,299.102,299.502,
	 299.902,300.303,300.703,301.104,301.504,301.904,302.305,302.705,303.106,303.506,303.906,304.307,304.707,305.107,305.508,305.908,
	 306.309,306.709,307.109,307.51,307.91,308.311,308.711,309.111,309.512,309.912,310.313,310.713,311.113,311.514,311.914,312.315,
	 312.715,313.115,313.516,313.916,314.316,314.717,315.117,315.518,315.918,316.318,316.719,317.119,317.52,317.92,318.32,318.721,319.121,
	 319.522,319.922,320.322,320.723,321.123,321.523,321.924,322.324,322.725,323.125,323.525,323.926,324.326,324.727,325.127,325.527,325.928,
	 326.328,326.729,327.129,327.529,327.93,328.33,328.731,329.131,329.531,329.932,330.332,330.732,331.133,331.533,331.934,332.334,332.734,
	 333.135,333.535,333.936,334.336,334.736,335.137,335.537,335.938,336.338,336.738,337.139,337.539,337.939,338.34,338.74,339.141,
	 339.541,339.941,340.342,340.742,341.143,341.543,341.943,342.344,342.744,343.145,343.545,343.945,344.346,344.746,345.147,345.547,
	 345.947,346.348,346.748,347.148,347.549,347.949,348.35,348.75,349.15,349.551,349.951,350.352,350.752,351.152,351.553,351.953,352.354,
	 352.754,353.154,353.555,353.955,354.355,354.756,355.156,355.557,355.957,356.357,356.758,357.158,357.559,357.959,358.359,358.76,
	 359.16,359.561,359.961,360.361,360.762,361.162,361.563,361.963,362.363,362.764,363.164,363.564,363.965,364.365,364.766,365.166,
	 365.566,365.967,366.367,366.768,367.168,367.568,367.969,368.369,368.77,369.17,369.57,369.971,370.371,370.772,371.172,371.572,371.973,
	 372.373,372.773,373.174,373.574,373.975,374.375,374.775,375.176,375.576,375.977,376.377,376.777,377.178,377.578,377.979,378.379,378.779,
	 379.18,379.58,379.98,380.381,380.781,381.182,381.582,381.982,382.383,382.783,383.184,383.584,383.984,384.385,384.785,385.186,385.586,
	 385.986,386.387,386.787,387.188,387.588,387.988,388.389,388.789,389.189,389.59,389.99,390.391,390.791,391.191,391.592,391.992,392.393,
	 392.793,393.193,393.594,393.994,394.395,394.795,395.195,395.596,395.996,396.396,
	 396.797,397.197,397.598,397.998,398.398,398.799,399.199,399.6,400.};



#endif

