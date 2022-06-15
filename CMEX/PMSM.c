#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME PMSM
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) 

static void mdlInitializeSizes(SimStruct *S)
{ 
	ssSetNumContStates(S, 6);                 
	if (!ssSetNumInputPorts(S, 1)) return;
	ssSetInputPortWidth(S, 0, 4);               
	ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortOverWritable(S, 0, 1);
	
	if (!ssSetNumOutputPorts(S, 1)) return;
	ssSetOutputPortWidth(S, 0, 8);              
	ssSetNumSampleTimes(S, 1);
    
	ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 	
	| SS_OPTION_DISCRETE_VALUED_OUTPUT));			
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
  	ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
  	ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
static void mdlInitializeConditions(SimStruct *S)
{
  	real_T 	*X0 = ssGetContStates(S);
  	int_T 	nStates = ssGetNumContStates(S);
  	int_T 	i;                                      

  	for (i=0; i < nStates; i++)                     
	{
		X0[i] = 0.0;
	}
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

	real_T Ialfa, Ibeta;
	real_T Wm, Te, Tetae;
	real_T Ia, Ib, Ic, Isq, Isd;
	
	real_T G	 = 0.8164965809;
	real_T H	 = 0.8660254038;
	real_T N	 = 4;
	real_T Fluxr = 16.8e-3;
    
	Wm 			= X[0];
	Tetae 		= X[1];
	Ialfa 		= X[2];
	Ibeta 		= X[3];
    Isd         = X[4];
    Isq         = X[5];
    
    real_T I_Alfa = Isd*cos(Tetae) - Isq*sin(Tetae);
    real_T I_Beta = Isd*sin(Tetae) + Isq*cos(Tetae);
	Ia = G*I_Alfa; 
	Ib = G*(-0.5*I_Alfa + H*I_Beta); 
	Ic = G*(-0.5*I_Alfa - H*I_Beta);
	Te 			= -N*Fluxr*sin(Tetae)*Ialfa + N*Fluxr*cos(Tetae)*Ibeta;

	Y[0]		= Te;
	Y[1]		= Wm;
	Y[2]		= Tetae;
	Y[3]		= Ia;
	Y[4]		= Ib;
	Y[5]		= Ic;
    Y[6]        = I_Alfa;
    Y[7]        = I_Beta;

} 

#define MDL_DERIVATIVES
static void mdlDerivatives(SimStruct *S)
{
  	real_T  *dX = ssGetdX(S); 										
  	real_T  *X = ssGetContStates(S); 							
  	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);	

  	real_T Va,Vb,Vc;
	real_T Te, TL;
	real_T Wm, Wm_dot;

	real_T Valfa, Vbeta, Vsd, Vsq;
	real_T Ialfa, Ibeta;
	real_T Ialfa_dot, Ibeta_dot;
	real_T Tetae, Tetae_dot;

	real_T Jr    = 3.24e-5;
	real_T N	 = 4;
	real_T Lr	 = 0.66;
	real_T Rr	 = 0.251;
	real_T Fluxr = 16.8e-3;
    
	real_T G	 = 0.8164965809;
	real_T H	 = 0.8660254038;
    real_T Isd, Isq, Isd_dot, Isq_dot;
	
  	Va 			= U(0);
  	Vb			= U(1);
	Vc          = U(2);
	TL 			= U(3);

  	Wm          = X[0];
  	Tetae		= X[1];
	Ialfa 		= X[2];
	Ibeta 		= X[3];
    Isd         = X[4];
    Isq         = X[5];

	Valfa		= G*(Va - 0.5*Vb - 0.5*Vc);
	Vbeta		= G*(H*Vb - H*Vc);

	Te 			= -N*N*Fluxr*sin(Tetae)*Ialfa + N*N*Fluxr*cos(Tetae)*Ibeta;
	Wm_dot		= (Te - TL)/Jr;
	Tetae_dot	= N*Wm;
    Vsd         = Valfa*cos(Tetae) + Vbeta*sin(Tetae);
    Vsq         = -Valfa*sin(Tetae) + Vbeta*cos(Tetae);
	Ialfa_dot	= (1/Lr)*(Valfa + N*Fluxr*Wm*sin(Tetae) - Rr*Ialfa);
	Ibeta_dot	= (1/Lr)*(Vbeta - N*Fluxr*Wm*cos(Tetae) - Rr*Ibeta);
    Isd_dot     = (Vsd - Rr*Isd + N*Wm*Lr*Isq)/Lr;
    Isq_dot     = (Vsq - Rr*Isq - N*Wm*(Lr*Isd+Fluxr))/Lr;

  	dX[0]		= Wm_dot;
  	dX[1]		= Tetae_dot;
	dX[2]		= Ialfa_dot;
	dX[3]		= Ibeta_dot;
    dX[4]       = Isd_dot;
    dX[5]       = Isq_dot;
}


static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this fILe being compILed as a MEX-fILe? */ 
#include "simulink.c" /*MEX-fILe interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 

