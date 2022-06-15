#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME SMO
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element])

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumContStates(S, 4); 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 9); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 2); 
    ssSetNumSampleTimes(S, 1); 

    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE); 
} 

static void mdlInitializeSampleTimes(SimStruct *S) { 
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME); 
    ssSetOffsetTime(S, 0, 0.0); 
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S) { 

    real_T *X0 = ssGetContStates(S); 
    int_T nStates = ssGetNumContStates(S); 
    int_T i; 

    /* initialize the states to 0.0 */ 
    for (i=0; i < nStates; i++) {X0[i] = 0.0;} 
} 

static void mdlOutputs(SimStruct *S, int_T tid) { 
	real_T *Y = ssGetOutputPortRealSignal(S,0); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 

    Y[0] = X[2];
    Y[1] = X[3];
} 

#define MDL_DERIVATIVES 
static void mdlDerivatives(SimStruct *S) { 
	real_T *dX = ssGetdX(S); 
	real_T *X = ssGetContStates(S); 
	InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
	
	// PMSM Dynamics
    real_T N	 = 4;
	real_T Ls	 = 0.66;
	real_T Rs	 = 0.251;
	real_T psi   = 16.8e-3;
    

    real_T Va = U(0);
    real_T Vb = U(1);
    real_T Vc = U(2);
    real_T Ia = U(3);
    real_T Ib = U(4);
    real_T Ic = U(5);
    real_T we = U(6);
    real_T fluks_sin = U(7);
    real_T fluks_cos = U(8);
    
    real_T G 	    = 0.8164966; //sqrt (2/3)
    real_T H 	    = 0.866025; //(sqrt 3)/2
    real_T ua,ub,ia,ib;
    ia = G*(Ia*1 - Ib*0.5 - Ic*0.5);
    ib = G*(Ia*0 + Ib*H - Ic*H);
    ua = G*(Va*1 - Vb*0.5 - Vc*0.5);
    ub = G*(Va*0 + Vb*H - Vc*H);
    
    //SLIDING MODE OBSERVER
    real_T i_hatA,i_hatB,i_hatA_dot,i_hatB_dot;
    i_hatA = X[0];
    i_hatB = X[1];
    
    i_hatA_dot = (1/Ls)*ua + (-Rs/Ls)*i_hatA - (we/Ls)*fluks_sin;
    i_hatB_dot = (1/Ls)*ub + (-Rs/Ls)*i_hatB - (we/Ls)*fluks_cos;
    
    dX[0] = i_hatA_dot;
    dX[1] = i_hatB_dot;
    
    //SLIDING SURFACE
    real_T i_flagA,i_flagB,i_flagA_dot,i_flagB_dot;
    i_flagA = i_hatA - ia;
    i_flagB = i_hatB - ib;
    i_flagA_dot = (-Rs/Ls)*i_flagA + (1/Ls)*we*fluks_sin - (we/Ls)*fluks_sin;
    i_flagB_dot = (-Rs/Ls)*i_flagB + (1/Ls)*we*fluks_cos - (we/Ls)*fluks_cos;
    dX[2] = i_flagA_dot;
    dX[3] = i_flagB_dot;
}

static void mdlTerminate(SimStruct *S) 
{} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 