#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME LESOPLL
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element])

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumContStates(S, 4); 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 6); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 3); 
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

    Y[0] = X[0];
    Y[1] = X[1];
    Y[2] = X[3];
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
    real_T J     = 3.24e-5;
    
    real_T fluks_cos = U(0);
    real_T fluks_sin = U(1);
    real_T Ia = U(2);
    real_T Ib = U(3);
    real_T Ic = U(4);
    real_T wo = U(5);
    real_T theta_hat,theta_hat_dot;
    real_T we_hat,we_hat_dot;
    real_T f,f_dot;
    real_T Ialfa, Ibeta;
    real_T G	 = 0.8164965809;
	real_T H	 = 0.8660254038;
    
    theta_hat = X[0];
    we_hat = X[1];
    f = X[2];
    Ialfa = G*(Ia*1 - Ib*0.5 - Ic*0.5);
    Ibeta = G*(Ia*0 + Ib*H - Ic*H);
    //PLL
    real_T error;
    error = (sin(theta_hat)*fluks_sin - cos(theta_hat)*fluks_sin)/psi;
    theta_hat_dot = we_hat - 3*wo*error;
    we_hat_dot = -N*psi*sin(theta_hat)*Ialfa/J + N*psi*cos(theta_hat)*Ibeta/J + f - 3*wo*wo*error;
    f_dot = -(wo*wo*wo)*error;
    
    dX[0] = theta_hat_dot;
    dX[1] = we_hat_dot;
    dX[2] = f_dot;
    dX[3] = error;
}

static void mdlTerminate(SimStruct *S) 
{} /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /* MEX-file interface mechanism */ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 