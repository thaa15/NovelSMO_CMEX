#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME SIGNS
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumDiscStates(S, 4); 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 3); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 2); 
    ssSetNumSampleTimes(S, 1); 

    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE 
    | SS_OPTION_DISCRETE_VALUED_OUTPUT));
} 

static void mdlInitializeSampleTimes(SimStruct *S){ 
    ssSetSampleTime(S, 0, 0.00001); 
    ssSetOffsetTime(S, 0, 0.0);
} 

#define MDL_INITIALIZE_CONDITIONS 
static void mdlInitializeConditions(SimStruct *S){ 
    real_T *X0 = ssGetRealDiscStates(S); 
    int_T nXStates = ssGetNumDiscStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    int_T i; 

    /* initialize the states to 0.0 */ 
    for (i=0; i < nXStates; i++) { 
        X0[i] = 0.0; 
    } 
} 

static void mdlOutputs(SimStruct *S, int_T tid){ 
    real_T *Y = ssGetOutputPortRealSignal(S,0); 
    real_T *X = ssGetRealDiscStates(S); 
   
    Y[0] = X[2];
    Y[1] = X[3];
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) { 
    real_T *X = ssGetRealDiscStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    real_T dt = 0.00001;
    real_T Sa,Sb,sgn_Sa,sgn_Sb;
    real_T i_hatA_input,i_hatB_input;
    real_T i_hatA_old,i_hatA_new,i_hatB_old,i_hatB_new;

	real_T Ls	 = 0.66;
	real_T Rs	 = 0.251;
    
    real_T k    = U(2);
    i_hatA_input = U(0);
    i_hatB_input = U(1);
    
    i_hatA_old = X[0];
    i_hatB_old = X[1];
    i_hatA_new = i_hatA_old + (dt*(i_hatA_input - i_hatA_old));
    i_hatB_new = i_hatB_old + (dt*(i_hatB_input - i_hatB_old));
    
    Sa = (Rs/Ls)*i_hatA_new + i_hatA_input;
    Sb = (Rs/Ls)*i_hatB_new + i_hatB_input;
    if(Sa > 0){
		sgn_Sa = 1;
	}
	else if(Sa < 0){
		sgn_Sa = -1;
	}
	else{
		sgn_Sa = 1e-6;
	}
    if(Sb > 0){
		sgn_Sb = 1;
	}
	else if(Sb < 0){
		sgn_Sb = -1;
	}
	else{
		sgn_Sb = 1e-6;
	}
    
    X[0] = i_hatA_new;
    X[1] = i_hatB_new;
    X[2] = k*sgn_Sa;
    X[3] = k*sgn_Sb;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 
