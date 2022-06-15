#define S_FUNCTION_LEVEL 2 
#define S_FUNCTION_NAME SPEEDCONTROLLER
#include "simstruc.h" 
#include <math.h> 

#define U(element) (*uPtrs[element]) /*Pointer to Input Port0*/ 

static void mdlInitializeSizes(SimStruct *S){ 
    ssSetNumDiscStates(S, 3); 
    if (!ssSetNumInputPorts(S, 1)) return; 
    ssSetInputPortWidth(S, 0, 4); 
    ssSetInputPortDirectFeedThrough(S, 0, 1); 
    ssSetInputPortOverWritable(S, 0, 1); 
    if (!ssSetNumOutputPorts(S, 1)) return; 
    ssSetOutputPortWidth(S, 0, 1); 
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
    real_T ref;
    
    ref = X[2];
    Y[0] = ref;
} 

#define MDL_UPDATE 
static void mdlUpdate(SimStruct *S, int_T tid) { 
    real_T *X = ssGetRealDiscStates(S); 
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0); 
    real_T dt = 0.00001;
    real_T Kp, Ki,W_ref, W_act, error,theta_m,theta_old;
    real_T integral, integral_old, ip;
    W_ref = U(0);
    Ki = U(1);
    Kp = U(2);
    theta_m = U(3);
    
//     PI Controller
    integral_old = X[0];
    theta_old = X[1];
    W_act = (theta_m - theta_old)/dt;
    error = W_ref - W_act;
    integral = integral_old + error*dt;
    ip = Ki*integral + Kp*W_act;
    
    X[0] = integral;
    X[1] = theta_m;
    X[2] = ip;
}

static void mdlTerminate(SimStruct *S) 
{ } /*Keep this function empty since no memory is allocated*/ 

#ifdef MATLAB_MEX_FILE 
/* Is this file being compiled as a MEX-file? */ 
#include "simulink.c" /*MEX-file interface mechanism*/ 
#else 
#include "cg_sfun.h" /*Code generation registration function*/ 
#endif 
