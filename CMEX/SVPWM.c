#define S_FUNCTION_LEVEL    2
#define S_FUNCTION_NAME     SVPWM
#include "simstruc.h"
#include <math.h>
#define U(element)    (*uPtrs[element]) /*Pointer to Input Port0*/

float DtR(float degrees){
    float radians = degrees * 3.14 / 180;
    return radians;
}

static void mdlInitializeSizes(SimStruct *S){
    if (!ssSetNumInputPorts(S, 1))
        return;
    ssSetInputPortWidth(S, 0, 3);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortOverWritable(S, 0, 1);
    if (!ssSetNumOutputPorts(S, 1))
        return;
    ssSetOutputPortWidth(S, 0, 3);
    ssSetNumSampleTimes(S, 1);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S){
    ssSetSampleTime(S, 0, 0.0001);
    ssSetOffsetTime(S, 0, 0.0);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *Y    = ssGetOutputPortRealSignal(S, 0);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
    int_T             i;
    real_T            t = ssGetT(S);
    real_T            halfVDC = (0.5 * U(0));
    real_T            Fc, Vs, vu, vv, vw, Valfa, Vbeta, Va, Vb, Vc, vab, vbc, vca, van, vbn, vcn;
    real_T            Xt, Yt, Zt;
    real_T            T, T1, Tm, Ta, Tb, Tc, Tcm1, Tcm2, Tcm3;
    real_T            A, B, C, sector;

    int_T             N;
    int_T             Tn;
    real_T            triangle;

    real_T K1 = 0.816497, K2 = 0.866025, K3 = 1.41421356;
    Fc = 1000;
	Vs = U(0);
	Valfa 	= U(1);
	Vbeta 	= U(2);
    T  	= 1 / Fc;
    Va 	= 0; Vb = 0; Vc = 0;

    if (Vbeta > 0)
        A = 1;
    else
        A = 0;
    
    if ((3 * Valfa - Vbeta) > 0)
        B = 1;
    else
        B = 0;
    
    if ((K2 * 2 * Valfa + Vbeta) < 0)
        C = 1;
    else
        C = 0;
    
    N = A + 2 * B + 4 * C;

    switch (N){
	    case 1: sector = 2; break;
	    case 2: sector = 6; break;
	    case 3: sector = 1; break;
	    case 4: sector = 4; break;
	    case 5: sector = 3; break;
	    case 6: sector = 5; break;
    }

    Zt = (T) *(-2 * K2 * Valfa + Vbeta) / (K3 * Vs);
    Yt = (T) * (2 * K2 * Valfa + Vbeta) / (K3 * Vs);
    Xt = (2 * T) * (Vbeta) / (K3 * Vs);
    
    switch (N){
	    case 1: T1 = Zt; Tm = Yt; break;
	    case 2: T1 = Yt; Tm = -Xt; break;
	    case 3: T1 = -Zt; Tm = Xt; break;
	    case 4: T1 = -Xt; Tm = Zt; break;
	    case 5: T1 = Xt; Tm = -Yt; break;
	    case 6: T1 = -Yt; Tm = -Zt; break;
    }
    
    if (T1 + Tm > 1 / Fc){
        T1 = T1 * (T / (T1 + Tm));
        Tm = Tm * (T / (T1 + Tm));
    }

    Ta = (T - T1 - Tm) / 4;
    Tb = Ta + T1 / 2;
    Tc = Tb + Tm / 2;
    switch (N){
	    case 1: Tcm1 = Tb; Tcm2 = Ta; Tcm3 = Tc; break;
	    case 2: Tcm1 = Ta; Tcm2 = Tc; Tcm3 = Tb; break;
	    case 3: Tcm1 = Ta; Tcm2 = Tb; Tcm3 = Tc; break;
	    case 4: Tcm1 = Tc; Tcm2 = Tb; Tcm3 = Ta; break;
	    case 5: Tcm1 = Tc; Tcm2 = Ta; Tcm3 = Tb; break;
	    case 6: Tcm1 = Tb; Tcm2 = Tc; Tcm3 = Ta; break;
	    default: break;
    }

    triangle = T / 2 - fabs(fmod(t, T) - T / 2);

    if (triangle > Tcm1)
        Va = Vs;
    else
        Va = -Vs;
    if (triangle > Tcm2)
        Vb = Vs;
    else
        Vb = -Vs;
    if (triangle > Tcm3)
        Vc = Vs;
    else
        Vc = -Vs;

    vab = Va - Vb;
	vbc = Vb - Vc;
	vca = Vc - Va;
    van = (vab - vca) / 3;
    vbn = (vbc - vab) / 3;
    vcn = (vca - vbc) / 3;

    Y[0] = van;
    Y[1] = vbn;
    Y[2] = vcn;
}

static void mdlTerminate(SimStruct *S){
}   /*Keep this function empty since no memory is allocated*/
#ifdef MATLAB_MEX_FILE
/* Is this file being compiled as a MEX-file? */
#include "simulink.c"
/* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /*Code generation registration function*/
#endif
