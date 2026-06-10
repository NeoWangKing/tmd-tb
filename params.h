// params.h
#ifndef PARAMS_H_
#define PARAMS_H_

#ifndef NN_MODEL
#define NN_MODEL 1
#endif

#ifndef TNN_MODEL
#define TNN_MODEL 1
#endif

#ifndef SOC_MODEL
#define SOC_MODEL 0
#endif

#if NN_MODEL
    #if TNN_MODEL
        static const double e1  =  0.683;
        static const double e2  =  1.707;

        static const double t11 = -0.146;
        static const double t12 =  0.114;
        static const double t13 = -0.506;
        static const double t22 =  0.085;
        static const double t23 =  0.162;
        static const double t33 =  0.073;

        static const double r11 =  0.060;
        static const double r12 =  0.236;
        static const double r21 = -0.067;
        static const double r22 =  0.016;
        static const double r23 =  0.087;

        static const double u11 = -0.038;
        static const double u12 = -0.046;
        static const double u13 = -0.001;
        static const double u22 =  0.266;
        static const double u23 = -0.176;
        static const double u33 = -0.150;
    #else
        // 仅最近邻模型参数
        static const double e1  =  1.046;
        static const double e2  =  2.104;
        static const double t11 = -0.184;
        static const double t12 =  0.401;
        static const double t13 =  0.507;
        static const double t22 =  0.218;
        static const double t23 =  0.338;
        static const double t33 =  0.057;
    #endif
#endif

#if SOC_MODEL
    // MoS2 自旋‑轨道耦合强度 (eV)
    static const double lambda = 0.073;
#endif

#endif // PARAMS_H_
