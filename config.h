#ifndef CONFIG_H_
#define CONFIG_H_

// 编译期开关（-D 覆盖）
#define NN_MODEL 1
#ifndef NNN_MODEL
#define NNN_MODEL 0
#endif
#ifndef SOC_MODEL
#define SOC_MODEL 0
#endif
#ifndef SOC_SPIN
#define SOC_SPIN 0
#endif

// 能带路径：0 = M-Γ-K-M'-K'（对 VASP），1 = Γ-K-M-Γ（对 wannier90 的 GW 能带）
#ifndef PATH_MODE
#define PATH_MODE 0
#endif

#include "params.h"     // 参数分支依赖 NN_MODEL / NNN_MODEL

#define PI 3.14159265358979323846

#endif
