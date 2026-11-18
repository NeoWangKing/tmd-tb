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

// 参数组：1 = 文献参数（Liu et al. PRB 88, 085433，拟合的是 DFT 能带）
//         2 = 本工作拟合到 GW 能带的参数（见 data/gw_params.txt）
#ifndef PARAM_SET
#define PARAM_SET 1
#endif

// 能带路径：0 = M-Γ-K-M'-K'（对 VASP），1 = Γ-K-M-Γ（对 wannier90 的 GW 能带）
#ifndef PATH_MODE
#define PATH_MODE 0
#endif

// 模型输出能带的顺序：0 = 按能量大小排序（本征值直接出来的顺序），
//                     1 = 拼成物理能带（两条导带交叉处不互换身份，见 src/bandtrack.h）
// 阶段二要跟拼接过的 GW 参考比、还要拟合，两边必须同一口径，所以 GW 路径默认开；
// 阶段一只是对 VASP 重画一遍师兄的图，两边都按原始排序输出，保持不变。
#ifndef TRACK_BANDS
#define TRACK_BANDS (PATH_MODE == 1)
#endif

#include "params.h"     // 参数分支依赖 NN_MODEL / NNN_MODEL

#define PI 3.14159265358979323846

#endif
