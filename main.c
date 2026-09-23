#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "config.h"
#include "tb.h"
#include "herm3.h"
#include "kpath.h"
#include "bandtrack.h"

// 输出表的列数（能带条数）
#if SOC_MODEL
#define TB_NBAND 6
#else
#define TB_NBAND 3
#endif

// 输入升序排列的能带 eig，以及能带总数 nband
// 返回价带顶能量 (VBM)
static double find_vbm_from_sorted_eigenvalues(double *eig, int nband) {
    // 找到相邻能级差最大的位置，作为带隙
    double max_gap = 0.0;
    int gap_idx = -1;
    for (int i = 0; i < nband-1; i++) {
        double gap = eig[i+1] - eig[i];
        if (gap > max_gap) {
            max_gap = gap;
            gap_idx = i;
        }
    }
    // gap_idx 是价带顶的索引，价带是从 0 到 gap_idx (包括)
    if (gap_idx < 0) return eig[0]; // 无带隙时取最低能带？不太可能
    return eig[gap_idx]; // VBM 是价带的最大值
}

int main(void)
{
    TB_Hop hop;
    tb_init(&hop);

    size_t num_segments = 0;
    const Segment *path = kpath_segments(&num_segments);

    double total_len = 0.0;
    double seg_len[num_segments];
    for (size_t s = 0; s < num_segments; s++) {
        seg_len[s] = kpath_seg_len(&path[s]);
        total_len += seg_len[s];
    }

#if SOC_MODEL
#if SOC_SPIN
    printf("# k-path    Up1    Up2    Up3    Dn1    Dn2    Dn3\n");
#else
    printf("# " K_AXIS_LABEL_S "    E1    E2    E3    E4    E5    E6\n");
#endif
#else
    printf("# " K_AXIS_LABEL "    E1 (eV)    E2 (eV)    E3 (eV)\n");
#endif
#if PATH_MODE == 0
    printf("# path: M → Γ → K → M' → K'\n");
    printf("# Symmetry points: M=0, Γ=%.6f, K=%.6f, M'=%.6f, K'=%.6f\n",
            (seg_len[0])/(2*PI), (seg_len[0]+seg_len[1])/(2*PI), (seg_len[0]+seg_len[1]+seg_len[2])/(2*PI), (total_len)/(2*PI));
#else
    printf("# path: Γ → K → M → Γ   (k 点与 wannier90_band.kpt 的 285 个点重合)\n");
    printf("# Symmetry points (Å^-1): Γ=0, K=%.6f, M=%.6f, Γ=%.6f\n",
            (seg_len[0])/A_ANG, (seg_len[0]+seg_len[1])/A_ANG, (total_len)/A_ANG);
#endif

    double cumulative = 0.0;
    double complex Hgamma[3][3];
    double kgx = 0.0, kgy = 0.0; // Γ点分数坐标
    tb_build_Hk(&hop, kgx, kgy, Hgamma);
    printf("\n");
    printf("# H at Gamma (0,0):\n");
    for (int i=0;i<3;i++) {
        printf("# ");
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(Hgamma[i][j]), cimag(Hgamma[i][j]));
        printf("\n");
    }

    double complex HK[3][3];
    double kKx, kKy;
    // 注意 tb_build_Hk 收的是笛卡尔 k，不能直接传分数坐标
    tb_frac_to_cart_k(2.0/3.0, 1.0/3.0, &kKx, &kKy);
    tb_build_Hk(&hop, kKx, kKy, HK);
    printf("\n");
    printf("# H at K (2/3,1/3):\n");
    for (int i=0;i<3;i++) {
        printf("# ");
        for (int j=0;j<3;j++) printf("(%f%+fi) ", creal(HK[i][j]), cimag(HK[i][j]));
        printf("\n");
    }

#if SOC_MODEL
    #if SOC_SPIN
    {
        double kx, ky;
        tb_frac_to_cart_k(KPATH_GAMMA[0], KPATH_GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        tb_build_Hk(&hop, kx, ky, H_gamma);

        double eig_up[3], eig_dn[3];
        tb_eig_soc_spin(H_gamma, eig_up, eig_dn);

        // 合并排序 6 个值
        double eig_all[6];
        for (int i = 0; i < 3; i++) eig_all[i]   = eig_up[i];
        for (int i = 0; i < 3; i++) eig_all[i+3] = eig_dn[i];
        // 排序
        for (int i = 0; i < 5; i++)
            for (int j = i+1; j < 6; j++)
                if (eig_all[i] > eig_all[j]) {
                    double tmp = eig_all[i];
                    eig_all[i] = eig_all[j];
                    eig_all[j] = tmp;
                }

        double vbm = find_vbm_from_sorted_eigenvalues(eig_all, 6);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
    #else
    {
        double kx, ky;
        tb_frac_to_cart_k(KPATH_GAMMA[0], KPATH_GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        tb_build_Hk(&hop, kx, ky, H_gamma);

        double eig[6];
        tb_eig_soc(H_gamma, eig);  // 已排序升序

        double vbm = find_vbm_from_sorted_eigenvalues(eig, 6);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
    #endif
#else
    {
        double kx, ky;
        tb_frac_to_cart_k(KPATH_GAMMA[0], KPATH_GAMMA[1], &kx, &ky);
        double complex H_gamma[3][3];
        tb_build_Hk(&hop, kx, ky, H_gamma);
        double eig[3];
        herm3_eigvals(H_gamma, eig); // 已排序
        double vbm = find_vbm_from_sorted_eigenvalues(eig, 3);
        printf("\n");
        printf("# TB_Gamma_VBM = %12.6f\n", vbm);
    }
#endif
#if SOC_MODEL
    // 打印 K 点的 6 个 SOC 本征值，并计算 VBM 劈裂
    {
        double kxK, kyK;
        tb_frac_to_cart_k(KPATH_K[0], KPATH_K[1], &kxK, &kyK);
        double complex HK_soc[3][3];
        tb_build_Hk(&hop, kxK, kyK, HK_soc);

        double eig_K[6];
        tb_eig_soc(HK_soc, eig_K);

        printf("\n# SOC eigenvalues at K (2/3, 1/3):\n");
        for (int b = 0; b < 6; b++) {
            printf("#   E[%d] = %12.6f eV\n", b, eig_K[b]);
        }

        // double gap_threshold = 0.5; // 假设带隙 > 0.5 eV
        double vbm1 = -100.0, vbm2 = -100.0;
        for (int b = 0; b < 6; b++) {
            if (eig_K[b] < 0.5) { // 取低于阈值的为价带
                if (eig_K[b] > vbm1) {
                    vbm2 = vbm1;
                    vbm1 = eig_K[b];
                } else if (eig_K[b] > vbm2) {
                    vbm2 = eig_K[b];
                }
            }
        }
        printf("# VBM1 = %12.6f eV, VBM2 = %12.6f eV\n", vbm1, vbm2);
        printf("# VBM SOC splitting = %12.6f eV  (expected 2*lambda = %12.6f eV)\n",
               vbm1 - vbm2, 2.0 * lambda);
    }
#endif
    // 能带表先整条路径算完存下来，最后一起输出。阶段二要先把排序本征值拼回物理能带，
    // 而拼接（斜率外推，见 bandtrack.h）要看前后若干个 k 点，没法边算边打。
    int np = 0;
    for (size_t s = 0; s < num_segments; s++) np += path[s].steps + 1;

    double *tab_x = malloc(sizeof *tab_x * np);
    double *tab_e = malloc(sizeof *tab_e * np * TB_NBAND);
    if (!tab_x || !tab_e) { fprintf(stderr, "[ERROR] out of memory\n"); return 1; }

    int ip = 0;
    for (size_t s = 0; s < num_segments; s++) {
        double start_x = path[s].start[0], start_y = path[s].start[1];
        double end_x   = path[s].end[0],   end_y   = path[s].end[1];
        int steps = path[s].steps;
        for (int i=0; i<=steps; i++) {
            double t = (double)i / steps;
            double fx = start_x + t * (end_x - start_x);
            double fy = start_y + t * (end_y - start_y);
            double kx, ky;
            tb_frac_to_cart_k(fx, fy, &kx, &ky);

            double complex H[3][3];
            tb_build_Hk(&hop, kx, ky, H);

            tab_x[ip] = (cumulative + t * seg_len[s]) / K_AXIS_UNIT;

#if SOC_MODEL
    #if SOC_SPIN
            double eig_up[3], eig_dn[3];
            tb_eig_soc_spin(H, eig_up, eig_dn);
            // 6 列：up1 up2 up3 dn1 dn2 dn3
            for (int b = 0; b < 3; b++) {
                tab_e[b*np + ip]     = eig_up[b];
                tab_e[(3+b)*np + ip] = eig_dn[b];
            }
    #else
            double eig[6];
            tb_eig_soc(H, eig);
            for (int b = 0; b < 6; b++) tab_e[b*np + ip] = eig[b];
    #endif
#else
            double eig[3];
            herm3_eigvals(H, eig);
            for (int b = 0; b < 3; b++) tab_e[b*np + ip] = eig[b];
#endif

            ip++;
        }
        cumulative += seg_len[s];

    }

#if TRACK_BANDS && !SOC_MODEL
    // 模型自己的两条导带在 M-Γ 之间会交叉，而 herm3_eigvals 给的是按能量排序的本征值，
    // 交叉处两条曲线会互换身份。这里拼回物理能带：输出的第 1/2/3 列就是第 1/2/3 条
    // 物理能带，跟 GW 参考（同样拼接过）一一对应。
    // SOC 打开时输出 6 列是两组本征值（或 up/dn 分列），不能放在一起排，阶段二也用不到。
    {
        int anchor = path[0].steps;     // K 点（Γ→K 段的末点），三条带在那里分得最开
        double *trk = malloc(sizeof *trk * np * TB_NBAND);
        if (!trk || band_track(np, TB_NBAND, tab_e, anchor, trk)) {
            fprintf(stderr, "[ERROR] 能带拼接失败\n");
            return 1;
        }
        for (int b = 0; b < TB_NBAND; b++)
            for (int i = 0; i < np; i++) tab_e[b*np + i] = trk[b*np + i];
        free(trk);
    }
#endif

    for (int i = 0; i < np; i++) {
#if SOC_MODEL
    #if SOC_SPIN
        printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                tab_x[i],
                tab_e[0*np+i], tab_e[1*np+i], tab_e[2*np+i],
                tab_e[3*np+i], tab_e[4*np+i], tab_e[5*np+i]);
    #else
        printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                tab_x[i],
                tab_e[0*np+i], tab_e[1*np+i], tab_e[2*np+i],
                tab_e[3*np+i], tab_e[4*np+i], tab_e[5*np+i]);
    #endif
#else
        printf("%12.6f %12.6f %12.6f %12.6f\n",
                tab_x[i], tab_e[0*np+i], tab_e[1*np+i], tab_e[2*np+i]);
#endif
    }

    free(tab_x);
    free(tab_e);
    return 0;
}
