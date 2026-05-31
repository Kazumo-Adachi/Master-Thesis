/* 鈴木式-侵食 */

/* 1. コメントとヘッダファイル */
/* 単位系は(cm,g,sec で) */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)


/* ========= ケース設定（ここだけ変えれば別ケースに対応） ========= */
#define FLUME_ID    1       /* 水路形状: 1=変STD, 2=変GEN, 3=曲 */
#define Q_ID        4       /* 流量: 1=Q_STD, 2=Q_HIGH, 3=Q_LOW, 4=Q_VERY */
/* ================================================================ */


/* ========= マクロ文字列化 ========= */
#define FLUME_ID_STR STR(FLUME_ID)
#define Q_ID_STR     STR(Q_ID)


/* ========= 水路形状の定義（FLUME_IDで切り替え） ========= */
#if   FLUME_ID == 1          /* 侵食_変勾配水路STD */
    #define TOP              147.0337684
    #define FLUME_NAME       "変STD"
    #define FLUME_FILE       "水路形状/侵食水路/E_変STD.txt"

#elif FLUME_ID == 2          /* 侵食_変勾配水路GEN */
    #define TOP              104.4483347
    #define FLUME_NAME       "変GEN"
    #define FLUME_FILE       "水路形状/侵食水路/E_変GEN.txt"

#elif FLUME_ID == 3          /* 侵食_曲線水路 */
    #define TOP              146.859271
    #define FLUME_NAME       "曲"
    #define FLUME_FILE       "水路形状/侵食水路/E_曲.txt"

#else
    #error "Unsupported FLUME_ID (1〜3のどれかを指定してください)"
#endif


/* ========= 流量条件の定義（Q_IDで切り替え） ========= */
#if   Q_ID == 1  /* Q_STD: 標準流量 */
    #define GRAIN_SIZE       0.20    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.20    /* 停止水深 stop [cm] */
    #define M_IN             150.0   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.0331  /* 給水濃度 co [-] */
    #define INFLOW_TIME      20      /* 給水している時間ステップ数 */
    #define GRAIN_NAME       "d020"
    #define Q_NAME           "Q150"

#elif Q_ID == 2  /* Q_HIGH: 大流量 */
    #define GRAIN_SIZE       0.20    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.20    /* 停止水深 stop [cm] */
    #define M_IN             300.0   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.0331  /* 給水濃度 co [-] */
    #define INFLOW_TIME      10      /* 給水している時間ステップ数 */
    #define GRAIN_NAME       "d020"
    #define Q_NAME           "Q300"

#elif Q_ID == 3  /* Q_LOW: 小流量 */
    #define GRAIN_SIZE       0.20    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.20    /* 停止水深 stop [cm] */
    #define M_IN             100.0   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.0331  /* 給水濃度 co [-] */
    #define INFLOW_TIME      30      /* 給水している時間ステップ数 */
    #define GRAIN_NAME       "d020"
    #define Q_NAME           "Q100"

#elif Q_ID == 4  /* Q_VERY: 超大流量 */
    #define GRAIN_SIZE       0.20    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.20    /* 停止水深 stop [cm] */
    #define M_IN             500.0   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.0331  /* 給水濃度 co [-] */
    #define INFLOW_TIME      6       /* 給水している時間ステップ数 */
    #define GRAIN_NAME       "d020"
    #define Q_NAME           "Q500"

#else
    #error "Unsupported Q_ID (1〜4のどれかを指定してください)"
#endif


/* ========= 給水水深 H_IN の決定（FLUME_ID × Q_ID） ========= */
#if   FLUME_ID == 1 && Q_ID == 1
    #define H_IN  2.74
#elif FLUME_ID == 1 && Q_ID == 2
    #define H_IN  3.88
#elif FLUME_ID == 1 && Q_ID == 3
    #define H_IN  2.24
#elif FLUME_ID == 1 && Q_ID == 4
    #define H_IN  5.01

/* ========= 未修正 ========= */
#elif FLUME_ID == 2 && Q_ID == 1
    #define H_IN  2.0
#elif FLUME_ID == 2 && Q_ID == 2
    #define H_IN  2.0
#elif FLUME_ID == 2 && Q_ID == 3
    #define H_IN  2.0
#elif FLUME_ID == 2 && Q_ID == 4
    #define H_IN  2.0

#elif FLUME_ID == 3 && Q_ID == 1
    #define H_IN  2.0
#elif FLUME_ID == 3 && Q_ID == 2
    #define H_IN  2.0
#elif FLUME_ID == 3 && Q_ID == 3
    #define H_IN  2.0
#elif FLUME_ID == 3 && Q_ID == 4
    #define H_IN  2.0

#else
    #error "Invalid FLUME_ID or Q_ID combination"
#endif


/* ========= 入出力パス定義 ========= */
#define INPUT_DIR        "txtファイル/"
#define OUTPUT_DIR       "鈴木侵食/SE_" GRAIN_NAME "_" FLUME_NAME "_" Q_NAME "/"
#define OUTPUT_PREFIX    "鈴木侵食_" GRAIN_NAME "_" FLUME_NAME "_" Q_NAME "_"


/* ========= 計算格子 ========= */
#define dx  5.0                               /* 空間間隔 */
#define dt  0.01F                            /* 時間間隔 */
#define X_GRID 202                            /* x方向グリッド数+1 (+1は境界条件) */
#define T_GRID ((int)(100/2/dt + 0.5))        /* t方向グリッド数/2 */

/* 数値発散防止 */
const double EPS = 1e-10;
const double h_min = 0.01;  /* 乾湿判定の水深閾値 */
const double c_min = 0.01;  /* 清水流判定の濃度閾値 */
#define Record ((int)(0.1/dt + 0.5))

/* CFL条件チェック用の定数 */
#define CFL_TARGET 0.85
#define CFL_CHECK_INTERVAL 100

/* 反復計算の最大回数 */
#define MAX_ITER 50

/* グローバル物性値（関数から参照） */
const double G_ro = 1.0;
const double G_sigma = 2.65;
const double G_e = 0.85;
const double G_tanfais = 0.6745;
const double G_cstar = 0.52;
const double G_g = 980.0;
const double G_kg = 0.0828;
const double G_kar = 0.40;

/* 符号関数 */
static inline double sgn(double x) {
    if (fabs(x) < EPS) return 0.0;
    return (x > 0.0) ? 1.0 : -1.0;
}

/* 二分法失敗カウント（グローバル） */
int g_bisection_fail_count = 0;

/* ========= Ce での fb_e を計算する関数（slope_E使用） ========= */
double calc_fb_e(double Ce, double he, double slope_E, double d, double cs) {
    double kkg_e, kf_e, kkf_e, etao_e, gg_e, r_e;
    double phiw_e, phis_e, phii_e, phi_e;
    double fb_e;

    if (Ce <= c_min) {
        /* 清水流 */
        kf_e = 0.025;
        etao_e = pow(kf_e, 0.5) * d * pow((1 - cs) / cs, 0.333);
        double denom_log = (1 + etao_e / he) * log(1 + he / etao_e) - 1;
        if (fabs(denom_log) < EPS) denom_log = EPS;
        fb_e = pow(G_kar / denom_log, 2);

    } else if (Ce / cs >= 1) {
        /* 土石流状態 */
        kkg_e = G_kg * G_sigma / G_ro * (1 - pow(G_e, 2)) * pow(Ce, 0.333);
        kf_e = 0.16;
        kkf_e = kf_e * pow(1 - Ce, 1.666) / pow(Ce, 0.666);
        fb_e = 6.25 * (kkg_e + kkf_e) * pow(d / he, 2);

    } else {
        /* 掃流状集合流動 */
        kkg_e = G_kg * G_sigma / G_ro * (1 - pow(G_e, 2)) * pow(Ce, 0.333);
        kf_e = 0.16;
        kkf_e = kf_e * pow(1 - Ce, 1.666) / pow(Ce, 0.666);
        etao_e = pow(kf_e, 0.5) * d * pow((1 - cs) / cs, 0.333);
        gg_e = (G_sigma / G_ro - 1) * Ce + 1;

        double log_arg = (etao_e / he + 1 - Ce / cs) * he / etao_e;
        if (log_arg <= 0) log_arg = EPS;

        phiw_e = pow(1 - Ce / cs, 0.5) * (etao_e / he + 1 - Ce / cs)
                 * log(log_arg) * d / he / G_kar
                 - pow(1 - Ce / cs, 1.5) * d / he / G_kar;

        r_e = 1 + cs / Ce * ((1 - pow(Ce / G_cstar, 0.2)) * gg_e - 1);

        if (fabs(r_e) <= 0.00001) {
            phis_e = pow(Ce / cs, 2) * pow(1 - Ce / cs, 0.5) / 2 / pow(kkf_e + kkg_e, 0.5);
            phii_e = Ce / cs * pow(1 - Ce / cs, 1.5) / pow(kkf_e + kkg_e, 0.5);
        } else {
            phis_e = (pow(1 - Ce / cs, 2.5) - pow(1 + (r_e - 1) * Ce / cs, 1.5) * (1 - (3 * r_e / 2 + 1) * Ce / cs))
                     * 4 / 15 / pow(kkg_e + kkf_e, 0.5) / pow(r_e, 2);
            phii_e = (pow(1 + (r_e - 1) * Ce / cs, 1.5) - pow(1 - Ce / cs, 1.5)) * (1 - Ce / cs)
                     * 2 / 3 / r_e / pow(kkg_e + kkf_e, 0.5);
        }
        phi_e = phiw_e + phis_e + phii_e;

        if (phi_e <= EPS) {
            fb_e = 6.25 * (kkg_e + kkf_e) * pow(d / he, 2);
        } else {
            fb_e = (1 - pow(Ce / G_cstar, 0.2)) * gg_e * pow(d / he / phi_e, 2);
        }
    }
    return fb_e;
}

/* ========= 式⑬の左辺 - 右辺 を計算（f(he) = 0 を解く） ========= */
/* (Ce*(σ-ρ)+ρ)*g*he*sin(slope_E) - tauy_e - ρ*fb_e*(M/he)² = 0 */
double eval_he_equation(double Ce, double he, double M, double slope_E, double d, double cs) {
    double fb_e = calc_fb_e(Ce, he, slope_E, d, cs);

    double ro_m_e = Ce * (G_sigma - G_ro) + G_ro;

    double tauy_e;
    if (Ce <= c_min) {
        tauy_e = 0.0;  // 清水流：降伏応力なし
    } else {
        double alpha_e = pow(Ce / G_cstar, 0.2);
        tauy_e = alpha_e * (G_sigma - G_ro) * Ce * G_g * he * cos(slope_E) * G_tanfais;
    }

    double lhs = ro_m_e * G_g * he * sin(slope_E);
    double rhs = tauy_e + G_ro * fb_e * M * M / (he * he);

    return lhs - rhs;
}

/* ========= 平衡水深 he を二分法で求める関数 ========= */
double calc_he_bisection(double Ce, double M, double slope_E, double d, double cs, double he_low_init, double he_high_init) {
    double he_low = he_low_init;
    double he_high = he_high_init;

    double f_low = eval_he_equation(Ce, he_low, M, slope_E, d, cs);
    double f_high = eval_he_equation(Ce, he_high, M, slope_E, d, cs);

    if (f_low * f_high > 0) {
        g_bisection_fail_count++;
        return -1.0;
    }

    for (int iter = 0; iter < MAX_ITER; iter++) {
        double he_mid = (he_low + he_high) / 2.0;
        double f_mid = eval_he_equation(Ce, he_mid, M, slope_E, d, cs);

        if (fabs(f_mid) < 1e-8 || (he_high - he_low) < 1e-8) {
            return he_mid;
        }

        if (f_low * f_mid < 0) {
            he_high = he_mid;
        } else {
            he_low = he_mid;
            f_low = f_mid;
        }
    }

    g_bisection_fail_count++;
    return -1.0;
}


/* ========= プログラムのスタート ========= */
int main(void) {
    /* ファイルの宣言（変数） */
    FILE *fhi, *fho;
    FILE *fzi, *fzo;
    FILE *fMi, *fMo;
    FILE *fui, *fuo;
    FILE *fci, *fco;
    FILE *fEi, *fEo;
    FILE *ftaui, *ftauo;

    /* ファイルの宣言（水路形状） */
    FILE *fflume;

    /* ファイルの宣言（その他） */
    FILE *fslopeo;
    FILE *fPexo;
    FILE *fPmaxo;
    FILE *fPmino;

    /* 離散化の添字 */
    int i, n;

    /* 配列の宣言（リープフロッグ） */
    double z1[X_GRID + 1] = { 0 };
    double h1[X_GRID + 1] = { 0 };
    double M1[X_GRID + 1] = { 0 };
    double u1[X_GRID + 1] = { 0 };
    double c1[X_GRID + 1] = { 0 };

    double z2[X_GRID + 1] = { 0 };
    double h2[X_GRID + 1] = { 0 };
    double M2[X_GRID + 1] = { 0 };
    double u2[X_GRID + 1] = { 0 };
    double c2[X_GRID + 1] = { 0 };

    double E[X_GRID + 1] = { 0 };
    double tau[X_GRID + 1] = { 0 };

    double ct[X_GRID + 1] = { 0 };
    double phi[X_GRID + 1] = { 0 };
    double phis[X_GRID + 1] = { 0 };
    double phiw[X_GRID + 1] = { 0 };
    double phii[X_GRID + 1] = { 0 };

    double flume[X_GRID + 1] = { 0 };
    double slope_arr[X_GRID + 1] = { 0 };

    /* 浅水方程式のパラメータ */
    double d;
    double Mo, ho, co, uo;
    double stop;

    /* 構成則のパラメータ */
    double tauy, kkg, kkf, gg, etao, fb, r;
    double cs;

    /* 物性・定数（ローカル） */
    double ro = 1.0;
    double sigma = 2.65;
    double e = 0.85;
    double tanfais = 0.6745;
    double cstar = 0.52;
    double g = 980.0;
    double beta = 1.00;
    double kg = 0.0828;
    double kf;
    double kar = 0.40;

    double T = 0.56;

    cs = cstar / 2.0;

    /* 供給条件の設定 */
    d = GRAIN_SIZE;
    stop = STOP_DEPTH;
    Mo = M_IN;
    co = C_IN;
    ho = H_IN;
    uo = Mo / ho;

    /* CFL条件チェック用の変数 */
    double max_char_speed = 0.0;
    double CFL_current = 0.0;
    double CFL_max_overall = 0.0;
    int CFL_violation_count = 0;
    double char_speed_plus, char_speed_minus, char_speed_local;
    double discriminant;

    /* ファイルオープン（入力） */
    fhi   = fopen(INPUT_DIR "hi.txt", "r");
    fzi   = fopen(FLUME_FILE, "r");
    fMi   = fopen(INPUT_DIR "Mi.txt", "r");
    fui   = fopen(INPUT_DIR "ui.txt", "r");
    fci   = fopen(INPUT_DIR "ci.txt", "r");
    fEi   = fopen(INPUT_DIR "Ei.txt", "r");
    ftaui = fopen(INPUT_DIR "taui.txt", "r");
    fflume = fopen(FLUME_FILE, "r");

    if (fhi == NULL || fzi == NULL || fMi == NULL || fui == NULL ||
        fci == NULL || fEi == NULL || ftaui == NULL || fflume == NULL) {
        printf("入力ファイルエラー\n");
        printf("  共通入力: %s\n", INPUT_DIR);
        printf("  水路形状: %s\n", FLUME_FILE);
        exit(1);
    }

    /* ファイルオープン（出力） */
    fho    = fopen(OUTPUT_DIR OUTPUT_PREFIX "h.txt", "w");
    fzo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "z.txt", "w");
    fMo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "M.txt", "w");
    fuo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "u.txt", "w");
    fco    = fopen(OUTPUT_DIR OUTPUT_PREFIX "c.txt", "w");
    fEo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "E.txt", "w");
    ftauo  = fopen(OUTPUT_DIR OUTPUT_PREFIX "tau.txt", "w");
    fslopeo = fopen(OUTPUT_DIR OUTPUT_PREFIX "slope.txt", "w");
    fPexo  = fopen(OUTPUT_DIR OUTPUT_PREFIX "Pex.txt", "w");
    fPmaxo = fopen(OUTPUT_DIR OUTPUT_PREFIX "Pmax.txt", "w");
    fPmino = fopen(OUTPUT_DIR OUTPUT_PREFIX "Pmin.txt", "w");

    if (fho == NULL || fzo == NULL || fMo == NULL || fuo == NULL ||
        fco == NULL || fEo == NULL || ftauo == NULL ||
        fslopeo == NULL || fPexo == NULL || fPmaxo == NULL) {
        printf("出力ファイルエラー\n");
        printf("  出力先: %s\n", OUTPUT_DIR);
        printf("  フォルダが存在するか確認してください\n");
        exit(1);
    }

    printf("========== 鈴木侵食モデル ==========\n");
    printf("  粒径:     %s (%.2f cm)\n", GRAIN_NAME, GRAIN_SIZE);
    printf("  水路形状: %s\n", FLUME_NAME);
    printf("  流量:     %s (Mo=%.1f cm^2/s)\n", Q_NAME, M_IN);
    printf("  出力先:   %s\n", OUTPUT_DIR);
    printf("====================================\n");

    printf("\n===== CFL条件チェック設定 =====\n");
    printf("目標クーラン数: %.2f\n", CFL_TARGET);
    printf("dx = %.2f [cm], dt = %.4f [s]\n", dx, dt);
    printf("beta = %.2f, g = %.1f [cm/s^2]\n", beta, g);
    printf("================================\n\n");

    /* 端点の初期条件 */
    h1[0] = ho;
    h1[X_GRID] = 0;
    M1[0] = Mo;
    M1[X_GRID] = 0;
    u1[0] = uo;
    u1[X_GRID] = 0;
    c1[0] = co;
    c1[X_GRID] = 0;
    ct[0] = co;
    ct[X_GRID] = 0;
    z1[0] = TOP;
    z1[X_GRID] = -1.29;

    h2[0] = ho;
    h2[X_GRID] = 0;
    M2[0] = Mo;
    M2[X_GRID] = 0;
    u2[0] = uo;
    u2[X_GRID] = 0;
    c2[0] = co;
    c2[X_GRID] = 0;
    z2[0] = TOP;
    z2[X_GRID] = -1.29;

    E[0] = 0;
    E[X_GRID] = 0;
    tau[0] = 0;
    tau[X_GRID] = 0;
    flume[0] = TOP;
    flume[X_GRID] = -1.29;

    /* ファイルから内部点の初期条件を読み込む */
    for (i = 1; i < X_GRID; i++) {
        fscanf(fhi, "%lf", &h1[i]); fscanf(fzi, "%lf", &z1[i]);
        fscanf(fMi, "%lf", &M1[i]); fscanf(fui, "%lf", &u1[i]);
        fscanf(fci, "%lf", &c1[i]); fscanf(fEi, "%lf", &E[i]);
        fscanf(ftaui, "%lf", &tau[i]); fscanf(fflume, "%lf", &flume[i]);
    }

    /* 流速係数の初期計算 */
    for (i = 0; i < X_GRID; i++) {
        if (h1[i] > h_min && c1[i] > c_min) {
            kkg = kg * sigma / ro * (1 - pow(e, 2)) * pow(c1[i], 0.333);
            kf = 0.16;
            kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
            etao = pow(kf, 0.5) * d * pow((1 - cs) / cs, 0.333);
            gg = (sigma / ro - 1) * c1[i] + 1;

            if (c1[i] / cs < 1) {
                phiw[i] = pow(1 - c1[i] / cs, 0.5) * (etao / h1[i] + 1 - c1[i] / cs)
                          * log((etao / h1[i] + 1 - c1[i] / cs) * h1[i] / etao) * d / h1[i] / kar
                          - pow(1 - c1[i] / cs, 1.5) * d / h1[i] / kar;
                r = 1 + cs / c1[i] * ((1 - pow(c1[i] / cstar, 0.2)) * gg - 1);
                if (fabs(r) <= 0.00001) {
                    phis[i] = pow(c1[i] / cs, 2) * pow(1 - c1[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5);
                    phii[i] = c1[i] / cs * pow(1 - c1[i] / cs, 1.5) / pow(kkf + kkg, 0.5);
                } else {
                    phis[i] = (pow(1 - c1[i] / cs, 2.5) - pow(1 + (r - 1) * c1[i] / cs, 1.5) * (1 - (3 * r / 2 + 1) * c1[i] / cs))
                              * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2);
                    phii[i] = (pow(1 + (r - 1) * c1[i] / cs, 1.5) - pow(1 - c1[i] / cs, 1.5)) * (1 - c1[i] / cs)
                              * 2 / 3 / r / pow(kkg + kkf, 0.5);
                }
                phi[i] = phiw[i] + phis[i] + phii[i];
            }
        } else {
            phiw[i] = 0; phis[i] = 0; phii[i] = 0; phi[i] = 0;
        }
    }

    /* ======================================================================= */
    /* メインループ */
    /* ======================================================================= */
    for (n = 0; n < T_GRID; n++) {

        /* 給水条件 */
        if (n < INFLOW_TIME / (2 * dt)) {
            h1[0] = ho; h2[0] = ho; M1[0] = Mo; M2[0] = Mo;
            u1[0] = uo; u2[0] = uo; c1[0] = co; c2[0] = co; ct[0] = ro*tan(2*M_PI/180)/(sigma-ro)/(tanfais-tan(2*M_PI/180));
        } else {
            h1[0] = 0; h2[0] = 0; M1[0] = 0; M2[0] = 0;
            u1[0] = 0; u2[0] = 0; c1[0] = 0; c2[0] = 0; ct[0] = 0;
        }

        /* CFL条件チェック */
        max_char_speed = 0.0;
        for (i = 1; i < X_GRID; i++) {
            if (h1[i] > h_min) {
                discriminant = g * h1[i] + beta * (beta - 1.0) * u1[i] * u1[i];
                if (discriminant > 0) {
                    double sqrt_disc = sqrt(discriminant);
                    char_speed_plus = beta * u1[i] + sqrt_disc;
                    char_speed_minus = beta * u1[i] - sqrt_disc;
                    char_speed_local = fmax(fabs(char_speed_plus), fabs(char_speed_minus));
                    if (char_speed_local > max_char_speed) max_char_speed = char_speed_local;
                }
            }
        }
        if (max_char_speed > EPS) {
            CFL_current = max_char_speed * dt / dx;
            if (CFL_current > CFL_max_overall) CFL_max_overall = CFL_current;
            if (CFL_current > CFL_TARGET) {
                CFL_violation_count++;
                if (CFL_violation_count == 1 || n % CFL_CHECK_INTERVAL == 0) {
                    printf("[WARNING] CFL violation at n=%d: CFL=%.4f > %.2f\n", n, CFL_current, CFL_TARGET);
                }
            }
        }

        /* ===== Step.A1 連続式 ===== */
        for (i = 1; i < X_GRID; i++) {
            double M_diff;
            if (u1[i] >= 0) {
                M_diff = M1[i] - M1[i-1];
            } else {
                M_diff = M1[i+1] - M1[i];
            }

            h2[i] = h1[i] + dt * (E[i] - M_diff / dx);

            if (h2[i] <= h_min) {
                h2[i] = 0; c2[i] = 0; z2[i] = z1[i];
            } else {
                if (c1[i] <= c_min) ct[i] = 0;
                else if (c1[i] / cs >= 1 || phis[i] <= EPS || phi[i] <= EPS) ct[i] = c1[i];
                else ct[i] = cs * phis[i] / phi[i];

                double cM_diff;
                if (u1[i] >= 0) {
                    cM_diff = ct[i] * M1[i] - ct[i-1] * M1[i-1];
                } else {
                    cM_diff = ct[i+1] * M1[i+1] - ct[i] * M1[i];
                }

                c2[i] = (c1[i] * h1[i] + dt * (E[i] * cstar - cM_diff / dx)) / h2[i];
                z2[i] = z1[i] - E[i] * dt;
            }
        }

        /* 境界値の更新 */
        z2[X_GRID] = 2*z2[X_GRID-1] - z2[X_GRID-2];
        h2[X_GRID] = h2[X_GRID-1];
        c2[X_GRID] = c2[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];

        /* ===== Step.B1+C1 構成則+運動方程式 ===== */
        for (i = 1; i < X_GRID; i++) {
            phi[i] = phis[i] = phiw[i] = phii[i] = 0.0;

            double slope = atan((z2[i] - z2[i + 1]) / dx);

            if (h2[i] < stop) {
                tau[i] = 0; M2[i] = 0; u2[i] = 0; E[i] = 0;
            } else {
                if (c2[i] <= c_min) {
                    kf = 0.025;
                    etao = pow(kf, 0.5) * d * pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h2[i]) * log(1 + h2[i] / etao) - 1), 2);
                    tauy = 0;
                } else if (c2[i] / cs >= 1) {
                    kkg = kg * sigma / ro * (1 - pow(e, 2)) * pow(c2[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c2[i], 1.666) / pow(c2[i], 0.666);
                    fb = 6.25 * (kkg + kkf) * pow(d / h2[i], 2);
                    tauy = pow(c2[i] / cstar, 0.2) * (sigma - ro) * c2[i] * g * h2[i] * cos(slope) * tanfais;
                } else {
                    kkg = kg * sigma / ro * (1 - pow(e, 2)) * pow(c2[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c2[i], 1.666) / pow(c2[i], 0.666);
                    etao = pow(kf, 0.5) * d * pow((1 - cs) / cs, 0.333);
                    gg = (sigma / ro - 1) * c2[i] + 1;
                    phiw[i] = pow(1 - c2[i] / cs, 0.5) * (etao / h2[i] + 1 - c2[i] / cs)
                              * log((etao / h2[i] + 1 - c2[i] / cs) * h2[i] / etao) * d / h2[i] / kar
                              - pow(1 - c2[i] / cs, 1.5) * d / h2[i] / kar;
                    r = 1 + cs / c2[i] * ((1 - pow(c2[i] / cstar, 0.2)) * gg - 1);
                    if (fabs(r) <= 0.00001) {
                        phis[i] = pow(c2[i] / cs, 2) * pow(1 - c2[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5);
                        phii[i] = c2[i] / cs * pow(1 - c2[i] / cs, 1.5) / pow(kkf + kkg, 0.5);
                    } else {
                        phis[i] = (pow(1 - c2[i] / cs, 2.5) - pow(1 + (r - 1) * c2[i] / cs, 1.5) * (1 - (3 * r / 2 + 1) * c2[i] / cs))
                                  * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2);
                        phii[i] = (pow(1 + (r - 1) * c2[i] / cs, 1.5) - pow(1 - c2[i] / cs, 1.5)) * (1 - c2[i] / cs)
                                  * 2 / 3 / r / pow(kkg + kkf, 0.5);
                    }
                    phi[i] = phiw[i] + phis[i] + phii[i];
                    if (phi[i] <= EPS) fb = 6.25 * (kkg + kkf) * pow(d / h2[i], 2);
                    else fb = (1 - pow(c2[i] / cstar, 0.2)) * gg * pow(d / h2[i] / phi[i], 2);
                    tauy = pow(c2[i] / cstar, 0.2) * (sigma - ro) * c2[i] * g * h2[i] * cos(slope) * tanfais;
                }

                tau[i] = tauy * sgn(u1[i]) + ro * fb * u1[i] * fabs(u1[i]);

                double ro_m = fmax(ro + (sigma - ro) * c2[i], ro);

                double uM_diff;
                if (u1[i] >= 0) {
                    uM_diff = u1[i] * M1[i] - u1[i-1] * M1[i-1];
                } else {
                    uM_diff = u1[i+1] * M1[i+1] - u1[i] * M1[i];
                }

                M2[i] = M1[i] - dt * (beta * uM_diff / dx + g * h2[i] * (h2[i+1] + z2[i+1] - h2[i] - z2[i]) / dx + tau[i] / ro_m);
                u2[i] = M2[i] / h2[i];

                /* ===== 鈴木式：侵食速度 E の計算 ===== */
                double slope_E = slope;
                if (slope_E < 0) {
                    slope_E = atan(((z2[i]+h2[i])-(z2[i+1]+h2[i+1]))/dx);
                }

                double tan_slope = tan(slope_E);
                double Ce = ro * tan_slope / (sigma - ro) / (tanfais - tan_slope);

                if (Ce < 0) Ce = 0;
                if (Ce > cstar) Ce = cstar - EPS;

                double he = calc_he_bisection(Ce, fabs(M2[i]), slope_E, d, cs, 0.001, 20.0);
                

                if (he < 0 || isnan(he) || isinf(he)) {
                    E[i] = -h2[i] * c2[i] / T / cstar;
                } else {
                    E[i] = (he * Ce - h2[i] * c2[i]) / T / cstar;
                }

                if (c2[i] <= c_min && E[i] < 0) E[i] = 0;

                /* Eの下限を設定 */
                double E_min1 = -h2[i] * c2[i] / dt / cstar;
                double E_min2 = -h2[i] / dt;
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                slope_arr[i] = slope;
            }
        }

        /* 境界値の更新 */
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];

        /* ===== Step.A2 連続式 ===== */
        for (i = 1; i < X_GRID; i++) {
            double M_diff;
            if (u2[i] >= 0) {
                M_diff = M2[i] - M2[i-1];
            } else {
                M_diff = M2[i+1] - M2[i];
            }

            h1[i] = h2[i] + dt * (E[i] - M_diff / dx);

            if (h1[i] <= h_min) {
                h1[i] = 0; c1[i] = 0; z1[i] = z2[i];
            } else {
                if (c2[i] <= c_min) ct[i] = 0;
                else if (c2[i] / cs >= 1 || phis[i] <= EPS || phi[i] <= EPS) ct[i] = c2[i];
                else ct[i] = cs * phis[i] / phi[i];

                double cM_diff;
                if (u2[i] >= 0) {
                    cM_diff = ct[i] * M2[i] - ct[i-1] * M2[i-1];
                } else {
                    cM_diff = ct[i+1] * M2[i+1] - ct[i] * M2[i];
                }

                c1[i] = (c2[i] * h2[i] + dt * (E[i] * cstar - cM_diff / dx)) / h1[i];
                z1[i] = z2[i] - E[i] * dt;
            }

            if (n % Record == Record - 1) { fprintf(fho, "%f ", h1[i]); }
            if (n % Record == Record - 1) { fprintf(fco, "%f ", c1[i]); }
            if (n % Record == Record - 1) { fprintf(fzo, "%f ", z1[i]); }
        }

        /* 境界値の更新 */
        z1[X_GRID] = 2*z1[X_GRID-1] - z1[X_GRID-2];
        h1[X_GRID] = h1[X_GRID-1];
        c1[X_GRID] = c1[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];

        /* ===== Step.B2+C2 ===== */
        for (i = 1; i < X_GRID; i++) {
            phi[i] = phis[i] = phiw[i] = phii[i] = 0.0;

            double slope = atan((z1[i] - z1[i + 1]) / dx);

            if (h1[i] < stop) {
                tau[i] = 0; M1[i] = 0; u1[i] = 0; E[i] = 0;
            } else {
                if (c1[i] <= c_min) {
                    kf = 0.025;
                    etao = pow(kf, 0.5) * d * pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h1[i]) * log(1 + h1[i] / etao) - 1), 2);
                    tauy = 0;
                } else if (c1[i] / cs >= 1) {
                    kkg = kg * sigma / ro * (1 - pow(e, 2)) * pow(c1[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
                    fb = 6.25 * (kkg + kkf) * pow(d / h1[i], 2);
                    tauy = pow(c1[i] / cstar, 0.2) * (sigma - ro) * c1[i] * g * h1[i] * cos(slope) * tanfais;
                } else {
                    kkg = kg * sigma / ro * (1 - pow(e, 2)) * pow(c1[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
                    etao = pow(kf, 0.5) * d * pow((1 - cs) / cs, 0.333);
                    gg = (sigma / ro - 1) * c1[i] + 1;
                    phiw[i] = pow(1 - c1[i] / cs, 0.5) * (etao / h1[i] + 1 - c1[i] / cs)
                              * log((etao / h1[i] + 1 - c1[i] / cs) * h1[i] / etao) * d / h1[i] / kar
                              - pow(1 - c1[i] / cs, 1.5) * d / h1[i] / kar;
                    r = 1 + cs / c1[i] * ((1 - pow(c1[i] / cstar, 0.2)) * gg - 1);
                    if (fabs(r) <= 0.00001) {
                        phis[i] = pow(c1[i] / cs, 2) * pow(1 - c1[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5);
                        phii[i] = c1[i] / cs * pow(1 - c1[i] / cs, 1.5) / pow(kkf + kkg, 0.5);
                    } else {
                        phis[i] = (pow(1 - c1[i] / cs, 2.5) - pow(1 + (r - 1) * c1[i] / cs, 1.5) * (1 - (3 * r / 2 + 1) * c1[i] / cs))
                                  * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2);
                        phii[i] = (pow(1 + (r - 1) * c1[i] / cs, 1.5) - pow(1 - c1[i] / cs, 1.5)) * (1 - c1[i] / cs)
                                  * 2 / 3 / r / pow(kkg + kkf, 0.5);
                    }
                    phi[i] = phiw[i] + phis[i] + phii[i];
                    if (phi[i] <= EPS) fb = 6.25 * (kkg + kkf) * pow(d / h1[i], 2);
                    else fb = (1 - pow(c1[i] / cstar, 0.2)) * gg * pow(d / h1[i] / phi[i], 2);
                    tauy = pow(c1[i] / cstar, 0.2) * (sigma - ro) * c1[i] * g * h1[i] * cos(slope) * tanfais;
                }

                tau[i] = tauy * sgn(u2[i]) + ro * fb * u2[i] * fabs(u2[i]);

                double ro_m = fmax(ro + (sigma - ro) * c1[i], ro);

                double uM_diff;
                if (u2[i] >= 0) {
                    uM_diff = u2[i] * M2[i] - u2[i-1] * M2[i-1];
                } else {
                    uM_diff = u2[i+1] * M2[i+1] - u2[i] * M2[i];
                }

                M1[i] = M2[i] - dt * (beta * uM_diff / dx + g * h1[i] * (h1[i+1] + z1[i+1] - h1[i] - z1[i]) / dx + tau[i] / ro_m);
                u1[i] = M1[i] / h1[i];

                /* ===== 鈴木式：侵食速度 E の計算 ===== */
                double slope_E = slope;
                if (slope_E < 0) {
                    slope_E = atan(((z1[i]+h1[i])-(z1[i+1]+h1[i+1]))/dx);
                }

                double tan_slope = tan(slope_E);
                double Ce = ro * tan_slope / (sigma - ro) / (tanfais - tan_slope);

                if (Ce < 0) Ce = 0;
                if (Ce > cstar) Ce = cstar - EPS;

                double he = calc_he_bisection(Ce, fabs(M1[i]), slope_E, d, cs, 0.001, 20.0);

                if (he < 0 || isnan(he) || isinf(he)) {
                    E[i] = -h1[i] * c1[i] / T / cstar;
                } else {
                    E[i] = (he * Ce - h1[i] * c1[i]) / T / cstar;
                }

                if (c1[i] <= c_min && E[i] < 0) E[i] = 0;

                /* Eの下限を設定 */
                double E_min1 = -h1[i] * c1[i] / dt / cstar;
                double E_min2 = -h1[i] / dt;
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                slope_arr[i] = slope;
            }

            if (n % Record == Record - 1) { fprintf(ftauo, "%f ", tau[i]); }
            if (n % Record == Record - 1) { fprintf(fMo,   "%f ", M1[i]); }
            if (n % Record == Record - 1) { fprintf(fuo,   "%f ", u1[i]); }
            if (n % Record == Record - 1) { fprintf(fEo,   "%f ", E[i]); }
            if (n % Record == Record - 1) { fprintf(fslopeo, "%f ", slope_arr[i]); }
        }

        /* 境界値の更新 */
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];

        /* 改行出力 */
        if (n % Record == Record - 1) {
            fprintf(fho, "\n");
            fprintf(fco, "\n");
            fprintf(fzo, "\n");
            fprintf(ftauo, "\n");
            fprintf(fMo, "\n");
            fprintf(fuo, "\n");
            fprintf(fEo, "\n");
            fprintf(fPexo, "\n");
            fprintf(fPmaxo, "\n");
            fprintf(fPmino, "\n");
            fprintf(fslopeo, "\n");
        }
    }

    /* ======================================================================= */
    /* CFL条件の最終サマリー */
    /* ======================================================================= */
    printf("\n===== CFL条件チェック結果 =====\n");
    printf("計算全体での最大クーラン数: %.4f\n", CFL_max_overall);
    printf("目標クーラン数: %.2f\n", CFL_TARGET);
    if (CFL_violation_count > 0) {
        printf("[WARNING] CFL違反が %d 回発生しました\n", CFL_violation_count);
        double dt_recommended = CFL_TARGET * dx / (CFL_max_overall * dx / dt);
        printf("安定性のため、dt <= %.6f [s] を推奨します\n", dt_recommended);
    } else {
        printf("[OK] CFL条件は全ステップで満たされています\n");
    }
    printf("================================\n");

    /* 二分法失敗カウントの表示 */
    printf("\n===== 二分法（he計算）結果 =====\n");
    if (g_bisection_fail_count > 0) {
        printf("[WARNING] 二分法で解が見つからなかった回数: %d\n", g_bisection_fail_count);
    } else {
        printf("[OK] 二分法は全て成功しました\n");
    }
    printf("================================\n");

    fclose(fhi); fclose(fho); fclose(fzi); fclose(fzo);
    fclose(fMi); fclose(fMo); fclose(fui); fclose(fuo);
    fclose(fci); fclose(fco); fclose(fEi); fclose(fEo);
    fclose(ftaui); fclose(ftauo); fclose(fPexo);
    fclose(fPmaxo); fclose(fPmino); fclose(fslopeo); fclose(fflume);

    printf("----- 計算完了（鈴木侵食 %s_%s_%s） -----\n", GRAIN_NAME, FLUME_NAME, Q_NAME);
    return 0;
}
