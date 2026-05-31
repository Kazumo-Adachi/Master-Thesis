/* いただいた經隆コード準拠-江頭式 */

/* 1. コメントとヘッダファイル */
/* 単位系は(cm,g,sec で) */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)


/* ========= ケース設定（ここだけ変えれば別ケースに対応） ========= */
#define FLUME_ID    1       /* 水路形状: 1=直線, 2=起伏小, 3=起伏大 */
#define SLOPE_VALUE    30       /* 勾配: 20度, 30度 */
/* ================================================================ */


/* ========= マクロ文字列化 ========= */
#define FLUME_ID_STR STR(FLUME_ID)

//* ========= 水路の端点の定義（SLOPE_VALUEで切り替え） ========= */
#if     SLOPE_VALUE == 20
    #define TOP              8350.283818
    #define BOTTOM           1409.538931
    #define FIXED_BED_FILE   "水路形状/20度水路/固定床.txt"

#elif   SLOPE_VALUE == 30
    #define TOP              11399.53811
    #define BOTTOM           1299.038106
    #define FIXED_BED_FILE   "水路形状/30度水路/固定床.txt"
#else
    #error "Invalid SLOPE_VALUE combination"
#endif


/* ========= 水路形状の定義（FLUME_IDで切り替え） ========= */
#if     FLUME_ID == 1 && SLOPE_VALUE == 20         
    #define FLUME_NAME       "20度_直線"
    #define FLUME_FILE       "水路形状/20度水路/直線.txt"

#elif   FLUME_ID == 2 && SLOPE_VALUE == 20   
    #define FLUME_NAME       "20度_起伏小"
    #define FLUME_FILE       "水路形状/20度水路/起伏小.txt"

#elif   FLUME_ID == 3 && SLOPE_VALUE == 20
    #define FLUME_NAME       "20度_起伏大"
    #define FLUME_FILE       "水路形状/20度水路/起伏大.txt"

#elif   FLUME_ID == 1 && SLOPE_VALUE == 30         
    #define FLUME_NAME       "30度_直線"
    #define FLUME_FILE       "水路形状/30度水路/直線.txt"

#elif   FLUME_ID == 2 && SLOPE_VALUE == 30
    #define FLUME_NAME       "30度_起伏小"
    #define FLUME_FILE       "水路形状/30度水路/起伏小.txt"

#elif   FLUME_ID == 3 && SLOPE_VALUE == 30 
    #define FLUME_NAME       "30度_起伏大"
    #define FLUME_FILE       "水路形状/30度水路/起伏大.txt"

#else
    #error "Invalid FLUME_ID or SLOPE_VALUE combination"
#endif


/* ========= 流量条件の定義 ========= */
#define GRAIN_SIZE       20.0    /* 粒径 d [cm] */
#define STOP_DEPTH       20.0    /* 停止水深 stop [cm] */
#define M_IN             1830.0  /* 給水流量 Mo [cm^2/s]  0.183[m^3/s]/1[m]->0.183[m^2/s]->1830[cm^2/s] */
#define C_IN             0.01     /* 給水濃度 co [-] */
//#define H_IN             100   /* 給水水深 ho [cm] */
#define H_IN  100 * 0.2 * sqrt(0.183)
#define INFLOW_TIME      1200     /* 給水している時間ステップ数 */




/* ========= 入出力パス定義 ========= */
#define INPUT_DIR        "txtファイル/"
#define OUTPUT_DIR       "FIX江頭/" FLUME_NAME "/"
#define OUTPUT_PREFIX    FLUME_NAME "_"


/* ========= 計算格子 ========= */
#define dx  100.0                            /* 空間間隔 */
#define dt  0.01F                            /* 時間間隔 */
#define X_GRID 201                           /* x方向グリッド数+1 (+1は境界条件) */
#define T_GRID 60000                         /* t方向グリッド数/2  総計算時間600s = T_GRID*2*dt */

/* 数値発散防止 */
const double EPS = 1e-10;
const double h_min = EPS;  /* 乾湿判定の水深閾値 */
const double c_min = EPS;  /* 清水流判定の濃度閾値 */
#define Record ((int)(0.1/dt + 0.5))

/* CFL条件チェック用の定数 */
#define CFL_TARGET 0.85           /* 目標クーラン数 */
#define CFL_CHECK_INTERVAL 100    /* CFLチェックの間隔（ステップ数） */

/* 符号関数sgnの追加 */
static inline double sgn(double x) {
    if (fabs(x) < EPS) return 0.0;
    return (x > 0.0) ? 1.0 : -1.0;
}


/* プログラムのスタート */
int main(void) {
    /* ファイルの宣言（変数） */
    FILE *fhi, *fho;
    FILE *fzi, *fzo;
    FILE *fMi, *fMo;
    FILE *fui, *fuo;
    FILE *fci, *fco;
    FILE *fEi, *fEo;
    FILE *ftaui, *ftauo;
    FILE *fslopeo;

    /* ファイルの宣言（水路形状） */
    FILE *fflume;

    /* ファイルの宣言（その他） */
    FILE *fM200o;   /* i=200 M 時系列 */

    /* 離散化の添字 */
    int i, n;

    /* 配列の宣言（リープフロッグ）- 全て X_GRID + 1 に統一 */
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

    double flume[X_GRID + 1] = { 0 };

    double slope_arr[X_GRID + 1] = {0};

    /* 浅水方程式のパラメータ */
    double d;
    double Mo, ho, co, uo;
    double stop;

    /* 構成則のパラメータ */
    double tauy, kkg, kkf, etao, fb;
    double cs;

    /* 物性・定数 */
    double ro = 1.0;
    double sigma = 2.65;
    double e = 0.775;
    double tanfais = 0.795436;
    double cstar = 0.70;
    double g = 980.0;
    double beta = 1.25;
    double kg = 0.0828;
    double kf;
    double kar = 0.40;

    double mu = 0.01;
    double T = 1.428571;

    cs = cstar / 2.0F;

    double limit_slope = atan(0.9*cstar*(sigma-ro)/(0.9*cstar*(sigma-ro)+ro)*tanfais);

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
    fhi = fopen(INPUT_DIR "hi.txt", "r");
    fzi = fopen(FLUME_FILE, "r");
    fMi = fopen(INPUT_DIR "Mi.txt", "r");
    fui = fopen(INPUT_DIR "ui.txt", "r");
    fci = fopen(INPUT_DIR "ci.txt", "r");
    fEi = fopen(INPUT_DIR "Ei.txt", "r");
    ftaui = fopen(INPUT_DIR "taui.txt", "r");
    fflume = fopen(FIXED_BED_FILE, "r");

    if (fhi == NULL || fzi == NULL || fMi == NULL || fui == NULL || 
        fci == NULL || fEi == NULL || ftaui == NULL || fflume == NULL) {
        printf("入力ファイルエラー\n");
        printf("  共通入力: %s\n", INPUT_DIR);
        printf("  水路形状: %s\n", FLUME_FILE);
        printf("  固定床:   %s\n", FIXED_BED_FILE);
        exit(1);
    }

    /* 出力ディレクトリの自動作成 */
    system("mkdir -p " OUTPUT_DIR);

    /* ファイルオープン（出力） */
    fho    = fopen(OUTPUT_DIR OUTPUT_PREFIX "h.txt", "w");
    fzo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "z.txt", "w");
    fMo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "M.txt", "w");
    fuo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "u.txt", "w");
    fco    = fopen(OUTPUT_DIR OUTPUT_PREFIX "c.txt", "w");
    fEo    = fopen(OUTPUT_DIR OUTPUT_PREFIX "E.txt", "w");
    ftauo  = fopen(OUTPUT_DIR OUTPUT_PREFIX "tau.txt", "w");
    fslopeo = fopen(OUTPUT_DIR OUTPUT_PREFIX "slope.txt", "w");
    
    fM200o  = fopen(OUTPUT_DIR OUTPUT_PREFIX "M_i200.txt", "w");

    if (fho == NULL || fzo == NULL || fMo == NULL || fuo == NULL ||
        fco == NULL || fEo == NULL || ftauo == NULL || fslopeo == NULL || fM200o == NULL) {
        printf("出力ファイルエラー\n");
        printf("  出力先: %s\n", OUTPUT_DIR);
        printf("  フォルダが存在するか確認してください\n");
        exit(1);
    }   
    fprintf(fM200o, "t[s] Q_i200[m^3/s]\n");

    printf("========== 江頭侵食モデル ==========\n");
    printf("  粒径:     %.2f cm\n", GRAIN_SIZE);
    printf("  水路形状: %s\n", FLUME_NAME);
    printf("  流量:     Mo=%.1f cm^2/s\n", M_IN);
    printf("  固定床:   %s\n", FIXED_BED_FILE);
    printf("  出力先:   %s\n", OUTPUT_DIR);
    printf("====================================\n");

    /* 計算開始前にCFL条件の初期情報を表示 */
    printf("\n===== CFL条件チェック設定 =====\n");
    printf("目標クーラン数: %.2f\n", CFL_TARGET);
    printf("dx = %.2f [cm], dt = %.4f [s]\n", dx, dt);
    printf("beta = %.2f, g = %.1f [cm/s^2]\n", beta, g);
    printf("特性速度: beta*u ± sqrt(g*h + beta*(beta-1)*u^2)\n");
    printf("================================\n\n");


    //=======================================================================

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
    z1[X_GRID] = BOTTOM;

    h2[0] = ho;
    h2[X_GRID] = 0;
    M2[0] = Mo;
    M2[X_GRID] = 0;
    u2[0] = uo;
    u2[X_GRID] = 0;
    c2[0] = co;
    c2[X_GRID] = 0;
    z2[0] = TOP;
    z2[X_GRID] = BOTTOM;

    E[0] = 0;
    E[X_GRID] = 0;
    tau[0] = 0;
    tau[X_GRID] = 0;

    /* ファイルから内部点の初期条件を読み込む */
    for (i = 1; i < X_GRID; i++) {
        fscanf(fhi, "%lf", &h1[i]); fscanf(fzi, "%lf", &z1[i]);
        fscanf(fMi, "%lf", &M1[i]); fscanf(fui, "%lf", &u1[i]);
        fscanf(fci, "%lf", &c1[i]); fscanf(fEi, "%lf", &E[i]);
        fscanf(ftaui, "%lf", &tau[i]); fscanf(fflume, "%lf", &flume[i]);
    }

    //======================================================================= 

    /* メインループ */
    for (n = 0; n < T_GRID; n++) {

        /* 給水条件 */
        if (n < INFLOW_TIME/(2*dt)) {
            h1[0] = ho; h2[0] = ho; M1[0] = Mo; M2[0] = Mo;
            u1[0] = uo; u2[0] = uo; c1[0] = co; c2[0] = co; ct[0] = co;
        } else {
            h1[0] = 0; h2[0] = 0; M1[0] = 0; M2[0] = 0;
            u1[0] = 0; u2[0] = 0; c1[0] = 0; c2[0] = 0; ct[0] = 0;
        }

        /* ===== CFL条件チェック ===== */
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
                    printf("[WARNING] CFL violation at n=%d: CFL=%.4f > %.2f (max_speed=%.2f cm/s)\n", n, CFL_current, CFL_TARGET, max_char_speed);
                    printf("         推奨dt <= %.6f [s]\n", CFL_TARGET * dx / max_char_speed);
                }
            }
        }

        /* Step.A1 連続式 */
        for (i = 1; i < X_GRID; i++) {
            /* 風上差分 */
            double M_diff;
            if(u1[i] >= 0){
                M_diff = M1[i] - M1[i-1];
            } else {
                M_diff = M1[i+1] - M1[i];
            }          

            h2[i] = h1[i] + dt * (E[i] - M_diff / dx);

            //if (h2[i] < h_min) { h2[i] = 0; c2[i] = 0; z2[i] = z1[i]; }
            if (h2[i] < h_min) { h2[i] = 0; c2[i] = 0; }

            else {
                if (c1[i] < c_min) ct[i] = 0;
                else {
                    ct[i] = c1[i];
                }

                /* 風上差分 */
                double cM_diff;
                if(u1[i] >= 0){
                    cM_diff = ct[i] * M1[i] - ct[i-1] * M1[i-1];
                } else {
                    cM_diff = ct[i+1] * M1[i+1] - ct[i] * M1[i];
                }       
                c2[i] = (c1[i] * h1[i] + dt * (E[i] * cstar - cM_diff / dx)) / h2[i];

                //if (z1[i] - E[i] * dt < flume[i]) { z2[i] = flume[i];}
                //else {
                   //z2[i] = z1[i] - E[i] * dt;
                //}
            }
        }

        /* 境界値の更新 */
        //z2[X_GRID] = 2*z2[X_GRID-1] - z2[X_GRID-2];
        h2[X_GRID] = h2[X_GRID-1];
        c2[X_GRID] = c2[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];


        /* Step.B1+C1 構成則+運動方程式 */
        for (i = 1; i < X_GRID; i++) {

            double slope = atan((z1[i] - z1[i + 1]) / dx);
           
            double theta_0 = atan(c2[i] * (sigma - ro) / (c2[i] * (sigma - ro) + ro) * tanfais);
            if (theta_0 < 0) theta_0=0;
            if (theta_0 >= limit_slope) theta_0 = limit_slope;

            if (h2[i] < stop) { tau[i] = 0; M2[i] = 0; u2[i] = 0; E[i] = 0; z2[i]=z1[i];}
            else {
                if(c2[i] < c_min) {
                    kf = 0.025; etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h2[i])*log(1 + h2[i] / etao) - 1), 2); tauy = 0;
                } else {
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c2[i], 0.333);
                    kf = 0.08;
                    kkf = kf * pow(1 - c2[i], 1.666) / pow(c2[i], 0.666);
                    fb = 6.25*(kkg + kkf)*pow(d / h2[i], 2);
                    tauy = pow(c2[i] / cstar, 0.2)*(sigma - ro)*c2[i] * g*h2[i] * cos(slope)*tanfais;
                }

                tau[i] = tauy*sgn(u1[i]) + ro * fb*u1[i]*fabs(u1[i]);

                double ro_m = fmax(ro + (sigma - ro)*c2[i], ro);

                /* 風上差分 */
                double uM_diff;
                if(u1[i] >= 0){
                    uM_diff = u1[i] * M1[i] - u1[i-1] * M1[i-1];
                } else {
                    uM_diff = u1[i+1] * M1[i+1] - u1[i] * M1[i];
                }          

                M2[i] = M1[i] - dt * (beta * uM_diff / dx + u1[i]*(cstar*(sigma-ro)+ro)*E[i] + g * h2[i] * (h2[i+1] + z1[i+1] - h2[i] - z1[i]) / dx + tau[i] / ro_m);

                u2[i] = M2[i] / h2[i];

                E[i] = fabs(u2[i]) * tan(slope - theta_0);

                /* 清水流の掘削防止 */
                if (c2[i] < c_min && E[i] < 0) E[i] = 0;

                /* Eの下限を設定(堆積しすぎ防止) */
                double E_min1 = -h2[i] * c2[i] / dt / cstar; // 濃度が非負
                double E_min2 = -h2[i] / dt;                 // 水深が非負
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                /* 固定床の掘りすぎ防止 */
                if (E[i]>0){
                    double dist = z1[i] - flume[i];
                    if(dist<0){E[i]=0;}
                    else if(dist < E[i]*dt){ E[i]=dist/dt;}
                    else{E[i]=E[i];}
                }

                slope_arr[i] = slope;
                z2[i] = z1[i] - E[i] * dt;
                
            }
        }

        /* 境界値の更新 */
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];

        /* Step.A2 連続式 */
        for (i = 1; i < X_GRID; i++) {
            /* 風上差分 */
            double M_diff;
            if(u2[i] >= 0){
                M_diff = M2[i] - M2[i-1];
            } else {
                M_diff = M2[i+1] - M2[i];
            }          

            h1[i] = h2[i] + dt * (E[i] - M_diff / dx);

            //if (h1[i] < h_min) { h1[i] = 0; c1[i] = 0; z1[i] = z2[i]; }
            if (h1[i] < h_min) { h1[i] = 0; c1[i] = 0; }

            else {
                if (c2[i] < c_min) ct[i] = 0;
                else {
                    ct[i] = c2[i];
                }

                /* 風上差分 */
                double cM_diff;
                if(u2[i] >= 0){
                    cM_diff = ct[i] * M2[i] - ct[i-1] * M2[i-1];
                } else {
                    cM_diff = ct[i+1] * M2[i+1] - ct[i] * M2[i];
                }    

                c1[i] = (c2[i] * h2[i] + dt * (E[i] * cstar - cM_diff / dx)) / h1[i];

                //if (z2[i] - E[i] * dt < flume[i]) { z1[i] = flume[i];}
                //else {
                    //z1[i] = z2[i] - E[i] * dt;
                //}
            }
            if(n % Record == Record - 1){fprintf(fho, "%f ", h1[i]);}
            if(n % Record == Record - 1){fprintf(fco, "%f ", c1[i]);}
            //if(n % Record == Record - 1){fprintf(fzo, "%f ", z1[i]);}
            if(n == T_GRID - 1){fprintf(fzo, "%f ", z1[i]);}
        }

        /* 境界値の更新 */
        //z1[X_GRID] = 2*z1[X_GRID-1] - z1[X_GRID-2];
        h1[X_GRID] = h1[X_GRID-1];
        c1[X_GRID] = c1[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];


        /* Step.B2+C2 */
        for (i = 1; i < X_GRID; i++) {
            double slope = atan((z2[i] - z2[i + 1]) / dx);
            double theta_0 = atan(c1[i] * (sigma - ro) / (c1[i] * (sigma - ro) + ro) * tanfais);
            if (theta_0 < 0) theta_0=0;
            if (theta_0 >= limit_slope) theta_0 = limit_slope;

            if (h1[i] < stop) { tau[i] = 0; M1[i] = 0; u1[i] = 0; E[i] = 0; z1[i]=z2[i];}
            else {
                if(c1[i] < c_min) {
                    kf = 0.025; etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h1[i])*log(1 + h1[i] / etao) - 1), 2); tauy = 0;
                } else {
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c1[i], 0.333);
                    kf = 0.08;
                    kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
                    fb = 6.25*(kkg + kkf)*pow(d / h1[i], 2);
                    tauy = pow(c1[i] / cstar, 0.2)*(sigma - ro)*c1[i] * g*h1[i] * cos(slope)*tanfais;
                }

                tau[i] = tauy*sgn(u2[i]) + ro * fb*u2[i]*fabs(u2[i]);

                double ro_m = fmax(ro + (sigma - ro)*c1[i], ro);

                /* 風上差分 */
                double uM_diff;
                if(u2[i] >= 0){
                    uM_diff = u2[i] * M2[i] - u2[i-1] * M2[i-1];
                } else {
                    uM_diff = u2[i+1] * M2[i+1] - u2[i] * M2[i];
                }          

                M1[i] = M2[i] - dt * (beta * uM_diff / dx + u2[i] * (cstar * (sigma - ro) + ro) * E[i] + g * h1[i] * (h1[i+1] + z2[i+1] - h1[i] - z2[i]) / dx + tau[i] / ro_m);

                u1[i] = M1[i] / h1[i];

                E[i] = fabs(u1[i]) * tan(slope - theta_0);

                /* 清水流の掘削防止 */
                if (c1[i] < c_min && E[i] < 0) E[i] = 0;

                /* Eの下限を設定(堆積しすぎ防止) */
                double E_min1 = -h1[i] * c1[i] / dt / cstar; // 濃度が非負
                double E_min2 = -h1[i] / dt;                 // 水深が非負
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                /* 固定床の掘りすぎ防止 */
                if (E[i]>0){
                    double dist = z2[i] - flume[i];
                    if(dist<0){E[i]=0;}
                    else if(dist < E[i]*dt){ E[i]=dist/dt;}
                    else{E[i]=E[i];}
                }   

                slope_arr[i] = slope;
                z1[i] = z2[i] - E[i] * dt;
            }

            if(n % Record == Record - 1){fprintf(ftauo, "%f ", tau[i]);}
            if(n % Record == Record - 1){fprintf(fMo,   "%f ", M1[i]);}
            if(n % Record == Record - 1){fprintf(fuo,   "%f ", u1[i]);}
            if(n % Record == Record - 1){fprintf(fEo,   "%f ", E[i]);}
            if(n % Record == Record - 1){fprintf(fslopeo, "%f ", slope_arr[i]);}
        }
        /* 境界値の更新 */
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];

        /* 2ループ目の終了 */
        if(n % Record == Record - 1){fprintf(fho, "\n");}
        if(n % Record == Record - 1){fprintf(fco, "\n");}
        //if(n % Record == Record - 1){fprintf(fzo, "\n");}
        if(n == T_GRID - 1){fprintf(fzo, "\n");}
        if(n % Record == Record - 1){fprintf(ftauo, "\n");}
        if(n % Record == Record - 1){fprintf(fMo, "\n");}
        if(n % Record == Record - 1){fprintf(fuo, "\n");}
        if(n % Record == Record - 1){fprintf(fEo, "\n");}
        if(n % Record == Record - 1){fprintf(fslopeo, "\n");}
        if(n % Record == Record - 1){fprintf(fM200o, "%.4f %.6f\n", 2.0*(n+1)*dt, M1[X_GRID] * 1e-4);}
    }

    //======================================================================= 

    /* ===== CFL条件の最終サマリー ===== */
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

    fclose(fhi); fclose(fho); fclose(fzi); fclose(fzo);
    fclose(fMi); fclose(fMo); fclose(fui); fclose(fuo);
    fclose(fci); fclose(fco); fclose(fEi); fclose(fEo);
    fclose(ftaui); fclose(ftauo); fclose(fslopeo); fclose(fflume); fclose(fM200o);

    printf("----- 計算完了（江頭侵食_%s） -----\n", FLUME_NAME);
    return 0;
}