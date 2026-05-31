/* 田澤実験-江頭式*/

/* 1. コメントとヘッダファイル */
/* 単位系は(cm,g,sec で) */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)


/* ========= ケース設定（ここだけ変えれば別ケースに対応できる） ========= */
#define CASE_ID          6      /* ケース番号 */


/* ========= ケースごとのパラメータ定義 ========= */
#define CASE_ID_STR STR(CASE_ID)
#if   CASE_ID == 1
    #define GRAIN_SIZE       0.67    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.67    /* 停止水深 stop [cm] */
    #define M_IN             116.6   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.14    /* 給水濃度 co [-] */
    #define H_IN             1.96    /* 給水水深 ho [cm] */
    #define INFLOW_TIME      30      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        38      /* 実験形状の枚数 */
    #define EXP_INTERVAL     1.0     /* 実験形状の間隔 */

#elif CASE_ID == 2
    #define GRAIN_SIZE       0.67    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.67    /* 停止水深 stop [cm] */
    #define M_IN             349.7   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.14    /* 給水濃度 co [-] */
    #define H_IN             3.14    /* 給水水深 ho [cm] */
    #define INFLOW_TIME      10      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        40      /* 実験形状の枚数 */
    #define EXP_INTERVAL     0.5     /* 実験形状の間隔 */

#elif CASE_ID == 3
    #define GRAIN_SIZE       0.28    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.28    /* 停止水深 stop [cm] */
    #define M_IN             126.5   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.21    /* 給水濃度 co [-] */
    #define H_IN             1.62    /* 給水水深 ho [cm] */
    #define INFLOW_TIME      30      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        37      /* 実験形状の枚数 */
    #define EXP_INTERVAL     1.0     /* 実験形状の間隔 */

#elif CASE_ID == 4
    #define GRAIN_SIZE       0.28    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.28    /* 停止水深 stop [cm] */
    #define M_IN             352.1   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.15    /* 給水濃度 co [-] */
    #define H_IN             2.32     /* 給水水深 ho [cm] */
    #define INFLOW_TIME      10      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        30      /* 実験形状の枚数 */
    #define EXP_INTERVAL     0.5     /* 実験形状の間隔 */

#elif CASE_ID == 5
    #define GRAIN_SIZE       0.14    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.14    /* 停止水深 stop [cm] */
    #define M_IN             116.3   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.14    /* 給水濃度 co [-] */
    #define H_IN             1.12    /* 給水水深 ho [cm] */
    #define INFLOW_TIME      30      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        36      /* 実験形状の枚数 */
    #define EXP_INTERVAL     1.0     /* 実験形状の間隔 */

#elif CASE_ID == 6
    #define GRAIN_SIZE       0.14    /* 粒径 d [cm] */
    #define STOP_DEPTH       0.14    /* 停止水深 stop [cm] */
    #define M_IN             369.9   /* 給水流量 Mo [cm^2/s] */
    #define C_IN             0.19  /* 給水濃度 co [-] |(369.9-300)/369.9 */
    #define H_IN             1.87    /* 給水水深 ho [cm] */
    #define INFLOW_TIME      10      /* 給水している時間ステップ数 n<この間だけ給水 */
    #define EXP_SNAPS        31      /* 実験形状の枚数 */
    #define EXP_INTERVAL     0.5     /* 実験形状の間隔 */

#else
    #error "Unsupported CASE_ID (1〜6のどれかを指定してください)"
#endif
/* ========================================================================= */


/* 計算格子 */
#define dx  5.0                               /* 空間間隔 */
#define dt  0.005F                             /* 時間間隔 */

/* スナップ1枚あたりの秒数（EXP_INTERVAL  / (2dt)）、0.5Fは浮動小数を回避 */
#define EXP_TIMEperSNAPS ((int)(EXP_INTERVAL/(2*dt)+0.5F))

#define X_GRID 193                           /* x方向グリッド数+1 (+1は境界条件) */
#define T_GRID (EXP_SNAPS*EXP_TIMEperSNAPS)  /* t方向グリッド数/2 */


/*  数値発散防止 */
const double EPS = 1e-10;
const double h_min = 0.01;  /* 乾湿判定の水深閾値 */
const double c_min = 0.01;  /* 清水流判定の濃度閾値 *//*  数値発散防止 */

#define Record ((int)(0.1/dt + 0.5))

/* CFL条件チェック用の定数 */
#define CFL_TARGET 0.85      /* 目標クーラン数 */
#define CFL_CHECK_INTERVAL 100  /* CFLチェックの間隔（ステップ数） */

/* 符号関数sgnの追加 */
static inline double sgn(double x) {
    if (fabs(x) < EPS) return 0.0;
    return (x > 0.0) ? 1.0 : -1.0;
}

/* 時系列データを収集する地点 */
#define POINT_1 126
#define POINT_2 134
#define POINT_3 138


/* 3. main関数とファイルポインタ */
/* プログラムのスタート */
int main(void) {
    /* ファイルの宣言（変数） *//* i:input o:output の意味 */
    FILE *fhi, *fho;
    FILE *fzi, *fzo;
    FILE *fMi, *fMo;
    FILE *fui, *fuo;
    FILE *fci, *fco;
    FILE *fEi, *fEo;
    FILE *ftaui, *ftauo;
    /* ファイルの宣言（θ用）追加 */
    FILE *fslopeo;    /* 勾配角度 slope [rad] */
    /* ファイルの宣言（水路形状） */
    FILE *fflume;
    /* ファイルの宣言（実験の堆積形状） */
    FILE *fzexp;

    /* ファイルの宣言（感度分析） */
    FILE *fRepo;
    FILE *fH_Pexo;
    FILE *fL_Pexo;
    FILE *fPP_Pexo;
    FILE *fPmaxo;
    FILE *fPmino;

    /* ファイルの宣言（堆積厚さD） */
    FILE *fDo;

    /* 感度分析用統合CSVファイル */
    FILE *fcsv;

    /* 時系列データ用CSVファイル（i=126, 134, 138） */
    FILE *fts_126;
    FILE *fts_134;
    FILE *fts_134_txt;
    FILE *fts_138;


    /* 離散化の添字 */
    int i; /* i:x方向の番号 Qoとかはi=0に当てる */
    int n; /* n:t方向の番号 初期値はn=0に当てる */

    /* 配列の宣言（リープフロッグ） */
    /* z:河床高,h:水深,M:流量,u:速度,c:濃度,E:侵食速度,tau:河床面せん断力 */
    double zansa[X_GRID + 1] = { 0 };
    double zexp[EXP_SNAPS*192+1] = { 0 };

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
    /* slope配列を追加 */
    double slope_arr[X_GRID + 1] = { 0 };

    double ct[X_GRID + 1] = { 0 }; 
    double phi[X_GRID + 1] = { 0 }; 
    double phis[X_GRID + 1] = { 0 }; 
    double phiw[X_GRID + 1] = { 0 }; 
    double phii[X_GRID + 1] = { 0 }; 

    double flume[X_GRID + 1] = { 0 }; 

    /* 感度分析の配列を追加 */
    double Rep_arr[X_GRID + 1] = { 0 }; 
    double H_Pex_arr[X_GRID + 1] = { 0 }; 
    double L_Pex_arr[X_GRID + 1] = { 0 }; 
    double PP_Pex_arr[X_GRID + 1] = { 0 }; 
    double Pmax_arr[X_GRID + 1] = { 0 }; 
    double Pmin_arr[X_GRID + 1] = { 0 }; 

    /* 堆積厚さD = z - flume */
    double D_arr[X_GRID + 1] = { 0 }; 



    /* 浅水方程式のパラメータ */
    double d;                  /* 粒径 */
    double Mo, ho, co, uo;     /* 給水の流量,水深,濃度,流速(oはオー) */
    double stop;               /* 土石流停止の基準 */

    /* 構成則のパラメータ（仮の初期値を設定） */
    double tauy, kkg, kkf, gg, etao, fb, r;
    double cs;

    /* 物性・定数 */
    double ro = 1.0;             /* 水の密度 */
    double sigma = 2.65;        /* 砂の密度 */
    double e = 0.85;            /* 反発係数 */
    double tanfais = 0.7954;    /* tanΦs */
    double cstar = 0.58;        /* 堆砂濃度 */
    double g = 980.0;            /* 重力加速度(cm/s/s) */
    double beta = 1.25;          /* 運動量補正係数 */
    double kg = 0.0828;         /* kg(構成則中の定数) */
    double kf;                  /* kf(構成則中の定数) */
    double kar = 0.40;          /* kar(構成則中カルマン定数) */
    double mu = 0.01;           /* μ（粘性係数） */

    cs = cstar / 2.0F;          /* 土石流と土砂流の境界 */


    /* 供給条件の設定 */
    d = GRAIN_SIZE;  
    stop = STOP_DEPTH;
    Mo = M_IN;
    co = C_IN;

    ho = H_IN;
    uo = Mo / ho;

    /* CFL条件チェック用の変数 */
    double max_char_speed = 0.0;      /* 最大特性速度 */
    double CFL_current = 0.0;         /* 現在のクーラン数 */
    double CFL_max_overall = 0.0;     /* 計算全体での最大クーラン数 */
    int CFL_violation_count = 0;      /* CFL違反回数 */
    double char_speed_plus, char_speed_minus, char_speed_local;
    double discriminant;

    /* 残差計算 */
    double b = 0.00;
    double RMSE = 0.00;


    /* hの初期値hi.txtと結果用のho.txtを読み込み */
    fhi = fopen("txtファイル/hi.txt", "r");
    fho = fopen("江頭/E_case" CASE_ID_STR "/江頭ho" CASE_ID_STR ".txt", "w");
    if (fhi == NULL || fho == NULL) {
        printf("error h\n");
        exit(1);
    }
    /* zの初期値zi.txt 結果用zo.txt */
    fzi = fopen("txtファイル/zi.txt", "r");
    fzo = fopen("江頭/E_case" CASE_ID_STR "/江頭zo" CASE_ID_STR ".txt", "w");
    if (fzi == NULL || fzo == NULL) {
        printf("error z\n");
      exit(1);
    }
    /* Eの初期値Ei.txt 結果用Eo.txt */
    fEi = fopen("txtファイル/Ei.txt", "r");
    fEo = fopen("江頭/E_case" CASE_ID_STR "/江頭Eo" CASE_ID_STR ".txt", "w");
    if (fEi == NULL || fEo == NULL) {
        printf("error E\n");
        exit(1);
    }
    /* tauの初期値taui.txt 結果用tauo.txt */
    ftaui = fopen("txtファイル/taui.txt", "r");
    ftauo = fopen("江頭/E_case" CASE_ID_STR "/江頭tau" CASE_ID_STR ".txt", "w");
    if (ftaui == NULL || ftauo == NULL) {
        printf("error tau\n");
        exit(1);
    }
    /* Mの初期値Mi.txt 結果用Mo.txt */
    fMi = fopen("txtファイル/Mi.txt", "r");
  fMo = fopen("江頭/E_case" CASE_ID_STR "/江頭Mo" CASE_ID_STR ".txt", "w");
  if (fMi == NULL || fMo == NULL) {
      printf("error M\n");
      exit(1);
  }
    /* uの初期値ui.txt 結果用uo.txt */
    fui = fopen("txtファイル/ui.txt", "r");
    fuo = fopen("江頭/E_case" CASE_ID_STR "/江頭uo" CASE_ID_STR ".txt", "w");
    if (fui == NULL || fuo == NULL) {
        printf("error u\n");
        exit(1);
    }
    /* cの初期値ci.txt 結果用co.txt */
    fci = fopen("txtファイル/ci.txt", "r");
    fco = fopen("江頭/E_case" CASE_ID_STR "/江頭co" CASE_ID_STR ".txt", "w");
    if (fci == NULL || fco == NULL) {
        printf("error c\n");
        exit(1);
    }

    /* slope(θ)の結果用ファイル */
    fslopeo = fopen("江頭/E_case" CASE_ID_STR "/江頭slopeo" CASE_ID_STR ".txt", "w");
    if (fslopeo == NULL) {
        printf("error slope\n");
        exit(1);
    }

    /* 水路形状flume.txt */
    fflume = fopen("txtファイル/flume.txt", "r");
    if (fflume == NULL) {
        printf("error flume\n");
        exit(1);
    }

    /* 実験結果zexp.txt */
    fzexp = fopen("expファイル/zexp" CASE_ID_STR ".txt", "r");
    if (fzexp == NULL) {
        printf("error zexp\n");
        exit(1);
    }


    /* 感度分析用ファイル */
    fRepo = fopen("江頭/E_case" CASE_ID_STR "/江頭Repo" CASE_ID_STR ".txt", "w");
    if (fRepo == NULL) {printf("error Rep\n");exit(1);}
    fH_Pexo = fopen("江頭/E_case" CASE_ID_STR "/江頭H_Pexo" CASE_ID_STR ".txt", "w");
    if (fH_Pexo == NULL) {printf("error H_Pex\n");exit(1);}
    fL_Pexo = fopen("江頭/E_case" CASE_ID_STR "/江頭L_Pexo" CASE_ID_STR ".txt", "w");
    if (fL_Pexo == NULL) {printf("error L_Pex\n");exit(1);}
    fPP_Pexo = fopen("江頭/E_case" CASE_ID_STR "/江頭PP_Pexo" CASE_ID_STR ".txt", "w");
    if (fPP_Pexo == NULL) {printf("error PP_Pex\n");exit(1);}
    fPmaxo = fopen("江頭/E_case" CASE_ID_STR "/江頭Pmaxo" CASE_ID_STR ".txt", "w");
    if (fPmaxo == NULL) {printf("error Pmax\n");exit(1);}
    fPmino = fopen("江頭/E_case" CASE_ID_STR "/江頭Pmino" CASE_ID_STR ".txt", "w");
    if (fPmino == NULL) {printf("error Pmin\n");exit(1);}

    /* 堆積厚さD = z - flume */
    fDo = fopen("江頭/E_case" CASE_ID_STR "/江頭Do" CASE_ID_STR ".txt", "w");
    if (fDo == NULL) {
        printf("error D\n");
        exit(1);
    }

    /* 感度分析用統合CSVファイル */
    fcsv = fopen("江頭/E_case" CASE_ID_STR "/sensitivity_data" CASE_ID_STR ".csv", "w");
    if (fcsv == NULL) {printf("error csv\n");exit(1);}
    /* CSVヘッダー出力 */
    fprintf(fcsv, "i,n,time,x,c,h,Rep,H_Pex,Pmax,Pmin,D\n");

    /* 時系列データ用CSVファイル（i=126, 134, 138） */
    fts_126 = fopen("江頭/E_case" CASE_ID_STR "/timeseries_i126_" CASE_ID_STR ".csv", "w");
    if (fts_126 == NULL) {printf("error timeseries i=126\n");exit(1);}
    fprintf(fts_126, "n,time,c,h,Rep,H_Pex,Pmax,Pmin,D\n");

    fts_134 = fopen("江頭/E_case" CASE_ID_STR "/timeseries_i134_" CASE_ID_STR ".csv", "w");
    if (fts_134 == NULL) {printf("error timeseries i=134\n");exit(1);}
    fprintf(fts_134, "n,time,c,h,Rep,H_Pex,Pmax,Pmin,D\n");

    fts_134_txt = fopen("江頭/E_case" CASE_ID_STR "/timeseries_i134_" CASE_ID_STR ".txt", "w");
    if (fts_134_txt == NULL) {printf("error timeseries txt i=134\n");exit(1);}

    fts_138 = fopen("江頭/E_case" CASE_ID_STR "/timeseries_i138_" CASE_ID_STR ".csv", "w");
    if (fts_138 == NULL) {printf("error timeseries i=138\n");exit(1);}
    fprintf(fts_138, "n,time,c,h,Rep,H_Pex,Pmax,Pmin,D\n");




    //=======================================================================

    /* ファイルに存在しない端点（i=0, i=X_GRID）の初期条件 */
    h1[0] = ho;
    h1[X_GRID] = 0; 
    M1[0] = Mo; 
    M1[X_GRID] = 0;
    u1[0] = uo;
    u1[X_GRID] = 0;
    c1[0] = co;
    c1[X_GRID] = 0;
    z1[0] = 186.2417709; 
    z1[X_GRID] = -0.17;

    h2[0] = ho;
    h2[X_GRID] = 0; 
    M2[X_GRID] = 0;
    u2[0] = uo;
    u2[X_GRID] = 0;
    c2[0] = co;
    c2[X_GRID] = 0;
    z2[0] = 186.2417709; 
    z2[X_GRID] = -0.17;

    E[0] = 0;
    E[X_GRID] = 0;
    tau[0] = 0;
    tau[X_GRID] = 0;
    flume[0] = 186.2417709;
    flume[X_GRID] = -0.17;


    /* ファイルから "内部点" の初期条件を読み込む */
    for (i = 1; i < X_GRID; i++) { /* 初期値を変数に入れる */
        fscanf(fhi, "%lf", &h1[i]);
        fscanf(fzi, "%lf", &z1[i]);
        fscanf(fMi, "%lf", &M1[i]);
        fscanf(fui, "%lf", &u1[i]);
        fscanf(fci, "%lf", &c1[i]);
        fscanf(fEi, "%lf", &E[i]);
        fscanf(ftaui, "%lf", &tau[i]);
        fscanf(fflume, "%lf", &flume[i]);
    }

    for (i = 1; i < EXP_SNAPS*192+1; i++) {
        fscanf(fzexp, "%lf", &zexp[i]);
    }

    /* 流速係数の初期計算（構成則パラメータ） */
    for (i = 0; i < X_GRID; i++) {
        if (h1[i] > h_min && c1[i] > c_min) {

            kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c1[i], 0.333);
            kf = 0.16;
            kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
            etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
            gg = (sigma / ro - 1)*c1[i] + 1;

            if(c1[i]/ cs < 1){
            phiw[i] = pow(1 - c1[i] / cs, 0.5)*(etao / h1[i] + 1 - c1[i] / cs)*log((etao / h1[i] + 1 - c1[i] / cs)*h1[i] / etao)*d / h1[i] / kar - pow(1 - c1[i] / cs, 1.5)*d / h1[i] / kar; 

            r = 1 + cs / c1[i] * ((1 - pow(c1[i] / cstar, 0.2))*gg - 1);
            if (fabs(r) <= 0.00001) { 
                phis[i] = pow(c1[i] / cs, 2)*pow(1 - c1[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5); 
                phii[i] = c1[i] / cs * pow(1 - c1[i] / cs, 1.5) / pow(kkf + kkg, 0.5); 
            }
            else {
                phis[i] = (pow(1 - c1[i] / cs, 2.5) - pow(1 + (r - 1)*c1[i] / cs, 1.5)*(1 - (3 * r / 2 + 1)*c1[i] / cs)) * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2); 
                phii[i] = (pow(1 + (r - 1)*c1[i] / cs, 1.5) - pow(1 - c1[i] / cs, 1.5))*(1 - c1[i] / cs) * 2 / 3 / r / pow(kkg + kkf, 0.5);
            }

            phi[i] = phiw[i] + phis[i] + phii[i];
            }
        }
        else {
            phiw[i] = 0;
            phis[i] = 0;
            phii[i] = 0;
            phi[i] = 0;
        }
    }

    /* 計算開始前にCFL条件の初期情報を表示 */
    printf("===== CFL条件チェック設定 =====\n");
    printf("目標クーラン数: %.2f\n", CFL_TARGET);
    printf("dx = %.2f [cm], dt = %.4f [s]\n", dx, dt);
    printf("beta = %.2f, g = %.1f [cm/s^2]\n", beta, g);
    printf("特性速度: beta*u ± sqrt(g*h + beta*(beta-1)*u^2)\n");
    printf("================================\n\n");

    //======================================================================= 

    /* メインループ */
    for (n = 0; n < T_GRID; n++) {

        /* 給水条件 */
        if (n < INFLOW_TIME/(2*dt)) {
            h1[0] = ho;
            h2[0] = ho;
            M1[0] = Mo;
            M2[0] = Mo;
            u1[0] = uo;
            u2[0] = uo;
            c1[0] = co;
            c2[0] = co;
            ct[0] = co;
        }
        else {
            h1[0] = 0;
            h2[0] = 0;
            M1[0] = 0;
            M2[0] = 0;
            u1[0] = 0;
            u2[0] = 0;
            c1[0] = 0;
            c2[0] = 0;
            ct[0] = 0;
        }

        /* ===== CFL条件チェック（各ステップ開始時） ===== */
        max_char_speed = 0.0;
        for (i = 1; i < X_GRID; i++) {
            if (h1[i] > h_min) {
                /* 特性速度: beta*u ± sqrt(g*h + beta*(beta-1)*u^2) */
                discriminant = g * h1[i] + beta * (beta - 1.0) * u1[i] * u1[i];
                if (discriminant > 0) {
                    double sqrt_disc = sqrt(discriminant);
                    char_speed_plus = beta * u1[i] + sqrt_disc;
                    char_speed_minus = beta * u1[i] - sqrt_disc;
                    /* 絶対値の最大を取る */
                    char_speed_local = fmax(fabs(char_speed_plus), fabs(char_speed_minus));
                    if (char_speed_local > max_char_speed) {
                        max_char_speed = char_speed_local;
                    }
                }
            }
        }

        /* クーラン数の計算 */
        if (max_char_speed > EPS) {
            CFL_current = max_char_speed * dt / dx;
            if (CFL_current > CFL_max_overall) {
                CFL_max_overall = CFL_current;
            }

            /* CFL違反チェック */
            if (CFL_current > CFL_TARGET) {
                CFL_violation_count++;
                /* 最初の違反時と定期的に警告を出力 */
                if (CFL_violation_count == 1 || n % CFL_CHECK_INTERVAL == 0) {
                    printf("[WARNING] CFL violation at n=%d: CFL=%.4f > %.2f (max_speed=%.2f cm/s)\n",
                           n, CFL_current, CFL_TARGET, max_char_speed);
                    /* 推奨dtの計算 */
                    double dt_recommended = CFL_TARGET * dx / max_char_speed;
                    printf("         推奨dt <= %.6f [s]\n", dt_recommended);
                }
            }
        }
        /* ===== CFL条件チェック終了 ===== */

        /* Step.A1 連続式 */
        for (i = 1; i < X_GRID; i++) {

            /* 風上差分 */
            double M_diff;
            if(u1[i] >= 0){M_diff = M1[i] - M1[i-1];} 
            else {M_diff = M1[i+1] - M1[i];}          

            h2[i] = h1[i] + dt * (E[i] - M_diff / dx);

            //ドライセル判定
            if (h2[i] <= h_min) {
                h2[i] = 0;
                c2[i] = 0; 
                z2[i] = z1[i];
            }

            //ドライセルでない
            else{
                if (c1[i] <= c_min) {ct[i] = 0;}
                else if (c1[i]/cs >= 1 ||phis[i] <= EPS||phi[i]<=EPS){ct[i] = c1[i];}
                else{ct[i] = cs * phis[i] / phi[i];}

                /* 風上差分 */
                double cM_diff;
                if(u1[i] >= 0){cM_diff = ct[i] * M1[i] - ct[i-1] * M1[i-1];}
                else {cM_diff = ct[i+1] * M1[i+1] - ct[i] * M1[i];}

                c2[i] = (c1[i] * h1[i] + dt * (E[i] * cstar - cM_diff / dx)) / h2[i];

                if (z1[i] - E[i] * dt < flume[i]) { z2[i] = flume[i];}
                else {
                    z2[i] = z1[i] - E[i] * dt;
                }
            }

            /* A1 の出力（zは出力しない） */
            //fprintf(fho, "%f ", h2[i]);
            //fprintf(fco, "%f ", c2[i]);
            //fprintf(fzo, "%f ", z2[i]); 
        }
        /* 境界値の更新 */
        z2[X_GRID] = 2*z2[X_GRID-1] - z2[X_GRID-2];
        h2[X_GRID] = h2[X_GRID-1];
        c2[X_GRID] = c2[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];

        /* Step.B1+C1 構成則+運動方程式 */
        for (i = 1; i < X_GRID; i++) {

            phi[i] = phis[i] =phiw[i] =phii[i] = 0.0;

            Rep_arr[i]=0.0; 
            H_Pex_arr[i]=0.0;
            L_Pex_arr[i]=0.0;
            PP_Pex_arr[i]=0.0;
            Pmax_arr[i]=0.0;
            Pmin_arr[i]=0.0;
            D_arr[i]=0.0;


            //int sgn_u = (u1[i] > 0.0) ? 1 : (u1[i] < 0.0) ? -1 : 0;
            double slope = atan((z2[i] - z2[i+1]) / dx);
            //if(slope<0) slope = atan((z2[i]+ h2[i] - z2[i+1]-h2[i+1]) / dx);

            double theta_0 =atan(c2[i] * (sigma - ro) / (c2[i] * (sigma - ro) + ro)*tanfais);

            /* 堆積厚さDの計算 */
            D_arr[i] = z2[i] - flume[i];

            if (h2[i] < stop ) {tau[i] = 0; M2[i] = 0; u2[i] = 0; E[i] = 0; slope_arr[i] = slope;}

            else{
                if(c2[i] <= c_min){
                    kf = 0.025;
                    etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h2[i])*log(1 + h2[i] / etao) - 1), 2);
                    tauy=0;
                }
                else if (c2[i] / cs >= 1) {
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c2[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c2[i], 1.666) / pow(c2[i], 0.666);
                    fb = 6.25*(kkg + kkf)*pow(d / h2[i], 2);
                    tauy = pow(c2[i] / cstar, 0.2)*(sigma - ro)*c2[i] * g*h2[i] * cos(slope)*tanfais;
                }
                else{
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c2[i], 0.333);
                    kf = 0.16;
                    kkf = kf * pow(1 - c2[i], 1.666) / pow(c2[i], 0.666);
                    etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
                    gg = (sigma / ro - 1)*c2[i] + 1;

                    phiw[i] = pow(1 - c2[i] / cs, 0.5)*(etao / h2[i] + 1 - c2[i] / cs)*log((etao / h2[i] + 1 - c2[i] / cs)*h2[i] / etao)*d / h2[i] / kar - pow(1 - c2[i] / cs, 1.5)*d / h2[i] / kar;

                    r = 1 + cs / c2[i] * ((1 - pow(c2[i] / cstar, 0.2))*gg - 1);
                    if (fabs(r) <= 0.00001) {
                        phis[i] = pow(c2[i] / cs, 2)*pow(1 - c2[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5);
                        phii[i] = c2[i] / cs * pow(1 - c2[i] / cs, 1.5) / pow(kkf + kkg, 0.5);
                    }
                    else {
                        phis[i] = (pow(1 - c2[i] / cs, 2.5) - pow(1 + (r - 1)*c2[i] / cs, 1.5)*(1 - (3 * r / 2 + 1)*c2[i] / cs)) * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2);
                        phii[i] = (pow(1 + (r - 1)*c2[i] / cs, 1.5) - pow(1 - c2[i] / cs, 1.5))*(1 - c2[i] / cs) * 2 / 3 / r / pow(kkg + kkf, 0.5);
                    }
                    phi[i] = phiw[i] + phis[i] + phii[i];
                    if(phi[i]<=EPS){fb = 6.25*(kkg + kkf)*pow(d / h2[i], 2);}
                    else{fb = (1 - pow(c2[i] / cstar, 0.2))*gg*pow(d / h2[i] / phi[i], 2);}
                    tauy = pow(c2[i] / cstar, 0.2)*(sigma - ro)*c2[i] * g*h2[i] * cos(slope)*tanfais; 
                }

                tau[i] = tauy*sgn(u1[i]) + ro * fb*u1[i]*fabs(u1[i]);

                double ro_m = fmax(ro + (sigma - ro)*c2[i], ro);

                /* 風上差分 */
                double uM_diff;
                if(u1[i] >= 0){uM_diff = u1[i] * M1[i] - u1[i-1] * M1[i-1];}
                else {uM_diff = u1[i+1] * M1[i+1] - u1[i] * M1[i];}  

                //摩擦なしの値
                //double M_adv=M1[i] - dt * (beta * uM_diff / dx + g *h2[i] *(h2[i+1]+z2[i+1]- h2[i]-z2[i]) / dx);
                //M2[i] = (M_adv -dt*tauy*sgn(u1[i]))/(1+dt*ro*fb*fabs(u1[i])/ro_m/h2[i]);

                M2[i] = M1[i] - dt * (beta * uM_diff / dx + g * h2[i] * (h2[i+1] + z2[i+1] - h2[i] - z2[i]) / dx + tau[i] / ro_m);

                u2[i] = M2[i] / h2[i];


                E[i] = fabs(u2[i]) * tan(slope - theta_0);    

                /* 清水流の掘削防止 */
                if (c2[i] <= c_min && E[i] < 0) E[i] = 0;

                /* Eの下限を設定(堆積しすぎ防止) */
                double E_min1 = -h2[i] * c2[i] / dt / cstar; // 濃度が非負
                double E_min2 = -h2[i] / dt;                 // 水深が非負
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                /* 固定床の掘りすぎ防止 */
                if (E[i]>0){
                    double dist = z2[i] - flume[i];
                    if(dist<0){E[i]=0;}
                    else if(dist < E[i]*dt){ E[i]=dist/dt;}
                    else{E[i]=E[i];}
                }

                Rep_arr[i]=(1-c2[i])*ro*d*fabs(E[i])/mu;
                H_Pex_arr[i]=(18*mu* c2[i] /d/ d* h2[i]*(-E[i])+ 27*ro* c2[i] /8 /d *h2[i]*(-E[i])*fabs(E[i])) /ro/g/cos(slope);
                double L_Pex;
                if(0.2 < c2[i]){L_Pex= 150*mu*c2[i]*c2[i] /d/ d /(1-c2[i]) * h2[i]*(-E[i])+7*ro* c2[i]/4/d *h2[i]*(-E[i])*fabs(E[i]);}
                else if (1000<= Rep_arr[i]){L_Pex= 0.33*ro* c2[i]/d/pow(1-c2[i], 1.65)*Rep_arr[i]*h2[i]*(-E[i])*fabs(E[i]);}
                else{L_Pex= 18*ro*c2[i]/d/pow(1-c2[i], 1.65)*(1+0.15*pow(Rep_arr[i], 0.687)) *h2[i]*(-E[i])*fabs(E[i]);}
                L_Pex_arr[i]=L_Pex /ro/g/cos(slope);
                PP_Pex_arr[i]=(150*mu*c2[i]*c2[i] /d/ d/(1-c2[i])* h2[i]*(-E[i]))/ro/g/cos(slope);
                Pmax_arr[i]=(c2[i]*(sigma - ro)*h2[i]);
                Pmin_arr[i]= -h2[i];

                /* slopeを配列に保存 */
                slope_arr[i] = slope;
            }         
            /* B1 の出力（値だけ、改行なし）*/
            //fprintf(ftauo, "%f ", tau[i]);
            //fprintf(fMo, "%f ", M2[i]);
            //fprintf(fuo, "%f ", u2[i]);
            //fprintf(fEo, "%f ", E[i]);
            //fprintf(fslopeo, "%f ", slope_arr[i]);
        }

        /* 1ループ目の終了（zは出力しない） */
        //fprintf(fho, "\n");
        //fprintf(fco, "\n"); 
        //fprintf(fzo, "\n"); 
        //fprintf(ftauo, "\n");
        //fprintf(fMo, "\n"); 
        //fprintf(fuo, "\n"); 
        //fprintf(fEo, "\n");
        //fprintf(fslopeo, "\n");

        /* 境界値の更新 */
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];

        /* Step.A2 連続式 */
        for (i = 1; i < X_GRID; i++) {
            /* 風上差分 */
            double M_diff;
            if(u2[i] >= 0){M_diff = M2[i] - M2[i-1];}
            else { M_diff = M2[i+1] - M2[i];}          

            h1[i] = h2[i] + dt * (E[i] - M_diff / dx);

            //ドライセル判定
            if (h1[i] <= h_min) {
                h1[i] = 0;
                c1[i] = 0;
                z1[i] = z2[i];
            }
            else{
                if (c2[i] <= c_min) {ct[i] = 0;}
                else if (c2[i]/cs >= 1 ||phis[i]<=EPS||phi[i] <= EPS){ct[i] = c2[i];}
                else{ct[i] = cs * phis[i] / phi[i];}

                /* 風上差分 */
                double cM_diff;
                if(u2[i] >= 0){cM_diff = ct[i] * M2[i] - ct[i-1] * M2[i-1];}
                else {cM_diff = ct[i+1] * M2[i+1] - ct[i] * M2[i];}

                c1[i] = (c2[i] * h2[i] + dt * (E[i] * cstar - cM_diff / dx)) / h1[i];

                if (z2[i] - E[i] * dt < flume[i]) {z1[i] = flume[i];}
                else {z1[i] = z2[i] - E[i] * dt;}
            }

            /* A2 の出力（zはexpに合わせて出力） */
            if(n % Record == Record - 1 ){fprintf(fho, "%f ", h1[i]);}
            if(n % Record == Record - 1 ){fprintf(fco, "%f ", c1[i]);}

            if(n % EXP_TIMEperSNAPS == EXP_TIMEperSNAPS-1 ){
                fprintf(fzo, "%f ", z1[i]);
            }
        }
        /* 境界値の更新 */
        z1[X_GRID] = 2*z1[X_GRID-1] - z1[X_GRID-2];
        h1[X_GRID] = h1[X_GRID-1];
        c1[X_GRID] = c1[X_GRID-1];
        ct[X_GRID] = ct[X_GRID-1];
        M2[X_GRID] = M2[X_GRID-1];
        u2[X_GRID] = u2[X_GRID-1];

        /* Step.B2+C2 連続式 */
        for (i = 1; i < X_GRID; i++) {
            phi[i] = phis[i] =phiw[i] =phii[i] = 0.0;
            //int sgn_u = (u2[i] > 0.0) ? 1 : (u2[i] < 0.0) ? -1 : 0;

            Rep_arr[i]=0.0; 
            H_Pex_arr[i]=0.0;
            L_Pex_arr[i]=0.0;
            PP_Pex_arr[i]=0.0;
            Pmax_arr[i]=0.0;
            Pmin_arr[i]=0.0;
            D_arr[i]=0.0;

            double slope = atan((z1[i] - z1[i+1]) / dx);
            //if(slope<0) slope = atan((z1[i]+ h1[i] - z1[i+1]-h1[i+1]) / dx);

            double theta_0 =atan(c1[i] * (sigma - ro) / (c1[i] * (sigma - ro) + ro)*tanfais);

            /* 堆積厚さDの計算 */
            D_arr[i] = z1[i] - flume[i];

            if (h1[i] < stop ) {tau[i] = 0; M1[i] = 0; u1[i] = 0; E[i] = 0; slope_arr[i] = slope;}

            else{
                if (c1[i] <= c_min) {
                    kf = 0.025;
                    etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333);
                    fb = pow(kar / ((1 + etao / h1[i])*log(1 + h1[i] / etao) - 1), 2);
                    tauy=0;
                }
                else if (c1[i] / cs >= 1) {
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c1[i], 0.333); 
                    kf = 0.16; 
                    kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
                    fb = 6.25*(kkg + kkf)*pow(d / h1[i], 2);
                    tauy = pow(c1[i] / cstar, 0.2)*(sigma - ro)*c1[i] * g*h1[i] * cos(slope)*tanfais; 
                }
                else{
                    kkg = kg * sigma / ro * (1 - pow(e, 2))*pow(c1[i], 0.333); 
                    kf = 0.16; 
                    kkf = kf * pow(1 - c1[i], 1.666) / pow(c1[i], 0.666);
                    etao = pow(kf, 0.5)*d*pow((1 - cs) / cs, 0.333); 
                    gg = (sigma / ro - 1)*c1[i] + 1;

                    phiw[i] = pow(1 - c1[i] / cs, 0.5)*(etao / h1[i] + 1 - c1[i] / cs)*log((etao / h1[i] + 1 - c1[i] / cs)*h1[i] / etao)*d / h1[i] / kar - pow(1 - c1[i] / cs, 1.5)*d / h1[i] / kar; 

                    r = 1 + cs / c1[i] * ((1 - pow(c1[i] / cstar, 0.2))*gg - 1);
                    if (fabs(r) <= 0.00001) { 
                        phis[i] = pow(c1[i] / cs, 2)*pow(1 - c1[i] / cs, 0.5) / 2 / pow(kkf + kkg, 0.5); 
                        phii[i] = c1[i] / cs * pow(1 - c1[i] / cs, 1.5) / pow(kkf + kkg, 0.5); 
                    }
                    else {
                        phis[i] = (pow(1 - c1[i] / cs, 2.5) - pow(1 + (r - 1)*c1[i] / cs, 1.5)*(1 - (3 * r / 2 + 1)*c1[i] / cs)) * 4 / 15 / pow(kkg + kkf, 0.5) / pow(r, 2); 
                        phii[i] = (pow(1 + (r - 1)*c1[i] / cs, 1.5) - pow(1 - c1[i] / cs, 1.5))*(1 - c1[i] / cs) * 2 / 3 / r / pow(kkg + kkf, 0.5);
                    }
                    phi[i] = phiw[i] + phis[i] + phii[i];
                    if(phi[i]<=EPS){ fb = 6.25*(kkg + kkf)*pow(d / h1[i], 2); }
                    else{ fb = (1 - pow(c1[i] / cstar, 0.2))*gg*pow(d / h1[i] / phi[i], 2);}
                    tauy = pow(c1[i] / cstar, 0.2)*(sigma - ro)*c1[i] * g*h1[i] * cos(slope)*tanfais;

                }
                tau[i] = tauy*sgn(u2[i]) + ro * fb*u2[i]*fabs(u2[i]);

                double ro_m = fmax(ro + (sigma - ro)*c1[i], ro);

                /* 風上差分 */
                double uM_diff;
                if(u2[i] >= 0){uM_diff = u2[i] * M2[i] - u2[i-1] * M2[i-1];} 
                else {uM_diff = u2[i+1] * M2[i+1] - u2[i] * M2[i];} 

                //摩擦なしの値
                //double M_adv=M2[i] - dt * (beta * uM_diff / dx + g *h1[i] *(h1[i+1]+z1[i+1]- h1[i]-z1[i]) / dx);
                //M1[i] = (M_adv -dt*tauy*sgn(u2[i]))/(1+dt*ro*fb*fabs(u2[i])/ro_m/h1[i]);

                M1[i] = M2[i] - dt * (beta * uM_diff / dx + g * h1[i] * (h1[i+1] + z1[i+1] - h1[i] - z1[i]) / dx + tau[i] / ro_m);

                u1[i] = M1[i] / h1[i];

                E[i] = fabs(u1[i]) * tan(slope - theta_0);

                /* 清水流の掘削防止 */
                if (c1[i] <= c_min && E[i] < 0) E[i] = 0;
                
                /* Eの下限を設定(堆積しすぎ防止) */
                double E_min1 = -h1[i] * c1[i] / dt / cstar; // 濃度が非負
                double E_min2 = -h1[i] / dt;                 // 水深が非負
                E[i] = fmax(E[i], fmax(E_min1, E_min2));

                /* 固定床の掘りすぎ防止 */
                if (E[i]>0){
                    double dist = z1[i] - flume[i];
                    if(dist<0){E[i]=0;}
                    else if(dist < E[i]*dt){ E[i]=dist/dt;}
                    else{E[i]=E[i];}
                }

                Rep_arr[i]=(1-c1[i])*ro*d*fabs(E[i])/mu;
                H_Pex_arr[i]=(18*mu* c1[i] /d/ d* h1[i]*(-E[i])+ 27*ro* c1[i] /8 /d *h1[i]*(-E[i])*fabs(E[i])) /ro/g/cos(slope);
                double L_Pex;
                if(0.2 < c1[i]){L_Pex= 150*mu*c1[i]*c1[i] /d/ d /(1-c1[i]) * h1[i]*(-E[i])+7*ro* c1[i]/4/d *h1[i]*(-E[i])*fabs(E[i]);}
                else if (1000<= Rep_arr[i]){L_Pex= 0.33*ro* c1[i]/d/pow(1-c1[i], 1.65)*Rep_arr[i]*h1[i]*(-E[i])*fabs(E[i]);}
                else{L_Pex= 18*ro*c1[i]/d/pow(1-c1[i], 1.65)*(1+0.15*pow(Rep_arr[i], 0.687)) *h1[i]*(-E[i])*fabs(E[i]);}
                L_Pex_arr[i]=L_Pex /ro/g/cos(slope);
                PP_Pex_arr[i]=(150*mu*c1[i]*c1[i] /d/ d/(1-c1[i])* h1[i]*(-E[i]))/ro/g/cos(slope);
                Pmax_arr[i]=(c1[i]*(sigma - ro)*h1[i]);
                Pmin_arr[i]= -h1[i];
                /* slopeを配列に保存 */
                slope_arr[i] = slope;
            }

            /* B2 の出力：値だけ */
            if(n % Record == Record - 1 ){fprintf(ftauo, "%f ", tau[i]);}
            if(n % Record == Record - 1 ){fprintf(fMo,   "%f ", M1[i]);}
            if(n % Record == Record - 1 ){fprintf(fuo,   "%f ", u1[i]);}
            if(n % Record == Record - 1 ){fprintf(fEo,   "%f ", E[i]);}
            if(n % Record == Record - 1){fprintf(fslopeo, "%f ", slope_arr[i]);}

            if(n % EXP_TIMEperSNAPS == EXP_TIMEperSNAPS-1 ){
                zansa[i] = pow(z1[i]-zexp[i+((n-(EXP_TIMEperSNAPS-1))/EXP_TIMEperSNAPS)*192], 2);
                b += zansa[i];
            }

            if(n % Record == Record - 1 ){fprintf(fRepo,   "%f ", Rep_arr[i]);}
            if(n % Record == Record - 1 ){fprintf(fH_Pexo,   "%f ", H_Pex_arr[i]);}
            if(n % Record == Record - 1 ){fprintf(fL_Pexo,   "%f ", L_Pex_arr[i]);}
            if(n % Record == Record - 1 ){fprintf(fPP_Pexo,   "%f ", PP_Pex_arr[i]);}
            if(n % Record == Record - 1 ){fprintf(fPmaxo,   "%f ", Pmax_arr[i]);}
            if(n % Record == Record - 1 ){fprintf(fPmino,   "%f ", Pmin_arr[i]);}

            /* 堆積厚さD */
            if(n % Record == Record - 1 ){fprintf(fDo, "%f ", D_arr[i]);}

            /* 感度分析用統合CSV出力 */
            if(n % Record == Record - 1 ){
                double time_sec = (n + 1) * 2.0 * dt;  /* 実時間[s] */
                double x_pos = i * dx;         /* 位置[cm] */
                fprintf(fcsv, "%d,%d,%.4f,%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                        i, n, time_sec, x_pos,
                        c1[i], h1[i],
                        Rep_arr[i], H_Pex_arr[i],
                        Pmax_arr[i], Pmin_arr[i], D_arr[i]);
            }

            /* 時系列データ出力（i=126, 134, 138） */
            if(n % Record == Record - 1 ){
                double time_sec = (n + 1) * 2.0 * dt;  /* 実時間[s] */
                if(i == POINT_1){
                    fprintf(fts_126, "%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            n, time_sec,
                            c1[i], h1[i],
                            Rep_arr[i], H_Pex_arr[i],
                            Pmax_arr[i], Pmin_arr[i], D_arr[i]);
                }
                if(i == POINT_2){
                    fprintf(fts_134, "%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            n, time_sec,
                            c1[i], h1[i],
                            Rep_arr[i], H_Pex_arr[i],
                            Pmax_arr[i], Pmin_arr[i], D_arr[i]);
                    fprintf(fts_134_txt, "%d %.4f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
                            n, time_sec,
                            c1[i], h1[i],
                            Rep_arr[i], H_Pex_arr[i],
                            Pmax_arr[i], Pmin_arr[i], D_arr[i]);
                }
                if(i == POINT_3){
                    fprintf(fts_138, "%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                            n, time_sec,
                            c1[i], h1[i],
                            Rep_arr[i], H_Pex_arr[i],
                            Pmax_arr[i], Pmin_arr[i], D_arr[i]);
                }
            }

        }
        /* 境界値の更新 */
        M1[X_GRID] = M1[X_GRID-1];
        u1[X_GRID] = u1[X_GRID-1];

        /* 2ループ目の終了 */
        if(n % Record == Record - 1 ){fprintf(fho, "\n");}
        if(n % Record == Record - 1 ){fprintf(fco, "\n");}

        if(n % EXP_TIMEperSNAPS == EXP_TIMEperSNAPS-1 ){
            fprintf(fzo, "\n");
        }

        if(n % Record == Record - 1 ){fprintf(ftauo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fMo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fuo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fEo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fslopeo, "\n");}

        if(n % Record == Record - 1 ){fprintf(fRepo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fH_Pexo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fL_Pexo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fPP_Pexo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fPmaxo, "\n");}
        if(n % Record == Record - 1 ){fprintf(fPmino, "\n");}

        /* 堆積厚さD */
        if(n % Record == Record - 1 ){fprintf(fDo, "\n");}


    }
    //======================================================================= 

    /* --- メインループ終了後に zexp との比較 --- */
    RMSE = pow(b / ((double)(X_GRID - 1)*(double)EXP_SNAPS), 0.5);  // = sqrt(b / 192.0)
    printf("江頭RMSE (final z vs experiment) = %f\n", RMSE);

    /* --- 最終時刻で z[i] - flume[i] > d となる範囲の幅を計算 --- */
    double D =0.0; 
    double i_D =0.0;
    for (i = 1; i < X_GRID; i++) {
        if (z1[i] - flume[i] >= d) {
            D += z1[i] - flume[i];
            i_D += i*(z1[i] - flume[i]);
        }
    }
    double X_G = i_D/D*dx;
    printf("総堆積量:D=%f, 堆積重心:X_G=%f\n", D, X_G);

    /* --- zexp（実験データ）でも同様に計算 --- */
    double D_exp = 0.0;
    double i_D_exp = 0.0;
    int exp_final_offset = (EXP_SNAPS - 1) * 192;  /* 最終スナップのオフセット */
    for (i = 1; i < X_GRID; i++) {
        double z_exp_final = zexp[i + exp_final_offset];
        if (z_exp_final - flume[i] >= d) {
            D_exp += z_exp_final - flume[i];
            i_D_exp += i * (z_exp_final - flume[i]);
        }
    }
    double X_G_exp = i_D_exp / D_exp * dx;
    printf("実験 総堆積量:D_exp=%f, 堆積重心:X_G_exp=%f\n", D_exp, X_G_exp);

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

    /* --- 最終時刻で z1[i] - flume[i] >= d となる範囲の開始/到達地点（計算結果） --- */
    int i_start = -1;
    int i_end   = -1;

    for (i = 1; i < X_GRID; i++) {
        if (z1[i] - flume[i] >= d) {
            if (i_start == -1) i_start = i;  // 最初に見つかったi
            i_end = i;                        // 見つかるたびに更新 → 最終的に最大i
        }
    }

    if (i_start == -1) {
        printf("計算 堆積なし (z1 - flume >= d を満たす点がありません)\n");
    } else {
        double X_start = i_start * dx;
        double X_end   = i_end   * dx;
        printf("計算 堆積開始地点:X_start=%f, 堆積到達地点:X_end=%f\n", X_start, X_end);
        /* ついでに幅も欲しければ */
        /* printf("計算 堆積範囲幅:W=%f\n", (i_end - i_start) * dx); */
    }

    /* --- zexp（実験データ）でも同様に開始/到達地点 --- */
    /* --- zexp（実験データ）開始/到達地点（exp_final_offset は既に定義済みのものを使う） --- */
    int i_start_exp = -1;
    int i_end_exp   = -1;

    for (i = 1; i < X_GRID; i++) {
        double z_exp_final = zexp[i + exp_final_offset];
        if (z_exp_final - flume[i] >= d) {
            if (i_start_exp == -1) i_start_exp = i;
            i_end_exp = i;
        }
    }

    if (i_start_exp == -1) {
        printf("実験 堆積なし (zexp - flume >= d)\n");
    } else {
        printf("実験 堆積開始地点:X_start_exp=%f, 堆積到達地点:X_end_exp=%f\n",
               i_start_exp * dx, i_end_exp * dx);
    }





    /* ファイルクローズ */
    fclose(fhi);
    fclose(fho);
    fclose(fzi);
    fclose(fzo);
    fclose(fMi);
    fclose(fMo);
    fclose(fui);
    fclose(fuo);
    fclose(fci);
    fclose(fco);
    fclose(fEi);
    fclose(fEo);
    fclose(ftaui);
    fclose(ftauo);
    fclose(fslopeo);
    fclose(fflume);
    fclose(fzexp);

    fclose(fRepo);
    fclose(fH_Pexo);
    fclose(fL_Pexo);
    fclose(fPP_Pexo);
    fclose(fPmaxo);
    fclose(fPmino);
    fclose(fcsv);      /* 感度分析用CSVのクローズ */

    /* 堆積厚さD */
    fclose(fDo);

    /* 時系列データ用CSVのクローズ */
    fclose(fts_126);
    fclose(fts_134);
    fclose(fts_134_txt);
    fclose(fts_138);


    printf("----- 計算完了しました（江頭case" CASE_ID_STR "） -----\n");
    printf("----- 感度分析用CSV: 江頭/E_case" CASE_ID_STR "/sensitivity_data" CASE_ID_STR ".csv -----\n");
    printf("----- 時系列データCSV: 江頭/E_case" CASE_ID_STR "/timeseries_i126_" CASE_ID_STR ".csv -----\n");
    printf("----- 時系列データCSV: 江頭/E_case" CASE_ID_STR "/timeseries_i134_" CASE_ID_STR ".csv -----\n");
    printf("----- 時系列データCSV: 江頭/E_case" CASE_ID_STR "/timeseries_i138_" CASE_ID_STR ".csv -----\n");
    printf("----- 堆積厚さD: 江頭/E_case" CASE_ID_STR "/江頭Do" CASE_ID_STR ".txt -----\n");
    return 0;
}