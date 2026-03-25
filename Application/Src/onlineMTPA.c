

#include "onlineMTPA.h"
#include <math.h>
#include <string.h>

/* ---------------- 全局可调变量（上位机可写） ---------------- */
volatile float g_ld = 0.030f;
volatile float g_ldd = 0.0f;
volatile float g_lq = 0.0165f;
volatile float g_lqq = 0.0f;
volatile float g_lc = 0.0f;

volatile float g_beta = 0.95f;
volatile float g_lambda = 1e-8f;
volatile float g_smoothFcHz = 5.0f;
volatile float g_outputFcHz10k = 5.0f;

volatile uint16_t g_binSamples = 40u; /* 2ms */
volatile uint16_t g_edgeSamples = 5u;
volatile uint16_t g_solveDecim = 1u; /* 每10ms解一次 */

volatile uint8_t g_fixLdLq = 0u;
volatile uint8_t g_identEnable = 1u;
volatile float g_identIsMin = 2.0f;
volatile float g_identWeMin = 5.0f; /* rad/s，低速冻结 */
volatile uint8_t g_freezeOnVsat = 1u;

volatile uint8_t g_freezeOnSteady = 1u;
volatile float g_steadyDidTh = 0.15f;
volatile float g_steadyDiqTh = 0.15f;
volatile uint8_t g_torqueGradFreezeSel = ONLINE_MTPA_TORQUE_GRAD_FREEZE_OFF;
volatile uint8_t Sel12 = 1u;
volatile uint32_t g_torqueGradFreezeSelDelayCnt = 5000u;
volatile float g_torqueGradTh = 0.02f;
volatile float g_torqueGradNormTh = 0.02f;
volatile float g_torqueGradFiltFcHz = 20.0f;
volatile float g_torqueGradRestartDeltaIs = 0.5f;
volatile uint8_t g_identFrozenByTorqueGrad = 0u;
volatile uint8_t g_identifying = 0u;
volatile float g_mtpaThetaRad = 0.0f;
volatile int8_t g_mtpaIqSign = 1;
volatile float g_torqueProxyNow = 0.0f;
volatile float g_torqueGradRawNow = 0.0f;
volatile float g_torqueGradFiltNow = 0.0f;
volatile float g_torqueGradNormNow = 0.0f;

volatile float g_gammaMin = 0.75f + 1e-3f;
volatile float g_gammaMax = 1.5707963f - 1e-3f;
volatile float g_gammaStepMaxDeg = 0.1f;
volatile float gamma_deg = 0.0f;
volatile uint8_t enable_45 = 0.0f;
/* ---------------- 内部状态 ---------------- */
typedef struct
{
  /* bin计数 */
  uint16_t cnt;       /* 当前bin内累计的10k点数 */
  uint32_t bins_acc;  /* 累计bin数 */
  uint16_t solve_cnt; /* 解算分频计数 */

  /* 10k积分累加 */
  float sum_bd, sum_bq; /* ∫(vd-Rs0*id)dt , ∫(vq-Rs0*iq)dt */
  float sum_wphid[5];   /* ∫ we*phi_d dt （5维） */
  float sum_wphiq[5];   /* ∫ we*phi_q dt （5维） */

  /* 起止均值 */
  float start_sum_id, start_sum_iq;
  uint16_t start_cnt;
  float end_sum_id, end_sum_iq;
  uint16_t end_cnt;

  /* EWLS 法方程 */
  float S[5][5];
  float t[5];

  /* 当前解 & 平滑解 */
  float xhat[5];
  float xhat_f[5];
  float xhat_out[5];

  /* 平滑系数（随 binSamples 变化） */
  float alpha;
  float alpha_out_10k;
} ident_state_t;

static ident_state_t st;
static float s_param_init[5];
static float s_gamma_prev = 0.0f;
static uint8_t s_gamma_prev_valid = 0u;
static uint8_t s_mtpa_theta_valid = 0u;
static float s_torqueGradAlpha10k = 1.0f;
static float s_torqueGradFiltState = 0.0f;
static uint8_t s_torqueGradFiltValid = 0u;
static float s_frozenIs = 0.0f;
static uint32_t s_identWarmCnt = 0u;
static uint8_t s_waitRecountByDeltaIs = 0u;
/* ---------------- 工具函数（小计算量） ---------------- */
static inline float clampf(float x, float a, float b)
{
  if (x < a) return a;
  if (x > b) return b;
  return x;
}

static inline float sqrf(float x)
{
  return x * x;
}

static void abort_current_bin(void)
{
  st.cnt = 0u;
  st.sum_bd = 0.0f;
  st.sum_bq = 0.0f;
  for (int i = 0; i < 5; i++)
  {
    st.sum_wphid[i] = 0.0f;
    st.sum_wphiq[i] = 0.0f;
  }
  st.start_sum_id = 0.0f;
  st.start_sum_iq = 0.0f;
  st.start_cnt = 0u;
  st.end_sum_id = 0.0f;
  st.end_sum_iq = 0.0f;
  st.end_cnt = 0u;
}

static void clear_ident_history(void)
{
  memset(st.S, 0, sizeof(st.S));
  memset(st.t, 0, sizeof(st.t));
  st.bins_acc = 0u;
  st.solve_cnt = 0u;
  abort_current_bin();
}

static void basis_5(float id, float iq, float phi_d[5], float phi_q[5])
{
  /* phi_d = [id, id^3, 0, 0, id*iq^2]
     phi_q = [0, 0, iq, iq^3, iq*id^2] */
  float id2 = id * id;
  float iq2 = iq * iq;

  phi_d[0] = id;
  phi_d[1] = id2 * id;
  phi_d[2] = 0.0f;
  phi_d[3] = 0.0f;
  phi_d[4] = id * iq2;

  phi_q[0] = 0.0f;
  phi_q[1] = 0.0f;
  phi_q[2] = iq;
  phi_q[3] = iq2 * iq;
  phi_q[4] = iq * id2;
}

static void update_alpha(void)
{
  /* alpha = 1 - exp(-2*pi*fc*dt_update) */
  float fc = g_smoothFcHz;
  uint16_t M = g_binSamples;
  if (M < 1u) M = 1u;
  if (M > ONLINE_MTPA_MAX_BIN_SAMPLES) M = ONLINE_MTPA_MAX_BIN_SAMPLES;

  uint16_t dec = g_solveDecim;
  if (dec < 1u) dec = 1u;
  float dt = (float)M * (float)dec * ONLINE_MTPA_TS_CTRL;
  if (fc <= 0.0f)
  {
    st.alpha = 1.0f;
    return;
  }
  float a = expf(-2.0f * 3.1415926f * fc * dt);
  st.alpha = 1.0f - a;
}

static void update_output_alpha_10k(void)
{
  float fc = g_outputFcHz10k;
  if (fc <= 0.0f)
  {
    st.alpha_out_10k = 1.0f;
    return;
  }
  float a = expf(-2.0f * 3.1415926f * fc * ONLINE_MTPA_TS_CTRL);
  st.alpha_out_10k = 1.0f - a;
}

static void update_torqueGrad_alpha_10k(void)
{
  float fc = g_torqueGradFiltFcHz;
  if (fc <= 0.0f)
  {
    s_torqueGradAlpha10k = 1.0f;
    return;
  }
  float a = expf(-2.0f * 3.1415926f * fc * ONLINE_MTPA_TS_CTRL);
  s_torqueGradAlpha10k = 1.0f - a;
}

static void update_torqueGrad_freeze_sel_by_count(uint8_t ident_active_now)
{
  if (!ident_active_now)
  {
    /* 辨识未激活时，清计数并关闭转矩偏导冻结 */
    s_identWarmCnt = 0u;
    g_torqueGradFreezeSel = ONLINE_MTPA_TORQUE_GRAD_FREEZE_OFF;
    return;
  }

  if (s_identWarmCnt < 0xFFFFFFFFu)
  {
    s_identWarmCnt++;
  }

  if (s_identWarmCnt > g_torqueGradFreezeSelDelayCnt)
  {
    uint8_t sel = Sel12;
    if ((sel != ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT) &&
        (sel != ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT_OVER_T))
    {
      sel = ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT;
    }
    g_torqueGradFreezeSel = sel;
  }
  else
  {
    g_torqueGradFreezeSel = ONLINE_MTPA_TORQUE_GRAD_FREEZE_OFF;
  }
}

static void publish_output_params(const float x[5])
{
  g_ld = x[0];
  g_ldd = x[1];
  g_lq = x[2];
  g_lqq = x[3];
  g_lc = x[4];
}

static void output_filter_step_10k(float Is2)
{
  float a = clampf(st.alpha_out_10k, 0.0f, 1.0f);
  float IsMin2 = sqrf(g_identIsMin);
  const float* xtarget = (Is2 < IsMin2) ? s_param_init : st.xhat_f;
  for (int i = 0; i < 5; i++)
  {
    st.xhat_out[i] += a * (xtarget[i] - st.xhat_out[i]);
  }
  publish_output_params(st.xhat_out);
}

/* Cholesky 求解 5x5 SPD: A x = b
   返回1成功，0失败（非正定/数值问题） */
static uint8_t solve_spd_chol5(float A[5][5], float b[5], float x[5])
{
  float L[5][5] = {0};

  for (int i = 0; i < 5; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      float sum = A[i][j];
      for (int k = 0; k < j; k++) sum -= L[i][k] * L[j][k];

      if (i == j)
      {
        if (sum <= 1e-12f) return 0;
        L[i][j] = sqrtf(sum);
      }
      else
      {
        L[i][j] = sum / L[j][j];
      }
    }
  }

  /* forward: L y = b */
  float y[5];
  for (int i = 0; i < 5; i++)
  {
    float sum = b[i];
    for (int k = 0; k < i; k++) sum -= L[i][k] * y[k];
    y[i] = sum / L[i][i];
  }

  /* backward: L^T x = y */
  for (int i = 4; i >= 0; i--)
  {
    float sum = y[i];
    for (int k = i + 1; k < 5; k++) sum -= L[k][i] * x[k];
    x[i] = sum / L[i][i];
  }

  return 1;
}

/* Cholesky 求解 3x3 SPD */
static uint8_t solve_spd_chol3(float A[3][3], float b[3], float x[3])
{
  float L[3][3] = {0};

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      float sum = A[i][j];
      for (int k = 0; k < j; k++) sum -= L[i][k] * L[j][k];

      if (i == j)
      {
        if (sum <= 1e-12f) return 0;
        L[i][j] = sqrtf(sum);
      }
      else
      {
        L[i][j] = sum / L[j][j];
      }
    }
  }

  float y[3];
  for (int i = 0; i < 3; i++)
  {
    float sum = b[i];
    for (int k = 0; k < i; k++) sum -= L[i][k] * y[k];
    y[i] = sum / L[i][i];
  }

  for (int i = 2; i >= 0; i--)
  {
    float sum = y[i];
    for (int k = i + 1; k < 3; k++) sum -= L[k][i] * x[k];
    x[i] = sum / L[i][i];
  }

  return 1;
}

/* ---------------- 公共API ---------------- */
void onlineMTPA_init(void)
{
  memset(&st, 0, sizeof(st));
  st.xhat[0] = g_ld;
  st.xhat[1] = g_ldd;
  st.xhat[2] = g_lq;
  st.xhat[3] = g_lqq;
  st.xhat[4] = g_lc;

  for (int i = 0; i < 5; i++) s_param_init[i] = st.xhat[i];
  for (int i = 0; i < 5; i++) st.xhat_f[i] = st.xhat[i];
  for (int i = 0; i < 5; i++) st.xhat_out[i] = st.xhat_f[i];

  update_alpha();
  update_output_alpha_10k();
  update_torqueGrad_alpha_10k();
  publish_output_params(st.xhat_out);
  s_gamma_prev = 0.0f;
  s_gamma_prev_valid = 0u;
  s_mtpa_theta_valid = 0u;
  s_torqueGradFiltState = 0.0f;
  s_torqueGradFiltValid = 0u;
  s_frozenIs = 0.0f;
  s_identWarmCnt = 0u;
  s_waitRecountByDeltaIs = 0u;
  g_identFrozenByTorqueGrad = 0u;
  g_identifying = 0u;
  g_mtpaThetaRad = 0.0f;
  g_mtpaIqSign = 1;
  g_torqueProxyNow = 0.0f;
  g_torqueGradRawNow = 0.0f;
  g_torqueGradFiltNow = 0.0f;
  g_torqueGradNormNow = 0.0f;
}

void onlineMTPA_reset_ident(void)
{
  memset(&st, 0, sizeof(st));
  st.xhat[0] = g_ld;
  st.xhat[1] = g_ldd;
  st.xhat[2] = g_lq;
  st.xhat[3] = g_lqq;
  st.xhat[4] = g_lc;

  for (int i = 0; i < 5; i++) s_param_init[i] = st.xhat[i];
  for (int i = 0; i < 5; i++) st.xhat_f[i] = st.xhat[i];
  for (int i = 0; i < 5; i++) st.xhat_out[i] = st.xhat_f[i];

  update_alpha();
  update_output_alpha_10k();
  update_torqueGrad_alpha_10k();
  publish_output_params(st.xhat_out);
  s_gamma_prev = 0.0f;
  s_gamma_prev_valid = 0u;
  s_mtpa_theta_valid = 0u;
  s_torqueGradFiltState = 0.0f;
  s_torqueGradFiltValid = 0u;
  s_frozenIs = 0.0f;
  s_identWarmCnt = 0u;
  s_waitRecountByDeltaIs = 0u;
  g_identFrozenByTorqueGrad = 0u;
  g_identifying = 0u;
  g_mtpaThetaRad = 0.0f;
  g_mtpaIqSign = 1;
  g_torqueProxyNow = 0.0f;
  g_torqueGradRawNow = 0.0f;
  g_torqueGradFiltNow = 0.0f;
  g_torqueGradNormNow = 0.0f;
}

onlineMTPA_status_t onlineMTPA_get_status(void)
{
  onlineMTPA_status_t s;
  s.bins_accumulated = st.bins_acc;
  s.valid = (st.bins_acc >= 10u) ? 1u : 0u;
  return s;
}

void onlineMTPA_flux_from_i(float id, float iq, float* psid, float* psiq)
{
  float ld = g_ld;
  float ldd = g_ldd;
  float lq = g_lq;
  float lqq = g_lqq;
  float lc = g_lc;

  float id2 = id * id;
  float iq2 = iq * iq;

  if (psid) *psid = ld * id + ldd * (id2 * id) + lc * id * iq2;
  if (psiq) *psiq = lq * iq + lqq * (iq2 * iq) + lc * iq * id2;
}

void onlineMTPA_incL_from_i(float id, float iq, float* Ldd, float* Lqq, float* Ldq)
{
  float ld = g_ld;
  float ldd = g_ldd;
  float lq = g_lq;
  float lqq = g_lqq;
  float lc = g_lc;

  float id2 = id * id;
  float iq2 = iq * iq;

  float Ldd_inc = ld + 3.0f * ldd * id2 + lc * iq2;
  float Lqq_inc = lq + 3.0f * lqq * iq2 + lc * id2;
  float Ldq_inc = 2.0f * lc * id * iq;

  if (Ldd) *Ldd = Ldd_inc;
  if (Lqq) *Lqq = Lqq_inc;
  if (Ldq) *Ldq = Ldq_inc;
}

float onlineMTPA_torque_proxy(float id, float iq)
{
  float psid, psiq;
  onlineMTPA_flux_from_i(id, iq, &psid, &psiq);
  return psid * iq - psiq * id;
}

void onlineMTPA_torque_grad_from_theta(float Is, float theta_rad, int8_t iq_sign, float* Tproxy,
                                       float* dTdTheta)
{
  float I = fabsf(Is);
  float c = COS(theta_rad);
  float s = SIN(theta_rad);
  float I2 = I * I;
  float c2 = c * c;
  float s2 = s * s;
  float sign = (iq_sign >= 0) ? 1.0f : -1.0f;

  float A = g_ld - g_lq;
  float B = g_ldd - g_lc;
  float C = g_lc - g_lqq;

  float g = A + B * I2 * c2 + C * I2 * s2;
  float T = sign * I2 * c * s * g;
  float dT = sign * I2 * ((c2 - s2) * g + 2.0f * I2 * c2 * s2 * (C - B));

  if (Tproxy) *Tproxy = T;
  if (dTdTheta) *dTdTheta = dT;
}

static void update_torque_grad_freeze(float Is, uint8_t base_gate_ok)
{
  uint8_t sel = g_torqueGradFreezeSel;

  if (sel == ONLINE_MTPA_TORQUE_GRAD_FREEZE_OFF)
  {
    g_identFrozenByTorqueGrad = 0u;
    /* 若正在等待“电流变化超阈值”后重启计数，则保留冻结电流参考 */
    if (!s_waitRecountByDeltaIs) s_frozenIs = 0.0f;
    s_torqueGradFiltState = 0.0f;
    s_torqueGradFiltValid = 0u;
    g_torqueGradRawNow = 0.0f;
    g_torqueGradFiltNow = 0.0f;
    g_torqueGradNormNow = 0.0f;
    g_torqueProxyNow = 0.0f;
    return;
  }

  if (g_identFrozenByTorqueGrad)
  {
    float dIs = (g_torqueGradRestartDeltaIs > 0.0f) ? g_torqueGradRestartDeltaIs : 0.0f;
    if (Is > (s_frozenIs + dIs))
    {
      g_identFrozenByTorqueGrad = 0u;
      s_frozenIs = 0.0f;
      s_torqueGradFiltState = 0.0f;
      s_torqueGradFiltValid = 0u;
      abort_current_bin();
    }
    return;
  }

  if ((!base_gate_ok) || (!s_mtpa_theta_valid) || (Is <= 1e-6f))
  {
    return;
  }

  float T = 0.0f;
  float dT = 0.0f;
  float dT_prev = s_torqueGradFiltState;
  uint8_t dT_prev_valid = s_torqueGradFiltValid;
  onlineMTPA_torque_grad_from_theta(Is, g_mtpaThetaRad, g_mtpaIqSign, &T, &dT);
  g_torqueProxyNow = T;
  g_torqueGradRawNow = dT;

  if (!s_torqueGradFiltValid)
  {
    s_torqueGradFiltState = dT;
    s_torqueGradFiltValid = 1u;
  }
  else
  {
    float a = clampf(s_torqueGradAlpha10k, 0.0f, 1.0f);
    s_torqueGradFiltState += a * (dT - s_torqueGradFiltState);
  }
  g_torqueGradFiltNow = s_torqueGradFiltState;
  g_torqueGradNormNow = fabsf(s_torqueGradFiltState) / fmaxf(fabsf(T), 1e-6f);

  float metric = 0.0f;
  float th = 0.0f;
  if (sel == ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT_OVER_T)
  {
    metric = g_torqueGradNormNow;
    th = g_torqueGradNormTh;
  }
  else
  {
    metric = fabsf(s_torqueGradFiltState);
    th = g_torqueGradTh;
  }

  if ((!dT_prev_valid) || (th <= 0.0f))
  {
    return;
  }

  /* 冻结只允许“负偏导变大并首次进入冻结区间”触发；
     正值进入阈值区或负值继续变小都不触发 */
  uint8_t freeze_hit = 0u;
  if (sel == ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT)
  {
    uint8_t in_zone_prev = (dT_prev >= -th) && (dT_prev <= th);
    uint8_t in_zone_now = (s_torqueGradFiltState >= -th) && (s_torqueGradFiltState <= th);
    if ((!in_zone_prev) && in_zone_now && (dT_prev < 0.0f) && (s_torqueGradFiltState > dT_prev))
    {
      freeze_hit = 1u;
    }
  }
  else
  {
    float prev_metric = fabsf(dT_prev) / fmaxf(fabsf(T), 1e-6f);
    uint8_t enter_zone = (prev_metric >= th) && (metric < th);
    if (enter_zone && (dT_prev < 0.0f) && (s_torqueGradFiltState > dT_prev))
    {
      freeze_hit = 1u;
    }
  }

  if (freeze_hit)
  {
    g_identFrozenByTorqueGrad = 1u;
    s_frozenIs = Is;
    s_waitRecountByDeltaIs = 1u;
    s_torqueGradFiltState = 0.0f;
    s_torqueGradFiltValid = 0u;
    clear_ident_history();
  }
}

/* 10kHz 辨识 step */
void onlineMTPA_step_10k(float vd, float vq, float id, float iq, float we, float Rs0, uint8_t flags)
{
  float Is2 = id * id + iq * iq;
  float Is = sqrtf(Is2);
  uint8_t base_gate_ok = 1u;

  output_filter_step_10k(Is2);
  update_torqueGrad_alpha_10k();
  g_identifying = 0u;

  if (Is < g_identIsMin) base_gate_ok = 0u;
  if (fabsf(we) < g_identWeMin) base_gate_ok = 0u;
  if (g_freezeOnVsat && (flags & ONLINE_MTPA_FLAG_VSAT)) base_gate_ok = 0u;
  if (flags & ONLINE_MTPA_FLAG_BAD_V) base_gate_ok = 0u;

  update_torque_grad_freeze(Is, base_gate_ok);

  /* 冻结后不立即重计数，需等待电流相对冻结点变化超过阈值 */
  if (s_waitRecountByDeltaIs)
  {
    float dIs = (g_torqueGradRestartDeltaIs > 0.0f) ? g_torqueGradRestartDeltaIs : 0.0f;
    if (fabsf(Is - s_frozenIs) > dIs)
    {
      s_waitRecountByDeltaIs = 0u;
      abort_current_bin();
    }
  }

  /* 解除冻结后，且满足辨识门控时开始计数；
     计满前保持 g_torqueGradFreezeSel=0，计满后置为 Sel12 */
  uint8_t ident_active_now =
      (g_identEnable && base_gate_ok && (!g_identFrozenByTorqueGrad) && (!s_waitRecountByDeltaIs))
          ? 1u
          : 0u;
  update_torqueGrad_freeze_sel_by_count(ident_active_now);

  if (!g_identEnable)
  {
    abort_current_bin();
    return;
  }
  if (s_waitRecountByDeltaIs)
  {
    abort_current_bin();
    return;
  }
  if (g_identFrozenByTorqueGrad)
  {
    abort_current_bin();
    return;
  }

  uint16_t M = g_binSamples;
  if (M < 1u) M = 1u;
  if (M > ONLINE_MTPA_MAX_BIN_SAMPLES) M = ONLINE_MTPA_MAX_BIN_SAMPLES;

  uint16_t edgeN = g_edgeSamples;
  if (edgeN < 1u) edgeN = 1u;
  if (edgeN > M) edgeN = M;

  /* 门控：低速/小电流/电压饱和冻结（直接不更新这一bin） */
  if (Is < g_identIsMin)
  {
    abort_current_bin();
    return;
  }
  if (fabsf(we) < g_identWeMin)
  {
    abort_current_bin();
    return;
  }
  if (g_freezeOnVsat && (flags & ONLINE_MTPA_FLAG_VSAT))
  {
    abort_current_bin();
    return;
  }
  if (flags & ONLINE_MTPA_FLAG_BAD_V)
  {
    abort_current_bin();
    return;
  }

  /* bin开始：清零累计 */
  g_identifying = 1u;

  if (st.cnt == 0u)
  {
    st.sum_bd = 0.0f;
    st.sum_bq = 0.0f;
    for (int i = 0; i < 5; i++)
    {
      st.sum_wphid[i] = 0.0f;
      st.sum_wphiq[i] = 0.0f;
    }
    st.start_sum_id = 0.0f;
    st.start_sum_iq = 0.0f;
    st.start_cnt = 0u;
    st.end_sum_id = 0.0f;
    st.end_sum_iq = 0.0f;
    st.end_cnt = 0u;
  }

  /* ∫(v - Rs0*i)dt */
  float yd = vd - Rs0 * id;
  float yq = vq - Rs0 * iq;
  st.sum_bd += yd * ONLINE_MTPA_TS_CTRL;
  st.sum_bq += yq * ONLINE_MTPA_TS_CTRL;

  /* ∫ we*phi dt */
  float phi_d[5], phi_q[5];
  basis_5(id, iq, phi_d, phi_q);
  float wdt = we * ONLINE_MTPA_TS_CTRL;
  for (int i = 0; i < 5; i++)
  {
    st.sum_wphid[i] += wdt * phi_d[i];
    st.sum_wphiq[i] += wdt * phi_q[i];
  }

  /* 起点均值：前 edgeN 点 */
  if (st.cnt < edgeN)
  {
    st.start_sum_id += id;
    st.start_sum_iq += iq;
    st.start_cnt++;
  }

  /* 终点均值：后 edgeN 点 */
  if (st.cnt >= (uint16_t)(M - edgeN))
  {
    st.end_sum_id += id;
    st.end_sum_iq += iq;
    st.end_cnt++;
  }

  st.cnt++;

  if (st.cnt < M) return;

  /* ----- bin结束：构造 Ad/Aq，更新EWLS ----- */
  float id_start = (st.start_cnt > 0u) ? (st.start_sum_id / (float)st.start_cnt) : id;
  float iq_start = (st.start_cnt > 0u) ? (st.start_sum_iq / (float)st.start_cnt) : iq;
  float id_end = (st.end_cnt > 0u) ? (st.end_sum_id / (float)st.end_cnt) : id;
  float iq_end = (st.end_cnt > 0u) ? (st.end_sum_iq / (float)st.end_cnt) : iq;

  if (g_freezeOnSteady)
  {
    float did = fabsf(id_end - id_start);
    float diq = fabsf(iq_end - iq_start);
    if ((did < g_steadyDidTh) && (diq < g_steadyDiqTh))
    {
      g_identifying = 0u;
      abort_current_bin();
      return;
    }
  }

  float phi_ds[5], phi_qs[5], phi_de[5], phi_qe[5];
  basis_5(id_start, iq_start, phi_ds, phi_qs);
  basis_5(id_end, iq_end, phi_de, phi_qe);

  float Ad[5], Aq[5];
  for (int i = 0; i < 5; i++)
  {
    float dphi_d = phi_de[i] - phi_ds[i];
    float dphi_q = phi_qe[i] - phi_qs[i];
    Ad[i] = dphi_d - st.sum_wphiq[i];
    Aq[i] = dphi_q + st.sum_wphid[i];
  }

  float bd = st.sum_bd;
  float bq = st.sum_bq;

  /* 当前bin贡献 AtA_new / Atb_new */
  float AtA_new[5][5] = {0};
  float Atb_new[5] = {0};

  for (int r = 0; r < 5; r++)
  {
    Atb_new[r] = Ad[r] * bd + Aq[r] * bq;
    for (int c = 0; c < 5; c++)
    {
      AtA_new[r][c] = Ad[r] * Ad[c] + Aq[r] * Aq[c];
    }
  }

  /* EWLS：S = beta*S + AtA_new ; t = beta*t + Atb_new */
  float beta = clampf(g_beta, 0.0f, 1.0f);
  for (int r = 0; r < 5; r++)
  {
    st.t[r] = beta * st.t[r] + Atb_new[r];
    for (int c = 0; c < 5; c++)
    {
      st.S[r][c] = beta * st.S[r][c] + AtA_new[r][c];
    }
  }

  st.bins_acc++;
  st.solve_cnt++;

  /* 解算间隔可调：每 g_solveDecim 个 bin 才解一次 */
  uint16_t dec = g_solveDecim;
  if (dec < 1u) dec = 1u;

  if (st.bins_acc >= 40u && (st.solve_cnt >= dec))
  {
    st.solve_cnt = 0u;

    /* 更新平滑系数（binSamples 或 fc 可被上位机改） */
    update_alpha();
    update_output_alpha_10k();

    /* Sreg = S + lambda*I */
    float Sreg[5][5];
    memcpy(Sreg, st.S, sizeof(Sreg));
    float lam = (g_lambda > 0.0f) ? g_lambda : 0.0f;
    for (int i = 0; i < 5; i++) Sreg[i][i] += lam;

    /* 固定 ld,lq 时：解 3x3（ldd,lqq,lc）更省 */
    if (g_fixLdLq)
    {
      /* x = [ld, ldd, lq, lqq, lc]
         固定 idxF=[0,2]，未知 idxU=[1,3,4] */
      float ld_fix = g_ld;
      float lq_fix = g_lq;

      float Suu[3][3];
      float rhs[3];

      /* rhs = t_u - S_uf * x_fixed */
      /* idxU:1,3,4 ; idxF:0,2 */
      int iu[3] = {1, 3, 4};
      int ifx[2] = {0, 2};
      float xfix[2] = {ld_fix, lq_fix};

      for (int r = 0; r < 3; r++)
      {
        int rr = iu[r];
        float tmp = st.t[rr];
        for (int k = 0; k < 2; k++)
        {
          int cc = ifx[k];
          tmp -= Sreg[rr][cc] * xfix[k];
        }
        rhs[r] = tmp;

        for (int c = 0; c < 3; c++)
        {
          int cc = iu[c];
          Suu[r][c] = Sreg[rr][cc];
        }
      }

      float xu[3];
      if (solve_spd_chol3(Suu, rhs, xu))
      {
        st.xhat[0] = ld_fix;
        st.xhat[2] = lq_fix;
        st.xhat[1] = xu[0]; /* ldd */
        st.xhat[3] = xu[1]; /* lqq */
        st.xhat[4] = xu[2]; /* lc  */
      }
      /* else: 保持上一次 st.xhat 不动 */
    }
    else
    {
      float xnew[5];
      if (solve_spd_chol5(Sreg, st.t, xnew))
      {
        for (int i = 0; i < 5; i++) st.xhat[i] = xnew[i];
      }
    }

    /* 参数平滑输出（低通） */
    float a = clampf(st.alpha, 0.0f, 1.0f);
    for (int i = 0; i < 5; i++)
    {
      st.xhat_f[i] += a * (st.xhat[i] - st.xhat_f[i]);
    }

    /* 写回全局参数（上位机也可读） */
    /* output params are published by output_filter_step_10k() */
  }

  /* bin复位 */
  abort_current_bin();
}

/* ---------------- MTPA 求解（解析） ---------------- */
static uint8_t cand_from_x(float x, float I, float* gamma, float* Tproxy)
{
  if (!isfinite(x) || fabsf(x) > 1.0f) return 0;

  float g = 0.5f * acosf(clampf(x, -1.0f, 1.0f));
  g = clampf(g, g_gammaMin, g_gammaMax);

  float id = I * COS(g);
  float iq = I * SIN(g);

  float psid, psiq;
  onlineMTPA_flux_from_i(id, iq, &psid, &psiq);
  float T = psid * iq - psiq * id;

  *gamma = g;
  *Tproxy = T;
  return 1;
}

void onlineMTPA_mtpa_for_Is(float Is_cmd, float* id_ref, float* iq_ref, float* gamma_deg,
                            uint8_t* ok)
{
  float I = fabsf(Is_cmd);
  int sgn = (Is_cmd >= 0.0f) ? 1 : -1;

  if (I < 1e-6f)
  {
    if (id_ref) *id_ref = 0.0f;
    if (iq_ref) *iq_ref = 0.0f;
    if (gamma_deg) *gamma_deg = 45.0f;
    if (ok) *ok = 1u;
    s_gamma_prev_valid = 0u;
    s_mtpa_theta_valid = 0u;
    g_mtpaThetaRad = 0.0f;
    g_mtpaIqSign = 1;
    return;
  }

  float ld = g_ld;
  float ldd = g_ldd;
  float lq = g_lq;
  float lqq = g_lqq;
  float lc = g_lc;

  float I2 = I * I;
  float I4 = I2 * I2;

  float K = (ld - lq) * I2 + 0.5f * (ldd - lqq) * I4;
  float M = (0.5f * (ldd + lqq) - lc) * I4;

  float gamma = 0.25f * 3.1415926f; /* fallback */
  uint8_t good = 1u;

  if (fabsf(M) < 1e-12f)
  {
    gamma = clampf(0.25f * 3.1415926f, g_gammaMin, g_gammaMax);
  }
  else
  {
    float disc = K * K + 8.0f * M * M;
    float sd = sqrtf(disc);

    float x1 = (-K + sd) / (4.0f * M);
    float x2 = (-K - sd) / (4.0f * M);

    float g1 = 0, T1 = 0, g2 = 0, T2 = 0;
    uint8_t ok1 = cand_from_x(x1, I, &g1, &T1);
    uint8_t ok2 = cand_from_x(x2, I, &g2, &T2);

    if (!ok1 && !ok2)
    {
      gamma = clampf(0.25f * 3.1415926f, g_gammaMin, g_gammaMax);
      good = 0u;
    }
    else if (ok1 && !ok2)
      gamma = g1;
    else if (!ok1 && ok2)
      gamma = g2;
    else
    {
      /* 规则：优先选择 >=45deg 的候选；否则取T更大 */
      float g45 = 0.25f * 3.1415926f;
      uint8_t s1 = (g1 >= g45);
      uint8_t s2 = (g2 >= g45);
      if (s1 && !s2)
        gamma = g1;
      else if (!s1 && s2)
        gamma = g2;
      else
        gamma = (T1 >= T2) ? g1 : g2;
    }
  }

  float stepMaxDeg = g_gammaStepMaxDeg;
  if ((stepMaxDeg > 0.0f) && s_gamma_prev_valid)
  {
    float dgamma = gamma - s_gamma_prev;
    float dgammaMax = stepMaxDeg * (3.1415926f / 180.0f);
    if (dgamma > dgammaMax)
      gamma = s_gamma_prev + dgammaMax;
    else if (dgamma < -dgammaMax)
      gamma = s_gamma_prev - dgammaMax;
    gamma = clampf(gamma, g_gammaMin, g_gammaMax);
  }
  s_gamma_prev = gamma;
  s_gamma_prev_valid = 1u;
  s_mtpa_theta_valid = 1u;
  g_mtpaThetaRad = gamma;
  g_mtpaIqSign = (sgn >= 0) ? 1 : -1;

  float idr = I * COS(gamma);
  float iqr = (float)sgn * I * SIN(gamma);

  if (id_ref) *id_ref = idr;
  if (iq_ref) *iq_ref = iqr;
  if (gamma_deg) *gamma_deg = gamma * (180.0f / 3.1415926f);
  if (ok) *ok = good;
}
