

#ifndef ONLINE_MTPA_H
#define ONLINE_MTPA_H

#include <stdint.h>
#include "arm_math.h" /* CMSIS-DSP math */  // IWYU pragma: export
#define COS(x) arm_cos_f32(x)
#define SIN(x) arm_sin_f32(x)

#ifdef __cplusplus
extern "C"
{
#endif

/* ---------------- 用户可改的宏 ---------------- */
#ifndef ONLINE_MTPA_MAX_BIN_SAMPLES
#define ONLINE_MTPA_MAX_BIN_SAMPLES (64u) /* 运行时允许的最大bin采样点数 */
#endif

/* 10kHz 调用周期：默认 1e-4s（你系统是10k电流环） */
#ifndef ONLINE_MTPA_TS_CTRL
#define ONLINE_MTPA_TS_CTRL (1.0e-4f)
#endif

/* flag 定义（step_10k 的 flags 输入） */
#define ONLINE_MTPA_FLAG_VSAT (0x01u)  /* 电压饱和/过调制，建议冻结辨识 */
#define ONLINE_MTPA_FLAG_BAD_V (0x02u) /* 电压不可信（例如未用实际施加电压） */

/* dT/dtheta 冻结判据选择 */
#define ONLINE_MTPA_TORQUE_GRAD_FREEZE_OFF (0x00u)
#define ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT (0x01u)
#define ONLINE_MTPA_TORQUE_GRAD_FREEZE_DT_OVER_T (0x02u)

  /* ---------------- 参数与配置（上位机可写） ---------------- */
  /* 识别得到的模型参数（单位见注释） */
  extern volatile float g_ld;  /* H */
  extern volatile float g_ldd; /* H/A^2 */
  extern volatile float g_lq;  /* H */
  extern volatile float g_lqq; /* H/A^2 */
  extern volatile float g_lc;  /* H/A^2 */

  /* EWLS 配置 */
  extern volatile float g_beta;           /* 遗忘因子 0~1，推荐0.98 */
  extern volatile float g_lambda;         /* 岭正则，推荐 1e-7~1e-6 */
  extern volatile float g_smoothFcHz;     /* 参数输出平滑截止(Hz)，推荐2~5 */
  extern volatile float g_outputFcHz10k;  /* 10kHz parameter output low-pass cutoff (Hz) */
  extern volatile uint16_t g_binSamples;  /* 一个bin包含多少个10k采样点，默认20 -> 2ms */
  extern volatile uint16_t g_edgeSamples; /* 起止均值点数(<=binSamples)，默认5 */
  extern volatile uint16_t g_solveDecim;  /* 每多少个bin求解一次，默认5 -> 10ms */
  extern volatile uint8_t g_fixLdLq;      /* 1:固定ld,lq；0:一起估 */

  /* 门控/冻结 */
  extern volatile uint8_t g_identEnable;  /* 1:允许辨识；0:冻结 */
  extern volatile float g_identIsMin;     /* A，电流幅值门限，小于冻结 */
  extern volatile float g_identWeMin;     /* rad/s，电角速度门限，小于冻结 */
  extern volatile uint8_t g_freezeOnVsat; /* 1:电压饱和冻结 */
  extern volatile uint8_t g_freezeOnSteady; /* 1: freeze identification when a bin is near steady state */
  extern volatile float g_steadyDidTh;               /* A, |id_end-id_start| threshold */
  extern volatile float g_steadyDiqTh;               /* A, |iq_end-iq_start| threshold */
  extern volatile uint8_t g_torqueGradFreezeSel;     /* 0:关 1:|dT/dtheta| 2:|dT/dtheta/T| */
  extern volatile uint8_t Sel12;                     /* 计数到阈值后，将 g_torqueGradFreezeSel 置为 Sel12(1或2) */
  extern volatile uint32_t g_torqueGradFreezeSelDelayCnt; /* 启动辨识后，切换到 Sel12 的10k计数阈值 */
  extern volatile float g_torqueGradTh;              /* 转矩偏导阈值 */
  extern volatile float g_torqueGradNormTh;          /* 归一化转矩偏导阈值 */
  extern volatile float g_torqueGradFiltFcHz;        /* 转矩偏导滤波截止频率(Hz) */
  extern volatile float g_torqueGradRestartDeltaIs;  /* A，重启辨识电流增量阈值 */
  extern volatile uint8_t g_identFrozenByTorqueGrad; /* 1: 因转矩偏导被冻结 */
  extern volatile uint8_t g_identifying;             /* 1: 当前周期正在辨识采集 */
  extern volatile float g_mtpaThetaRad;              /* MTPA 输出电流角(rad) */
  extern volatile int8_t g_mtpaIqSign;               /* MTPA 输出 iq 符号(+1/-1) */
  extern volatile float g_torqueProxyNow;            /* 当前角度下转矩代理值 */
  extern volatile float g_torqueGradRawNow;          /* 当前角度下 dT/dtheta 原始值 */
  extern volatile float g_torqueGradFiltNow;         /* 滤波后的 dT/dtheta */
  extern volatile float g_torqueGradNormNow;         /* |dT/dtheta|/|T| */

  /* MTPA 配置 */
  extern volatile float gamma_deg;
  extern volatile float g_gammaMin;        /* rad，默认1e-3 */
  extern volatile float g_gammaMax;        /* rad，默认pi/2-1e-3 */
  extern volatile float g_gammaStepMaxDeg; /* deg/call, <=0 means disabled */
  extern volatile uint8_t enable_45;

  /* ---------------- 状态查询 ---------------- */
  typedef struct
  {
    uint32_t bins_accumulated; /* 自 reset 起累计的 bin 数 */
    uint8_t valid;             /* 1: 已达到最小bin数，参数可用 */
  } onlineMTPA_status_t;

  /* ---------------- API ---------------- */
  void onlineMTPA_init(void);
  void onlineMTPA_reset_ident(void);

  /* 10kHz 调用：把你控制回路已有的 vd/vq/id/iq/we 直接传入
   * vd,vq: dq电压（尽量是实际施加电压，V）
   * id,iq: dq电流（A）
   * we:   电角速度（rad/s）
   * Rs0:  用于构造 y = v - Rs0*i 的基准电阻（Ω），不在此文件内估计
   * flags: ONLINE_MTPA_FLAG_*，用于冻结门控等
   */
  void onlineMTPA_step_10k(float vd, float vq, float id, float iq, float we, float Rs0,
                           uint8_t flags);

  /* 查询状态 */
  onlineMTPA_status_t onlineMTPA_get_status(void);

  /* 用当前参数计算磁链（psi(i) 模型） */
  void onlineMTPA_flux_from_i(float id, float iq, float* psid, float* psiq);

  /* 计算增量电感/互感（基于 psi(i) 模型的雅可比） */
  void onlineMTPA_incL_from_i(float id, float iq, float* Ldd, float* Lqq, float* Ldq);

  /* 转矩代理（不含 3/2*p 系数），用于比较/选择候选点 */
  float onlineMTPA_torque_proxy(float id, float iq);

  /* 基于当前辨识模型参数，在固定 Is 下按电流角直接解析计算 T 与 dT/dtheta */
  void onlineMTPA_torque_grad_from_theta(float Is, float theta_rad, int8_t iq_sign, float* Tproxy,
                                         float* dTdTheta);

  /* 固定电流幅值 Is 的 MTPA 求解（解析+候选点筛选）
   * Is_cmd: 电流幅值指令（A），允许带符号：Is_cmd<0 -> 输出iq_ref<0
   * 输出:
   *  gamma_deg: 电流角（deg，0~90）
   *  id_ref, iq_ref
   *  ok: 1有效，0退化
   */
  void onlineMTPA_mtpa_for_Is(float Is_cmd, float* id_ref, float* iq_ref, float* gamma_deg,
                              uint8_t* ok);

#ifdef __cplusplus
}
#endif
#endif
