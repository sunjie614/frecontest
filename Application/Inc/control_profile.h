#ifndef CONTROL_PROFILE_H
#define CONTROL_PROFILE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CONTROL_PROFILE_HIST_BINS (256u)
#define CONTROL_PROFILE_HIST_BIN_CYCLES (512u)

typedef enum
{
  /* 离线 LUT 查表 MTPA 路径。 */
  MTPA_METHOD_LUT = 0u,
  /* 当前解析在线 MTPA 路径，默认方法。 */
  MTPA_METHOD_ANALYTIC_ONLINE = 1u,
  /* 基于局部模型的固定迭代数值 MTPA 路径。 */
  MTPA_METHOD_LOCAL_NUMERIC = 2u,
  MTPA_METHOD_COUNT = 3u
} MTPA_Method_t;

typedef enum
{
  /* 每个 ADC 控制中断周期，不含 justfloat 通信耗时。 */
  PROFILE_METRIC_ALL_ISR = 0u,
  /* 执行 1 kHz 速度环那一拍的 ADC 控制中断周期。 */
  PROFILE_METRIC_SPEED_CYCLE_ISR = 1u,
  /* 触发在线最小二乘求解那一拍的 ADC 控制中断周期。 */
  PROFILE_METRIC_LS_SOLVE_ISR = 2u,
  PROFILE_METRIC_COUNT = 3u
} ControlProfileMetric_t;

typedef struct
{
  /* 预热丢弃结束后的有效统计次数。 */
  volatile uint32_t count;
  /* 统计开始前还需要丢弃的样本数。 */
  volatile uint32_t warmup_left;
  /* 最近一次记录的 ISR 周期，单位为 DWT CPU cycles。 */
  volatile uint32_t latest_cycles;
  /* 已记录 ISR 周期最小值，单位为 DWT CPU cycles。 */
  volatile uint32_t min_cycles;
  /* 已记录 ISR 周期最大值，单位为 DWT CPU cycles。 */
  volatile uint32_t max_cycles;
  /* 由周期直方图估算的 99 分位 ISR 周期。 */
  volatile uint32_t p99_cycles;
  /* 已记录样本的平均 ISR 周期，单位为 cycles，A2L 直接读取。 */
  volatile uint32_t mean_cycles;
  /* 已记录样本的平均 ISR 时间，单位为 us*100，例如 1234 表示 12.34us。 */
  volatile uint32_t mean_us_x100;
  /* 已记录周期总和；mean_cycles = sum_cycles / count。 */
  volatile uint64_t sum_cycles;
} ControlProfileStats_t;

/* A2L 可写的模式请求值：0=LUT，1=解析在线 MTPA，2=局部模型数值 MTPA。 */
extern volatile uint8_t g_mtpa_method_select;
/* 自动切换服务处理后，FOC 和 ISR 统计实际正在使用的方法。 */
extern volatile uint8_t g_prof_active_method;
/* A2L 可写：置 1 后，速度环 MTPA 使用 g_prof_force_is_cmd 作为测试电流指令。 */
extern volatile uint8_t g_prof_force_is_cmd_enable;
/* A2L 可写：profiling 测试用电流幅值指令，单位 A，允许正负号。 */
extern volatile float g_prof_force_is_cmd;
/* 仅在执行 1 kHz 速度环的 ADC ISR 周期内置 1。 */
extern volatile uint8_t g_prof_speed_cycle_hit;
/* 仅在触发在线最小二乘参数求解的 ADC ISR 周期内置 1。 */
extern volatile uint8_t g_prof_ls_solve_hit;
/* 最近一次原始 ADC ISR 周期，包含 justfloat 通信，单位为 DWT cycles。 */
extern volatile uint32_t g_prof_last_isr_cycles;
/* 最近一次 ADC 控制 ISR 周期，不含 justfloat 通信，单位为 DWT cycles。 */
extern volatile uint32_t g_prof_last_ctrl_cycles;
/* 最近一次 ADC ISR 内 justfloat 通信耗时，单位为 DWT cycles。 */
extern volatile uint32_t g_prof_last_comm_cycles;
/* 当前 profiling 计时源：0=DWT->CYCCNT，1=TIMER1 自动兜底。 */
extern volatile uint8_t g_prof_timer_src;
/* DWT 是否可用且本次 ISR 内计数有变化：1=可用，0=不可用或未走。 */
extern volatile uint8_t g_prof_dwt_ok;
/* 最近一次读取到的 CoreDebug->DEMCR 原值，用于确认 TRCENA 是否打开。 */
extern volatile uint32_t g_prof_dwt_demcr;
/* 最近一次读取到的 DWT->CTRL 原值，用于确认 CYCCNTENA 是否打开。 */
extern volatile uint32_t g_prof_dwt_ctrl;
/* 最近一次读取到的 DWT->CYCCNT 原值，用于确认计数器是否在变化。 */
extern volatile uint32_t g_prof_dwt_cyccnt;
/* TIMER1 每个 tick 折算的 CPU cycles，当前 TIM1_Init() 下应为 6。 */
extern volatile uint32_t g_prof_tmr1_tick_cycles;
/* A2L 可见的主统计表：g_isr_profile[method][metric]。 */
extern ControlProfileStats_t g_isr_profile[MTPA_METHOD_COUNT][PROFILE_METRIC_COUNT];

/* 清空所有方法/指标统计，并按 g_mtpa_method_select 设置当前生效方法。 */
void ControlProfile_ResetAll(void);
/* 清空某一种方法的统计；通常由模式切换服务自动调用。 */
void ControlProfile_ResetMethod(uint8_t method);
/* 在主循环调用：应用 g_mtpa_method_select 的变化，并清空新方法统计。 */
void ControlProfile_ServiceMethodSelect(void);
/* 读取经过合法性处理的当前生效方法，供控制代码和 ISR 分类使用。 */
uint8_t ControlProfile_GetActiveMethod(void);
/* 为指定方法/指标记录一次 ISR 耗时样本。 */
void ControlProfile_Record(uint8_t method, ControlProfileMetric_t metric, uint32_t cycles);
/* 120 MHz 下把 DWT cycles 转换为 us * 100，例如 1234 表示 12.34 us。 */
uint32_t ControlProfile_CyclesToUsX100(uint32_t cycles);

#ifdef __cplusplus
}
#endif

#endif /* CONTROL_PROFILE_H */
