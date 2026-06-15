#include "control_profile.h"

#include <string.h>

#define CONTROL_PROFILE_SPEED_WARMUP (100u)
#define CONTROL_PROFILE_P99_UPDATE_PERIOD (64u)
#define SYSTEM_CORE_CLOCK_MHZ (120u)

volatile uint8_t g_mtpa_method_select = MTPA_METHOD_ANALYTIC_ONLINE;
volatile uint8_t g_prof_active_method = MTPA_METHOD_ANALYTIC_ONLINE;
volatile uint8_t g_prof_force_is_cmd_enable = 0u;
volatile float g_prof_force_is_cmd = 10.0f;
volatile uint8_t g_prof_speed_cycle_hit = 0u;
volatile uint8_t g_prof_ls_solve_hit = 0u;
volatile uint32_t g_prof_last_isr_cycles = 0u;
volatile uint32_t g_prof_last_ctrl_cycles = 0u;
volatile uint32_t g_prof_last_comm_cycles = 0u;
volatile uint8_t g_prof_timer_src = 0u;
volatile uint8_t g_prof_dwt_ok = 0u;
volatile uint32_t g_prof_dwt_demcr = 0u;
volatile uint32_t g_prof_dwt_ctrl = 0u;
volatile uint32_t g_prof_dwt_cyccnt = 0u;
volatile uint32_t g_prof_tmr1_tick_cycles = 6u;

ControlProfileStats_t g_isr_profile[MTPA_METHOD_COUNT][PROFILE_METRIC_COUNT];
static uint32_t s_profile_hist[MTPA_METHOD_COUNT][PROFILE_METRIC_COUNT][CONTROL_PROFILE_HIST_BINS];
static uint8_t s_profile_seen_method = MTPA_METHOD_ANALYTIC_ONLINE;

static uint8_t sanitize_method(uint8_t method)
{
  if (method >= MTPA_METHOD_COUNT) return MTPA_METHOD_ANALYTIC_ONLINE;
  return method;
}

static uint32_t hist_bin_to_cycles(uint32_t bin)
{
  return (bin + 1u) * CONTROL_PROFILE_HIST_BIN_CYCLES;
}

static uint32_t cycles_to_hist_bin(uint32_t cycles)
{
  uint32_t bin = cycles / CONTROL_PROFILE_HIST_BIN_CYCLES;
  if (bin >= CONTROL_PROFILE_HIST_BINS) bin = CONTROL_PROFILE_HIST_BINS - 1u;
  return bin;
}

static uint32_t cycles_to_us_x100(uint32_t cycles)
{
  return (uint32_t)(((uint64_t)cycles * 100u + (SYSTEM_CORE_CLOCK_MHZ / 2u)) /
                    SYSTEM_CORE_CLOCK_MHZ);
}

static void update_p99(uint8_t method, ControlProfileMetric_t metric)
{
  ControlProfileStats_t* stats = &g_isr_profile[method][metric];
  uint32_t target = (stats->count * 99u + 99u) / 100u;
  uint32_t cumulative = 0u;
  uint32_t p99_estimate = 0u;

  if (target == 0u)
  {
    stats->p99_cycles = 0u;
    return;
  }

  for (uint32_t i = 0u; i < CONTROL_PROFILE_HIST_BINS; i++)
  {
    cumulative += s_profile_hist[method][metric][i];
    if (cumulative >= target)
    {
      p99_estimate = hist_bin_to_cycles(i);
      if (p99_estimate > stats->max_cycles)
      {
        p99_estimate = stats->max_cycles;
      }
      stats->p99_cycles = p99_estimate;
      return;
    }
  }

  p99_estimate = hist_bin_to_cycles(CONTROL_PROFILE_HIST_BINS - 1u);
  if (p99_estimate > stats->max_cycles)
  {
    p99_estimate = stats->max_cycles;
  }
  stats->p99_cycles = p99_estimate;
}

void ControlProfile_ResetMethod(uint8_t method)
{
  method = sanitize_method(method);

  for (uint32_t metric = 0u; metric < PROFILE_METRIC_COUNT; metric++)
  {
    ControlProfileStats_t* stats = &g_isr_profile[method][metric];
    stats->count = 0u;
    stats->warmup_left = (metric == PROFILE_METRIC_SPEED_CYCLE_ISR) ? CONTROL_PROFILE_SPEED_WARMUP : 0u;
    stats->latest_cycles = 0u;
    stats->min_cycles = 0xFFFFFFFFu;
    stats->max_cycles = 0u;
    stats->p99_cycles = 0u;
    stats->mean_cycles = 0u;
    stats->mean_us_x100 = 0u;
    stats->sum_cycles = 0u;
    memset(s_profile_hist[method][metric], 0, sizeof(s_profile_hist[method][metric]));
  }
}

void ControlProfile_ResetAll(void)
{
  for (uint8_t method = 0u; method < MTPA_METHOD_COUNT; method++)
  {
    ControlProfile_ResetMethod(method);
  }
  s_profile_seen_method = sanitize_method(g_mtpa_method_select);
  g_prof_active_method = s_profile_seen_method;
}

void ControlProfile_ServiceMethodSelect(void)
{
  uint8_t requested = sanitize_method(g_mtpa_method_select);
  if (requested != g_mtpa_method_select)
  {
    g_mtpa_method_select = requested;
  }

  if (requested != s_profile_seen_method)
  {
    ControlProfile_ResetMethod(requested);
    s_profile_seen_method = requested;
  }

  g_prof_active_method = requested;
}

uint8_t ControlProfile_GetActiveMethod(void)
{
  return sanitize_method(g_prof_active_method);
}

void ControlProfile_Record(uint8_t method, ControlProfileMetric_t metric, uint32_t cycles)
{
  method = sanitize_method(method);
  if (metric >= PROFILE_METRIC_COUNT) return;

  ControlProfileStats_t* stats = &g_isr_profile[method][metric];

  if (stats->warmup_left > 0u)
  {
    stats->warmup_left--;
    return;
  }

  stats->latest_cycles = cycles;
  stats->count++;
  stats->sum_cycles += cycles;
  stats->mean_cycles = (uint32_t)(stats->sum_cycles / stats->count);
  stats->mean_us_x100 = cycles_to_us_x100(stats->mean_cycles);
  if (cycles < stats->min_cycles) stats->min_cycles = cycles;
  if (cycles > stats->max_cycles) stats->max_cycles = cycles;

  s_profile_hist[method][metric][cycles_to_hist_bin(cycles)]++;
  if ((stats->count == 1u) || ((stats->count % CONTROL_PROFILE_P99_UPDATE_PERIOD) == 0u))
  {
    update_p99(method, metric);
  }
}

uint32_t ControlProfile_CyclesToUsX100(uint32_t cycles)
{
  return cycles_to_us_x100(cycles);
}
