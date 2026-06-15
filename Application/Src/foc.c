#include "foc.h"
#include "MTPA.h"
#include "control_profile.h"
#include "hardware_interface.h"
#include "identification.h"
#include "onlineMTPA.h"
#include "position_sensor.h"

typedef struct
{
  float a;       // 反馈系数（= 极点位置）
  float y_last;  // 上一次输出值
} LowPassFilter_t;

Motor_Parameter_t Motor;
FOC_Parameter_t FOC;
VF_Parameter_t VF;
IF_Parameter_t IF;
PID_Controller_t Id_PID;
PID_Controller_t Iq_PID;
PID_Controller_t Speed_PID;
RampGenerator_t Speed_Ramp;
InvPark_t Inv_Park;
Clarke_t Clarke;
MTPA_Coefficients_t MTPA;
float Electric_Power = 0.0F;
float Electric_Te = 0.0F;
LowPassFilter_t Power_Filter = {.a = 0.95442F, .y_last = 0.0F};

float theta_mech = 0.0F;
float theta_elec = 0.0F;
float theta_factor = 0.0F;  // Sensor data to mechanic angle conversion factor

float Speed_Ref = 0.0F;
FluxExperiment_t Experiment = {0};
/* 在线辨识支路滤波观测量（仅观测/传输，不参与控制） */
volatile float g_identInputLpfFcHz = 800.0f;
volatile float g_identUdFilt = 0.0f;
volatile float g_identUqFilt = 0.0f;
volatile float g_identIdFilt = 0.0f;
volatile float g_identIqFilt = 0.0f;

static inline float Get_Theta(float Freq, float Theta);
static inline void Parameter_Init(void);
static inline void Theta_Process(void);
static inline float wrap_theta_2pi(float theta);
static inline float RampGenerator(RampGenerator_t* ramp);
static inline void PID_Controller(float Ref, float Feedback, PID_Controller_t* PID_Controller);
static inline void ClarkeTransform(float_t Ia, float_t Ib, float_t Ic, Clarke_t* out);
static inline void ParkTransform(float_t Ialpha, float_t Ibeta, float_t theta,
                                 FOC_Parameter_t* out);
static inline void InvParkTransform(float_t Udaxis, float_t Uqaxis, float_t theta, InvPark_t* out);
static inline void SVPWM_Generate(float Ualpha, float Ubeta, float inv_Vdc, FOC_Parameter_t* foc);
static inline float Cal_Power(FOC_Parameter_t* foc);
static inline float LowPassFilter_Update(LowPassFilter_t* filter, float x);
static inline float sat(float x);
static void LocalModelNumericMTPA(float Is_cmd, float* id_ref, float* iq_ref, float* gamma_deg);

static float LocalModelTorqueAtGamma(float Is, float gamma)
{
  float id = Is * COS(gamma);
  float iq = Is * SIN(gamma);
  return onlineMTPA_torque_proxy(id, iq);
}

static void LocalModelNumericMTPA(float Is_cmd, float* id_ref, float* iq_ref, float* gamma_deg_out)
{
  const uint8_t iter_count = 12u;
  const float pi = 3.1415926f;
  const float gr = 0.61803398875f;
  float Is = fabsf(Is_cmd);
  float sign = (Is_cmd >= 0.0f) ? 1.0f : -1.0f;

  if (Is < 1e-6f)
  {
    if (id_ref) *id_ref = 0.0f;
    if (iq_ref) *iq_ref = 0.0f;
    if (gamma_deg_out) *gamma_deg_out = 45.0f;
    return;
  }

  float lo = g_gammaMin;
  float hi = g_gammaMax;
  if (lo < 0.0f) lo = 0.0f;
  if (hi > 0.5f * pi) hi = 0.5f * pi;
  if (hi <= lo)
  {
    lo = 1e-3f;
    hi = 0.5f * pi - 1e-3f;
  }

  float c = hi - (hi - lo) * gr;
  float d = lo + (hi - lo) * gr;
  float tc = LocalModelTorqueAtGamma(Is, c);
  float td = LocalModelTorqueAtGamma(Is, d);

  for (uint8_t i = 0u; i < iter_count; i++)
  {
    if (tc < td)
    {
      lo = c;
      c = d;
      tc = td;
      d = lo + (hi - lo) * gr;
      td = LocalModelTorqueAtGamma(Is, d);
    }
    else
    {
      hi = d;
      d = c;
      td = tc;
      c = hi - (hi - lo) * gr;
      tc = LocalModelTorqueAtGamma(Is, c);
    }
  }

  float gamma = 0.5f * (lo + hi);
  if (id_ref) *id_ref = Is * COS(gamma);
  if (iq_ref) *iq_ref = sign * Is * SIN(gamma);
  if (gamma_deg_out) *gamma_deg_out = gamma * (180.0f / pi);
}
// SECTION - FOC Main
void FOC_Main(void)
{
  Theta_Process();
  ClarkeTransform(FOC.Ia, FOC.Ib, FOC.Ic, &Clarke);

  switch (FOC.Mode)
  {
    case INIT:
    {
      Parameter_Init();
      FOC.Mode = IDLE;
      MTPA_build_table(mtpa_table, MTPA_TABLE_POINTS, 0.0f, 50.0f); /* T 从 0 到 50, 共 51 点 */
      if (Experiment.Complete == true)
      {
        float ad0 = 0.0F, add = 0.0F, aq0 = 0.0F, aqq = 0.0F, adq = 0.0F;
        Get_Identification_Results(&Experiment, &ad0, &add, &aq0, &aqq, &adq);
        MTPA_Get_Parameter(ad0, add, aq0, aqq, adq);
        MTPA_build_table(mtpa_table, MTPA_TABLE_POINTS, 0.0f, 50.0f); /* T 从 0 到 50, 共 51 点 */
      }
      Experiment_Init(&Experiment, FOC.Ts, 512, 2, 20, 2, 10, 1, 200);
      onlineMTPA_init();

      break;
    }
    case IDLE:
    {
      STOP = 1;
      break;
    }
    case VF_MODE:
    {
      VF.Theta = Get_Theta(VF.Freq, VF.Theta);
      FOC.Theta = VF.Theta;
      FOC.Ud_ref = VF.Vref_Ud;
      FOC.Uq_ref = VF.Vref_Uq;
      ParkTransform(Clarke.Ialpha, Clarke.Ibeta, FOC.Theta, &FOC);
      break;
    }
    // SECTION - IF Mode
    case IF_MODE:
    {
      if (IF.Sensor_State == Enable)
      {
        IF.Theta = FOC.Theta;
      }
      else
      {
        IF.Theta = Get_Theta(IF.IF_Freq, IF.Theta);
      }
      FOC.Theta = IF.Theta;

      ParkTransform(Clarke.Ialpha, Clarke.Ibeta, FOC.Theta, &FOC);

      static volatile float Is_test = 0;
      static volatile float Is_theta_test = 0;
      static volatile bool Test_flag = 0;
      float Id_test = Is_test * cosf(Is_theta_test * M_2PI / 360.0f);
      float Iq_test = Is_test * sinf(Is_theta_test * M_2PI / 360.0f);
      if (Test_flag == 1)
      {
        IF.Id_ref = Id_test;
        IF.Iq_ref = Iq_test;
      }

      FOC.Id_ref = IF.Id_ref;
      FOC.Iq_ref = IF.Iq_ref;

      PID_Controller(FOC.Id_ref, FOC.Id, &Id_PID);
      PID_Controller(FOC.Iq_ref, FOC.Iq, &Iq_PID);

      FOC.Ud_ref = Id_PID.output;
      FOC.Uq_ref = Iq_PID.output;

      break;
    }
    // !SECTION
    // SECTION - Speed Mode
    case Speed:
    {
      ParkTransform(Clarke.Ialpha, Clarke.Ibeta, FOC.Theta, &FOC);

      static uint16_t Speed_Count = 0;
      Speed_Count++;
      if (Speed_Count > 9)
      {
        g_prof_speed_cycle_hit = 1u;
        Speed_Count = 0;
        Speed_Ramp.target = Speed_Ref;

        PID_Controller(RampGenerator(&Speed_Ramp), FOC.Speed, &Speed_PID);

        float id = FOC.Id;
        float iq = FOC.Iq;
        static volatile float Ldd_est = 0;
        static volatile float Lqq_est = 0;
        static volatile float Ldq_est = 0;
        static volatile float id_modleref = 0;
        static volatile float iq_modleref = 0;

        uint8_t ok;

        float Ldd_tmp = Ldd_est;
        float Lqq_tmp = Lqq_est;
        float Ldq_tmp = Ldq_est;
        onlineMTPA_incL_from_i(id, iq, &Ldd_tmp, &Lqq_tmp, &Ldq_tmp);
        Ldd_est = Ldd_tmp;
        Lqq_est = Lqq_tmp;
        Ldq_est = Ldq_tmp;
        // float32_t Is2 = FOC.Id * FOC.Id + FOC.Iq * FOC.Iq;
        // static volatile float32_t IS_cmd = 0;
        // SQRT(Is2, &IS_cmd);
        float IS_cmd = Speed_PID.output;
        if (g_prof_force_is_cmd_enable != 0u)
        {
          IS_cmd = g_prof_force_is_cmd;
        }
        if (enable_45 != 0.0f)
        {
          FOC.Iq_ref = IS_cmd * 0.70710678f;  // cos(45°) ≈ 0.7071
          if (IS_cmd >= 0)
          {
            FOC.Id_ref = IS_cmd * 0.70710678f;  //
          }
          else
          {
            FOC.Id_ref = -IS_cmd * 0.70710678f;  //
          }
        }
        else
        {
          float id_ref_calc = id_modleref;
          float iq_ref_calc = iq_modleref;
          float gamma_calc = gamma_deg;
          switch (ControlProfile_GetActiveMethod())
          {
            case MTPA_METHOD_LUT:
            {
              float iq_abs = fabsf(IS_cmd);
              float iq_out = 0.0f;
              MTPA_interp_by_Iq(mtpa_table, MTPA_TABLE_POINTS, iq_abs, &id_ref_calc, &iq_out);
              iq_ref_calc = IS_cmd;
              ok = 1u;
              break;
            }
            case MTPA_METHOD_LOCAL_NUMERIC:
            {
              LocalModelNumericMTPA(IS_cmd, &id_ref_calc, &iq_ref_calc, &gamma_calc);
              ok = 1u;
              break;
            }
            case MTPA_METHOD_ANALYTIC_ONLINE:
            default:
            {
              onlineMTPA_mtpa_for_Is(IS_cmd, &id_ref_calc, &iq_ref_calc, &gamma_calc, &ok);
              break;
            }
          }
          id_modleref = id_ref_calc;
          iq_modleref = iq_ref_calc;
          gamma_deg = gamma_calc;
          FOC.Iq_ref = iq_modleref;
          FOC.Id_ref = id_modleref;
        }
      }
      float we = FOC.Speed * M_2PI / 30.0f;  // 电气角速度 (rad/s)
      float ia = FOC.Ia, ib = FOC.Ib, ic = FOC.Ic;
      static volatile float K_dabc = 10.454F;  // 这个系数需要根据实际电压和电流范围进行调整
      float dua = 0.13635 * ia + sat(ia / 0.816936) * K_dabc;
      float dub = 0.13635 * ib + sat(ib / 0.816936) * K_dabc;
      float duc = 0.13635 * ic + sat(ic / 0.816936) * K_dabc;
      float theta = FOC.Theta;
      float COS_theta = COS(theta);
      float SIN_theta = SIN(theta);

      /*float Udin = FOC.Ud_ref;
      float Uqin = FOC.Uq_ref;
      float theta = FOC.Theta;
      // 1. 逆 Park 变换 (d,q -> α,β)

      float Ualphain = Udin * COS_theta - Uqin * SIN_theta;
      float Ubetain = Udin * SIN_theta + Uqin * COS_theta;
      // 2. 逆 Clark 变换 (α,β -> a,b,c)  —— 等幅值变换
      float Ua = Ualphain;
      float Ub = -0.5f * Ualphain + 0.86602540378f * Ubetain;  // 0.8660254 = sqrt(3)/2
      float Uc = -0.5f * Ualphain - 0.86602540378f * Ubetain;*/
      float da = 1.0f - FOC.Tcm1, db = 1.0f - FOC.Tcm2, dc = 1.0f - FOC.Tcm3, Vdc = FOC.Udc;
      static volatile float Uai = 0, Ubi = 0, Uci = 0, Udi = 0, Uqi = 0;
      Uai = (Vdc / 3.0f) * (2.0f * da - db - dc);
      Ubi = (Vdc / 3.0f) * (2.0f * db - da - dc);
      Uci = (Vdc / 3.0f) * (2.0f * dc - da - db);

      float Ua_in = -dua + Uai;
      float Ub_in = -dub + Ubi;
      float Uc_in = -duc + Uci;
      // 1. Clark 变换 (通用等幅值)
      float Ualpha = (2.0f / 3.0f) * (Ua_in - 0.5f * Ub_in - 0.5f * Uc_in);
      float Ubeta =
          (2.0f / 3.0f) * (0.86602540378f * Ub_in - 0.86602540378f * Uc_in);  // √3/2 ≈ 0.8660254
      // 2. Park 变换
      Udi = Ualpha * COS_theta + Ubeta * SIN_theta;
      Uqi = -Ualpha * SIN_theta + Ubeta * COS_theta;
      static volatile float Rs0 =
          0.65f;  // 用于构造 y = v - Rs0*i 的基准电阻（Ω），不在此文件内估计

      /* 仅对在线辨识输入做一阶低通，不影响控制环其它路径 */
      static LowPassFilter_t Ud_ident_lpf = {.a = 0.0f, .y_last = 0.0f};
      static LowPassFilter_t Uq_ident_lpf = {.a = 0.0f, .y_last = 0.0f};
      static LowPassFilter_t Id_ident_lpf = {.a = 0.0f, .y_last = 0.0f};
      static LowPassFilter_t Iq_ident_lpf = {.a = 0.0f, .y_last = 0.0f};
      static uint8_t ident_lpf_inited = 0u;

      float Udi_ident = Udi;
      float Uqi_ident = Uqi;
      float Id_ident = FOC.Id;
      float Iq_ident = FOC.Iq;

      float fc_ident = g_identInputLpfFcHz;
      if (fc_ident > 0.0f)
      {
        /* 一阶离散低通：y = a*y(k-1) + (1-a)*x，a=1/(1+2*pi*fc*Ts) */
        float Ts_ident = (FOC.Ts > 0.0f) ? FOC.Ts : T_10kHz;
        float k = 2.0f * 3.1415926f * fc_ident * Ts_ident;
        float a = 1.0f / (1.0f + k);

        if (a < 0.0f) a = 0.0f;
        if (a > 1.0f) a = 1.0f;

        Ud_ident_lpf.a = a;
        Uq_ident_lpf.a = a;
        Id_ident_lpf.a = a;
        Iq_ident_lpf.a = a;

        if (!ident_lpf_inited)
        {
          Ud_ident_lpf.y_last = Udi_ident;
          Uq_ident_lpf.y_last = Uqi_ident;
          Id_ident_lpf.y_last = Id_ident;
          Iq_ident_lpf.y_last = Iq_ident;
          ident_lpf_inited = 1u;
        }

        Udi_ident = LowPassFilter_Update(&Ud_ident_lpf, Udi_ident);
        Uqi_ident = LowPassFilter_Update(&Uq_ident_lpf, Uqi_ident);
        Id_ident = LowPassFilter_Update(&Id_ident_lpf, Id_ident);
        Iq_ident = LowPassFilter_Update(&Iq_ident_lpf, Iq_ident);
      }
      else
      {
        ident_lpf_inited = 0u;
      }

      /* 输出滤波后的辨识输入供上位机观测 */
      g_identUdFilt = Udi_ident;
      g_identUqFilt = Uqi_ident;
      g_identIdFilt = Id_ident;
      g_identIqFilt = Iq_ident;

      onlineMTPA_step_10k(Udi_ident, Uqi_ident, Id_ident, Iq_ident, we, Rs0,
                          0);  // Rs0=0.65Ω, flags=0

      /*FOC.Iq_ref = Speed_PID.output;  // Iq_ref = Speed_PID.output

      // FOC.Iq_ref = IQtest;
      // IQtest=IQtest+0.0001;
      // if(IQtest>IQtestMax) IQtest=0;

      if (FOC.Iq_ref > 0)
      {
        float x = FOC.Iq_ref;
        FOC.Id_ref =
            ((((0.000004986 * x - 0.0003467) * x + 0.009454) * x - 0.1289) * x + 1.286) * x -
            0.2316;  // MTPA
        // FOC.Id_ref = (((-0.00000387*x+0.000307)*x-0.00688)*x+0.303)*x+0.113; // FC-MTPA
      }
      else
      {
        float x = -FOC.Iq_ref;
        FOC.Id_ref =
            ((((0.000004986 * x - 0.0003467) * x + 0.009454) * x - 0.1289) * x + 1.286) * x -
            0.2316;  // MTPA
        // FOC.Id_ref = (((-0.00000387*x+0.000307)*x-0.00688)*x+0.303)*x+0.113; // FC-MTPA
      }*/

      //  float Iq_meas = FOC.Iq_ref; // 从传感器或速度环估计得到的 Iq 目标
      //  float Id_mtpa;
      //  float Iq_out;

      /*MTPA插值*/
      //  if (FOC.Iq_ref>=0)
      //     {
      //     MTPA_interp_by_Iq(mtpa_table, MTPA_TABLE_POINTS, Iq_meas, &Id_mtpa, &Iq_out);
      //     FOC.Id_ref = Id_mtpa;
      //     }
      //     else
      //     {
      //        MTPA_interp_by_Iq(mtpa_table, MTPA_TABLE_POINTS, -Iq_meas, &Id_mtpa, &Iq_out);
      //       FOC.Id_ref = Id_mtpa;
      //     }
      /*MTPA插值*/
      // FOC.Id_ref = IQtest;
      // if(FOC.Iq_ref>0)
      //   {FOC.Id_ref =FOC.Iq_ref;} // 限制 Id_ref <= 0
      // else
      //   {FOC.Id_ref = -FOC.Iq_ref;}

      PID_Controller(FOC.Id_ref, FOC.Id, &Id_PID);
      PID_Controller(FOC.Iq_ref, FOC.Iq, &Iq_PID);

      FOC.Ud_ref = Id_PID.output;
      FOC.Uq_ref = Iq_PID.output;

      break;
    }
    // !SECTION
    case EXIT:
    {
      STOP = 1;
      FOC.Id_ref = 0.0F;  // Id_ref = 0
      FOC.Iq_ref = 0.0F;  // Iq_ref = Speed_PID.output
      FOC.Ud_ref = 0.0F;
      FOC.Uq_ref = 0.0F;
      break;
    }
    // SECTION - Identify Mode
    case Identify:
    {
      ParkTransform(Clarke.Ialpha, Clarke.Ibeta, FOC.Theta, &FOC);

      // //SquareWaveGenerater(&VoltageInjector, &FOC);
      // HighFrequencySquareWaveGenerater(&VoltageInjector);

      // FOC.Ud_ref = VoltageInjector.Vd;
      // FOC.Uq_ref = VoltageInjector.Vq;
      // //FOC.Theta = VoltageInjector.Theta;  // 保持当前 Theta
      if (STOP == 0)
      {
        Experiment_Step(&Experiment, FOC.Id, FOC.Iq, &FOC.Ud_ref, &FOC.Uq_ref);
      }

      break;
    }
    // !SECTION
    default:
    {
      STOP = 1;
      FOC.Mode = IDLE;
      break;
    }
  }

  Electric_Power = Cal_Power(&FOC);

  InvParkTransform(FOC.Ud_ref, FOC.Uq_ref, FOC.Theta, &Inv_Park);
  SVPWM_Generate(Inv_Park.Ualpha, Inv_Park.Ubeta, FOC.inv_Udc, &FOC);
}
// !SECTION

static inline float Cal_Power(FOC_Parameter_t* foc)
{
  float power = 0.0F;
  power = 1.5 * (foc->Ud_ref * foc->Id + foc->Uq_ref * foc->Iq) - foc->Id * foc->Id * Motor.Rs -
          foc->Iq * foc->Iq * Motor.Rs;
  power = LowPassFilter_Update(&Power_Filter, power);
  Electric_Te =
      power <= 10.0F ? 0 : power * 60.0F / (M_2PI * foc->Speed + 0.0001F);  // 计算电磁转矩
  return power;
}

void FOC_UpdateMainFrequency(float f, float Ts, float PWM_ARR)
{
  FOC.f = f;
  FOC.Ts = Ts;
  FOC.PWM_ARR = PWM_ARR;
}
// SECTION - Parameter Init
void Parameter_Init(void)
{
  memset(&VF, 0, sizeof(VF_Parameter_t));
  // memset(&FOC, 0, sizeof(FOC_Parameter_t));
  memset(&Id_PID, 0, sizeof(PID_Controller_t));
  memset(&Iq_PID, 0, sizeof(PID_Controller_t));
  memset(&Speed_PID, 0, sizeof(PID_Controller_t));
  memset(&Inv_Park, 0, sizeof(InvPark_t));
  memset(&Speed_Ramp, 0, sizeof(RampGenerator_t));
  memset(&Motor, 0, sizeof(Motor_Parameter_t));
  memset(&MTPA, 0, sizeof(MTPA_Coefficients_t));

  STOP = 1;

  Motor.Rs = 0.65F;
  Motor.Ld = 0.18F;
  Motor.Lq = 0.12F;
  Motor.Flux = 0.1F;
  Motor.Pn = 2.0F;
  Motor.Position_Scale = 10000 - 1;
  Motor.Resolver_Pn = 1.0F;
  Motor.inv_MotorPn = 1.0F / 2.0F;  // Pn
  Motor.Position_Offset = 6838.0F;

  MTPA.A = 0.00061141F;
  MTPA.B = -0.014627F;
  MTPA.C = 0.34737F;
  MTPA.D = 0.068985F;

#ifdef Resolver_Position
  theta_factor = M_2PI / ((Motor.Position_Scale + 1) * Motor.Resolver_Pn);
#endif
#ifdef Encoder_Position
  theta_factor = M_2PI / (float)(Motor.Position_Scale + 1);
#endif
  Speed_PID.Kp = 0.023F;
  Speed_PID.Ki = 0.08F;
  Speed_PID.Kd = 0.0F;
  Speed_PID.MaxOutput = 0.7F * FOC.I_Max;  // Maximum Iq
  Speed_PID.MinOutput = -0.7F * FOC.I_Max;
  Speed_PID.IntegralLimit = 0.7F * FOC.I_Max;
  Speed_PID.previous_error = 0.0F;
  Speed_PID.integral = 0.0F;
  Speed_PID.output = 0.0F;
  Speed_PID.Ts = 10 * FOC.Ts;

  Speed_Ramp.slope = 250.0F;  // limit to 50 rpm/s
  Speed_Ramp.limit_min = -3000.0F;
  Speed_Ramp.limit_max = 3000.0F;
  Speed_Ramp.value = 0.0F;
  Speed_Ramp.target = 0.0F;
  Speed_Ramp.Ts = 10 * FOC.Ts;

  Id_PID.Kp = 73.8274273F;
  Id_PID.Ki = 408.40704496F;
  Id_PID.Kd = 0.0F;
  Id_PID.MaxOutput = 300.0F;  // Maximum Udc/sqrt(3)
  Id_PID.MinOutput = -300.0F;
  Id_PID.IntegralLimit = 300.0F;
  Id_PID.previous_error = 0.0F;
  Id_PID.integral = 0.0F;
  Id_PID.output = 0.0F;
  Id_PID.Ts = FOC.Ts;

  Iq_PID.Kp = 27.646015F;
  Iq_PID.Ki = 408.40704496F;
  Iq_PID.Kd = 0.0F;
  Iq_PID.MaxOutput = 350.0F;
  Iq_PID.MinOutput = -350.0F;
  Iq_PID.IntegralLimit = 350.0F;
  Iq_PID.previous_error = 0.0F;
  Iq_PID.integral = 0.0F;
  Iq_PID.output = 0.0F;
  Iq_PID.Ts = FOC.Ts;

  VF.Vref_Ud = 0.0F;
  VF.Vref_Uq = 0.0F;
  VF.Freq = 0.0F;
  VF.Theta = 0.0F;

  IF.Id_ref = 0.0F;
  IF.Iq_ref = 0.0F;
  IF.IF_Freq = 0.0F;
  IF.Theta = 0.0F;
  IF.Sensor_State = Disable;
}
// !SECTION
// SECTION - PID Controller
static inline void PID_Controller(float Ref, float Feedback, PID_Controller_t* PID_Controller)
{
  float difference = Ref - Feedback;
  float integral = PID_Controller->integral;
  float derivative = difference - PID_Controller->previous_error;

  if (STOP == 1)
  {
    difference = 0.0F;
    integral = 0.0F;
    derivative = 0.0F;
  }
  // Proportional term
  float P = PID_Controller->Kp * difference;
  // Integral term
  integral += PID_Controller->Ki * difference * PID_Controller->Ts;
  // Derivative term
  float D = PID_Controller->Kd * derivative;
  // Calculate output
  float output_value = P + integral + D;
  // Clamp output to limits
  if (output_value > PID_Controller->MaxOutput)
  {
    output_value = PID_Controller->MaxOutput;
  }
  else if (output_value < PID_Controller->MinOutput)
  {
    output_value = PID_Controller->MinOutput;
  }
  if (integral > PID_Controller->IntegralLimit)
  {
    integral = PID_Controller->IntegralLimit;
  }
  else if (integral < -PID_Controller->IntegralLimit)
  {
    integral = -PID_Controller->IntegralLimit;
  }
  // Update integral and previous error for next iteration
  PID_Controller->integral = integral;
  PID_Controller->previous_error = difference;
  // Return the output value
  PID_Controller->output = output_value;
}
// !SECTION

static inline float wrap_theta_2pi(float theta)
{
  if (theta >= M_2PI)
  {
    theta -= M_2PI;
  }
  else if (theta < 0.0F)
  {
    theta += M_2PI;
  }
  return theta;
}

LowPassFilter_t Speed_Filter = {.a = 0.98442F, .y_last = 0.0F};

static inline float LowPassFilter_Update(LowPassFilter_t* filter, float x)
{
  float y = filter->a * filter->y_last + (1.0F - filter->a) * x;
  filter->y_last = y;
  return y;
}

// SECTION - Theta Process
static inline void Theta_Process(void)
{
  int32_t pos_delta = (int32_t)((float)FOC.Position - Motor.Position_Offset);
  if (pos_delta < 0)
  {
    pos_delta += (Motor.Position_Scale + 1);
  }

  theta_mech = (float)pos_delta * theta_factor;

  theta_elec = theta_mech * Motor.Pn;
  FOC.Theta = wrap_theta_2pi(theta_elec);

  static float last_theta = 0.0F;
  float delta_theta = theta_mech - last_theta;

  if (delta_theta > M_PI)
  {
    delta_theta -= M_2PI;
  }
  else if (delta_theta < -M_PI)
  {
    delta_theta += M_2PI;
  }

  last_theta = theta_mech;

  float Speed = delta_theta * FOC.f * 60.0F / M_2PI;
  FOC.Speed = LowPassFilter_Update(&Speed_Filter, Speed);
}
// !SECTION

// SECTION - Ramp Generator
static inline float RampGenerator(RampGenerator_t* ramp)
{
  if (STOP == 1)
  {
    ramp->value = 0.0F;
    ramp->target = 0.0F;
  }

  float delta = ramp->target - ramp->value;
  float step = ramp->slope * ramp->Ts;

  if (delta > step)
  {
    ramp->value += step;
  }
  else if (delta < -step)
  {
    ramp->value -= step;
  }
  else
  {
    ramp->value = ramp->target;
  }

  if (ramp->value > ramp->limit_max)
  {
    ramp->value = ramp->limit_max;
  }
  if (ramp->value < ramp->limit_min)
  {
    ramp->value = ramp->limit_min;
  }

  return ramp->value;
}
// !SECTION

static inline void ClarkeTransform(float_t Ia, float_t Ib, float_t Ic, Clarke_t* out)
{
#if (defined(TWO_PHASE_CURRENT_SENSING))
  out->Ialpha = Ia;
  out->Ibeta = 0.57735026919F * (Ia + 2.0F * Ib);  // 0.57735026919F 1/√3
#elif (defined(THREE_PHASE_CURRENT_SENSING))
  out->Ialpha = 0.66666666667F * Ia - 0.33333333333F * (Ib + Ic);
  out->Ibeta = 0.57735026919F * (Ib - Ic);  // 0.57735026919F 1/√3
#endif
}

static inline void ParkTransform(float_t Ialpha, float_t Ibeta, float_t theta, FOC_Parameter_t* out)
{
  float cos_theta = COS(theta);
  float sin_theta = SIN(theta);
  out->Id = Ialpha * cos_theta + Ibeta * sin_theta;
  out->Iq = -Ialpha * sin_theta + Ibeta * cos_theta;
}

static inline void InvParkTransform(float_t Ud, float_t Uq, float_t theta, InvPark_t* out)
{
  float cos_theta = COS(theta);
  float sin_theta = SIN(theta);
  out->Ualpha = Ud * cos_theta - Uq * sin_theta;
  out->Ubeta = Ud * sin_theta + Uq * cos_theta;
}
static inline float sat(float x)
{
  if (x > 1.0F) return 1.0F;
  if (x < -1.0F) return -1.0F;
  return x;
}
static inline float Get_Theta(float Freq, float Theta)
{
  // 电角度递推：θ += ω·Ts，ω = 2π·f
  Theta += M_2PI * Freq * FOC.Ts;
  if (Theta > M_2PI)
  {
    Theta -= M_2PI;
  }
  else if (Theta < 0.0F)
  {
    Theta += M_2PI;
  }
  return Theta;
}

static inline void SVPWM_Generate(float Ualpha, float Ubeta, float inv_Vdc, FOC_Parameter_t* foc)
{
  uint8_t sector = 0;
  float Vref1 = Ubeta;
  float Vref2 = (+SQRT3 * Ualpha - Ubeta) * 0.5F;
  float Vref3 = (-SQRT3 * Ualpha - Ubeta) * 0.5F;

  // 判断扇区（1~6）
  if (Vref1 > 0) sector += 1;
  if (Vref2 > 0) sector += 2;
  if (Vref3 > 0) sector += 4;

  // Clarke to T1/T2 projection
  float X = SQRT3 * Ubeta * inv_Vdc;
  float Y = (+1.5F * Ualpha + SQRT3_2 * Ubeta) * inv_Vdc;
  float Z = (-1.5F * Ualpha + SQRT3_2 * Ubeta) * inv_Vdc;

  float T1 = 0.0F, T2 = 0.0F;

  switch (sector)
  {
    case 1:
      T1 = Z;
      T2 = Y;
      break;
    case 2:
      T1 = Y;
      T2 = -X;
      break;
    case 3:
      T1 = -Z;
      T2 = X;
      break;
    case 4:
      T1 = -X;
      T2 = Z;
      break;
    case 5:
      T1 = X;
      T2 = -Y;
      break;
    case 6:
      T1 = -Y;
      T2 = -Z;
      break;
    default:
      T1 = 0.0F;
      T2 = 0.0F;
      break;
  }

  // 过调制处理
  float T_sum = T1 + T2;
  if (T_sum > 1.0F)
  {
    T1 /= T_sum;
    T2 /= T_sum;
  }

  // 中心对称调制时间计算
  float T0 = (1.0F - T1 - T2) * 0.5F;
  float Ta = T0;
  float Tb = T0 + T1;
  float Tc = Tb + T2;

  float Tcm1 = 0.0F;
  float Tcm2 = 0.0F;
  float Tcm3 = 0.0F;

  // 扇区映射到ABC换相点
  switch (sector)
  {
    case 1:
      Tcm1 = Tb;
      Tcm2 = Ta;
      Tcm3 = Tc;
      break;
    case 2:
      Tcm1 = Ta;
      Tcm2 = Tc;
      Tcm3 = Tb;
      break;
    case 3:
      Tcm1 = Ta;
      Tcm2 = Tb;
      Tcm3 = Tc;
      break;
    case 4:
      Tcm1 = Tc;
      Tcm2 = Tb;
      Tcm3 = Ta;
      break;
    case 5:
      Tcm1 = Tc;
      Tcm2 = Ta;
      Tcm3 = Tb;
      break;
    case 6:
      Tcm1 = Tb;
      Tcm2 = Tc;
      Tcm3 = Ta;
      break;
    default:
      Tcm1 = 0.5F;
      Tcm2 = 0.5F;
      Tcm3 = 0.5F;
      break;
  }

  foc->Tcm1 = Tcm1;
  foc->Tcm2 = Tcm2;
  foc->Tcm3 = Tcm3;
}
