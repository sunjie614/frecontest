/*!
    \file    gd32f30x_it.c
    \brief   interrupt service routines

   \version 2024-12-20, V3.0.1, firmware for GD32F30x
*/

/*
    Copyright (c) 2024, GigaDevice Semiconductor Inc.

    Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
       list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.
    3. Neither the name of the copyright holder nor the names of its contributors
       may be used to endorse or promote products derived from this software without
       specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/

#include "gd32f30x_it.h"
#include "can.h"
#include "control_profile.h"
#include "foc.h"
#include "gd32f30x.h"
#include "hardware_interface.h"
#include "justfloat.h"
#include "onlineMTPA.h"
#include "systick.h"

extern volatile uint16_t STOP;
float Id_max = 0.0F;
float Theta_max = 0.0F;

/*!
    \brief      this function handles NMI exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void NMI_Handler(void) {}

/*!
    \brief      this function handles HardFault exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void HardFault_Handler(void)
{
  /* if Hard Fault exception occurs, go to infinite loop */
  while (1)
  {
  }
}

/*!
    \brief      this function handles MemManage exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void MemManage_Handler(void)
{
  /* if Memory Manage exception occurs, go to infinite loop */
  while (1)
  {
  }
}

/*!
    \brief      this function handles BusFault exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void BusFault_Handler(void)
{
  /* if Bus Fault exception occurs, go to infinite loop */
  while (1)
  {
  }
}

/*!
    \brief      this function handles UsageFault exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void UsageFault_Handler(void)
{
  /* if Usage Fault exception occurs, go to infinite loop */
  while (1)
  {
  }
}

/*!
    \brief      this function handles DebugMon exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void DebugMon_Handler(void) {}

/*!
    \brief      this function handles SysTick exception
    \param[in]  none
    \param[out] none
    \retval     none
*/
void SysTick_Handler(void)
{
  delay_decrement();
  systick_ms++;
}

void DMA0_Channel3_IRQHandler(void)
{
  if (dma_interrupt_flag_get(DMA0, DMA_CH3, DMA_INT_FLAG_FTF))
  {
    dma_interrupt_flag_clear(DMA0, DMA_CH3, DMA_INT_FLAG_FTF);
    Peripheral_SCISendCallback();
  }
}

void USBD_LP_CAN0_RX0_IRQHandler(void)
{
  can_receive_message_struct rx_msg;
  can_message_receive(CAN0, CAN_FIFO0, &rx_msg);
  CAN_Buffer_Put(&rx_msg);
}
uint32_t systick_cnt = 0;
float IQtest = 0;
float IQtestMax = 20.5F;
uint32_t systick_cnt1 = 0;

static uint32_t Timer1_DeltaTicks(uint32_t start, uint32_t end)
{
  uint32_t period = (TIMER_CAR(TIMER1) & TIMER_CAR_CARL) + 1u;
  start &= TIMER_CNT_CNT;
  end &= TIMER_CNT_CNT;

  if (end >= start)
  {
    return end - start;
  }

  return (period - start) + end;
}

static uint32_t Timer1_TicksToCycles(uint32_t ticks)
{
  uint32_t tick_cycles = (TIMER_PSC(TIMER1) & TIMER_PSC_PSC) + 1u;
  g_prof_tmr1_tick_cycles = tick_cycles;
  return ticks * tick_cycles;
}

void ADC0_1_IRQHandler(void)
{
  uint32_t dwt_isr_start = DWT->CYCCNT;
  uint32_t tmr1_isr_start = TIMER_CNT(TIMER1);
  uint32_t dwt_comm_cycles = 0u;
  uint32_t tmr1_comm_ticks = 0u;
  uint8_t control_isr_handled = 0u;
  uint8_t method = ControlProfile_GetActiveMethod();
  g_prof_speed_cycle_hit = 0u;
  g_prof_ls_solve_hit = 0u;

  uint32_t cnt_start = tmr1_isr_start;
  if (adc_interrupt_flag_get(ADC0, ADC_INT_FLAG_EOIC))
  {
    control_isr_handled = 1u;
    adc_interrupt_flag_clear(ADC0, ADC_INT_FLAG_EOIC);

    Peripheral_UpdateCurrent();
    Peripheral_GateState();
    Peripheral_UpdateUdc();
    Peripheral_UpdatePosition();

    switch (FOC.Mode)
    {
      case INIT:
      {
        Peripheral_InitProtectParameter();
        Peripheral_GetSystemFrequency();
        Peripheral_CalibrateADC();
        if (FOC.Udc > 200.0F)
        {
          Peripheral_EnableHardwareProtect();
        }
        Protect.Flag = No_Protect;

        break;
      }
      case IDLE:
      case VF_MODE:
      case IF_MODE:
      case Speed:
      {
        break;
      }
      case Identify:
      {
        if (FOC.Id > Id_max)
        {
          Id_max = FOC.Id;
          Theta_max = FOC.Theta;
        }

        // float DMA_Buffer[5];
        // DMA_Buffer[0] = VoltageInjector.Vd;
        // DMA_Buffer[1] = VoltageInjector.Vq;
        // DMA_Buffer[2] = FOC.Id;
        // DMA_Buffer[3] = FOC.Iq;
        // DMA_Buffer[4] = (float)VoltageInjector.Count;
        // justfloat(DMA_Buffer, 5);
        break;
      }
      case EXIT:
      {
        Peripheral_DisableHardwareProtect();
        break;
      }
      default:
        break;
    }

    FOC_Main();
    onlineMTPA_conf_t conf = onlineMTPA_get_confidence();
    float psid_model = 0.0f;
    float psiq_model = 0.0f;
    onlineMTPA_flux_from_i(FOC.Id, FOC.Iq, &psid_model, &psiq_model);
    float DMA_Buffer[16];
    const uint8_t dma_float_count = (uint8_t)(sizeof(DMA_Buffer) / sizeof(DMA_Buffer[0]));
    DMA_Buffer[0] = FOC.Ia;
    DMA_Buffer[1] = gamma_deg;
    // DMA_Buffer[1] = FOC.Speed;
    //  DMA_Buffer[2] = FOC.Uq_ref;
    //  DMA_Buffer[3] = FOC.Id;
    //  DMA_Buffer[4] = FOC.Iq;
    //  DMA_Buffer[5] = FOC.Ud_ref;
    // DMA_Buffer[2] = g_ld;
    // DMA_Buffer[3] = g_lq;
    DMA_Buffer[2] = g_identIdFilt;
    DMA_Buffer[3] = g_identIqFilt;
    // DMA_Buffer[4] = g_ldd;
    // DMA_Buffer[5] = g_lqq;
    DMA_Buffer[4] = g_mtpaKNow;
    DMA_Buffer[5] = g_mtpaMNow;
    DMA_Buffer[6] = g_lc;
    DMA_Buffer[7] = FOC.Id;
    DMA_Buffer[8] = FOC.Iq;
    DMA_Buffer[9] = FOC.Speed;
    DMA_Buffer[10] = g_identUdFilt;
    DMA_Buffer[11] = g_torqueGradFiltNow;
    DMA_Buffer[12] = g_identUqFilt;
    DMA_Buffer[13] = psid_model;
    DMA_Buffer[14] = psiq_model;
    DMA_Buffer[15] = conf.Cp;
    uint32_t dwt_comm_start = DWT->CYCCNT;
    uint32_t tmr1_comm_start = TIMER_CNT(TIMER1);
    justfloat(DMA_Buffer, dma_float_count);
    dwt_comm_cycles += DWT->CYCCNT - dwt_comm_start;
    tmr1_comm_ticks += Timer1_DeltaTicks(tmr1_comm_start, TIMER_CNT(TIMER1));

    Peripheral_SetPWMChangePoint();
  }
  uint32_t cnt_end = TIMER_CNT(TIMER1);

  if (cnt_end > cnt_start)
  {
    systick_cnt1 = cnt_end - cnt_start;
    if (systick_cnt1 > systick_cnt)
    {
      systick_cnt = systick_cnt1;
    }
  }

  if (control_isr_handled)
  {
    uint32_t dwt_end = DWT->CYCCNT;
    uint32_t dwt_raw_cycles = dwt_end - dwt_isr_start;
    uint32_t tmr1_raw_cycles = Timer1_TicksToCycles(Timer1_DeltaTicks(tmr1_isr_start, cnt_end));
    uint32_t tmr1_comm_cycles = Timer1_TicksToCycles(tmr1_comm_ticks);
    uint8_t dwt_ready = (((CoreDebug->DEMCR & CoreDebug_DEMCR_TRCENA_Msk) != 0u) &&
                         ((DWT->CTRL & DWT_CTRL_CYCCNTENA_Msk) != 0u) &&
                         (dwt_raw_cycles != 0u)) ?
                            1u :
                            0u;
    uint32_t raw_cycles = dwt_ready ? dwt_raw_cycles : tmr1_raw_cycles;
    uint32_t comm_cycles = dwt_ready ? dwt_comm_cycles : tmr1_comm_cycles;
    uint32_t control_cycles = raw_cycles;
    if (control_cycles >= comm_cycles)
    {
      control_cycles -= comm_cycles;
    }

    g_prof_timer_src = dwt_ready ? 0u : 1u;
    g_prof_dwt_ok = dwt_ready;
    g_prof_dwt_demcr = CoreDebug->DEMCR;
    g_prof_dwt_ctrl = DWT->CTRL;
    g_prof_dwt_cyccnt = dwt_end;
    g_prof_last_isr_cycles = raw_cycles;
    g_prof_last_ctrl_cycles = control_cycles;
    g_prof_last_comm_cycles = comm_cycles;

    ControlProfile_Record(method, PROFILE_METRIC_ALL_ISR, control_cycles);
    if (g_prof_speed_cycle_hit)
    {
      ControlProfile_Record(method, PROFILE_METRIC_SPEED_CYCLE_ISR, control_cycles);
    }
    if (g_prof_ls_solve_hit)
    {
      ControlProfile_Record(method, PROFILE_METRIC_LS_SOLVE_ISR, control_cycles);
    }
  }
}

void EXTI4_IRQHandler(void)
{
  if (RESET != exti_interrupt_flag_get(EXTI_4))
  {
    TIMER_SWEVG(TIMER0) |= TIMER_SWEVG_BRKG;
    STOP = 1;
    exti_interrupt_flag_clear(EXTI_4);
  }
}

extern Protect_Flags Protect_Flag;
extern bool Software_BRK;

void TIMER0_BRK_IRQHandler(void)
{
  if (timer_interrupt_flag_get(TIMER0, TIMER_INT_FLAG_BRK))
  {
    // 清除 Break 中断标志
    timer_interrupt_flag_clear(TIMER0, TIMER_INT_FLAG_BRK);
    STOP = 1;
    if (Software_BRK == true)
    {
      Protect.Flag |= Hardware_Fault;
      timer_interrupt_disable(TIMER0, TIMER_INT_BRK);  // 禁用BRK中断
      timer_primary_output_config(TIMER0, DISABLE);
    }
  }
}

/* */
void TIMER3_IRQHandler(void)
{
  if (timer_interrupt_flag_get(TIMER3, TIMER_INT_FLAG_CH2))
  {
    // 清除 CH2 中断标志
    timer_interrupt_flag_clear(TIMER3, TIMER_INT_FLAG_CH2);

    TIMER_CNT(TIMER3) = 0;
  }
}
