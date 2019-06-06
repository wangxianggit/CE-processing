/******************************/
/*      Setup Definition      */
/******************************/
#include <stdint.h>
#include "CAENDigitizerType.h"

// module number
#define MAXNB_V1724    27
#define MAXNB_V1730    0
#define MAXNB          (MAXNB_V1724 + MAXNB_V1730)

// define number of samples to be recorded
#define NSAMPLE_V1724   1500
#define NSAMPLE_V1730   1600

// maximum channel number
#define MAXNCH_V1724   8
#define MAXNCH_V1730   16
#define MAXNCH         (MAXNB_V1730? MAXNCH_V1730 : MAXNCH_V1724)

// sampling rate
#define RATE_V1724     (100E6)     // 100 MHz
#define CLOCK_V1724    (10E-9)     // 10 ns
#define RATE_V1730     (500E6)     // 500 MHz
#define CLOCK_V1730    (2E-9)      // 2 ns

// V830 Scaler
#define USE_V830_SCALER
#define V830_BASE      0xC0000000
#define V830_CH_NUM    7

// Wave histo fill interval in roody
#define WAVE_FILL_INTERVAL    5    // second

// Read all registers when initialization
#define READ_ALL_REGISTERS

// display waveform or not in roody
#define DISPLAY_WAVEFORM

// Sampling CLK propagated to TRG_OUT (used for clock synchronization)
//#define CLK_TRG_OUT

// Disable the self-trigger for all channels, used to find the delay time of a run for each boards or external triggered DAQ
//#define SELF_TRIGGER_DISABLE

// display all debug message
#define DEBUG_MSG_ONLINE      1    // 1 for display message, 0 for not
#define DEBUG_MSG_OFFLINE     1    // 1 for display message, 0 for not

// define variables for the interrupt configuration
#define VME_BUS_NUM           1    // the number of VME_Bus_Link
#define VME_INTERRUPT_LEVEL   1    // must be 1 for direct connection through CONET
#define MAX_NUM_EVENTS_BLT    0    // if set to be 0, disable the interrupt configuration
#define INTERRUPT_MODE        CAEN_DGTZ_IRQ_MODE_ROAK
#define INTERRUPT_TIMEOUT     10   // unit: ms

// define the structure of digitizer papameters
typedef struct {
	CAEN_DGTZ_ConnectionType LinkType;
	int LinkNum;
	int ConetNode;
	uint32_t VMEBaseAddress;
	int EventAggr;
	uint32_t ChannelMask;
	uint16_t DynamicRange_V1730[MAXNCH_V1730];
	uint32_t RecordLength_V1724;
	uint32_t RecordLength_V1730[MAXNCH_V1730];
	CAEN_DGTZ_DPP_AcqMode_t DPPAcqMode;
	CAEN_DGTZ_DPP_SaveParam_t SaveParam;
	CAEN_DGTZ_AcqMode_t AcqMode;
	CAEN_DGTZ_IOLevel_t IOlev;
	CAEN_DGTZ_TriggerMode_t SWTrgMode;
	CAEN_DGTZ_TriggerMode_t ExtTrgMode;
	CAEN_DGTZ_TriggerMode_t SelfTrgMode;
	CAEN_DGTZ_RunSyncMode_t RunSyncMode;
	int VirtualProbe1;
	int VirtualProbe2;
	int DigitalProbe;
	CAEN_DGTZ_PulsePolarity_t PulsePolarity[MAXNCH_V1730];
	uint32_t DPPPreTriggerSize[MAXNCH_V1730];
	float ChannelDCOffset[MAXNCH_V1730];
	uint32_t Run_Start_Stop_Delay;
	CAEN_DGTZ_AnalogMonitorOutputMode_t AnalogMonOutput;
} DigitizerParams_t;


