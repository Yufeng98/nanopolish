set -v
vtune -collect-with runsa -start-paused -knob event-config=INST_RETIRED.ANY,CPU_CLK_UNHALTED.THREAD,MEM_UOPS_RETIRED.ALL_LOADS:sa=100000,MEM_LOAD_UOPS_RETIRED.L1_HIT:sa=100000,MEM_LOAD_UOPS_RETIRED.L2_HIT:sa=100,MEM_LOAD_UOPS_RETIRED.L3_HIT:sa=100,MEM_LOAD_UOPS_RETIRED.L1_MISS:sa=100,MEM_LOAD_UOPS_RETIRED.L2_MISS:sa=100,MEM_LOAD_UOPS_RETIRED.L3_MISS:sa=100,OFFCORE_RESPONSE:request=ALL_DATA_RD:response=LLC_MISS.REMOTE_DRAM:sa=1000,OFFCORE_RESPONSE:request=ALL_DATA_RD:response=LLC_MISS.LOCAL_DRAM:sa=1000,OFFCORE_RESPONSE:request=ALL_READS:response=LLC_MISS.LOCAL_DRAM:sa=1000,OFFCORE_RESPONSE:request=ALL_READS:response=LLC_MISS.REMOTE_DRAM:sa=1000,BR_INST_RETIRED.ALL_BRANCHES,BR_MISP_RETIRED.ALL_BRANCHES,DTLB_LOAD_MISSES.MISS_CAUSES_A_WALK ./$1