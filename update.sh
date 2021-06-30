scp src/nanopolish_call_methylation.cpp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src
scp src/nanopolish_raw_loader.cpp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src
scp src/nanopolish_squiggle_read.h yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src
scp src/nanopolish_squiggle_read.cpp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src
scp src/hmm/nanopolish_emissions.h yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/hmm
scp src/common/nanopolish_common.h yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/common
# scp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/nanopolish_squiggle_read.cpp src/
# scp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/nanopolish_raw_loader.cpp src/
# scp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/nanopolish_squiggle_read.h src/
# scp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/hmm/nanopolish_emissions.h src/hmm/
# scp yufenggu@mbit7.eecs.umich.edu:~/nanopolish/src/common/nanopolish_common.h src/commom/

# list=("db0badd6-b0c9-4e7e-b6ca-bc0deb8741d8" "d763e7f0-c28f-413d-aa80-345eeb054f62" "f77701dc-077a-47e8-827c-772550a48464" "cf1561b0-b570-41b6-bf54-c116a7fb3604" "bf5f4c8d-e1fb-475a-87bc-0c9196e8e469" "93c765dc-31df-4e52-89ae-55c00d636b70" "b6244afc-d75f-443e-b156-b2c4d1a6ed1d" "15335b85-1f8f-4b30-aac6-2291a032f9cb" "3cfd1cea-f662-4ef2-b4c8-12cfdf25dbe8" "5073afa1-3236-49bd-9ad7-5077f007e533" "4d929162-76ea-4bdd-84ff-002bd38e85a7" "371ea24d-1d18-4412-ac39-57988be2b08c" "383cc4b4-84e0-4cea-a9d7-b806793f088f" "84a2c397-6399-47be-9c8d-220ed27536d5" "fd7e7616-edac-4cbb-841a-5a9bc3dc1f48" "82cabb64-9154-4d3d-a94f-6b7d0877612f" "488f4ebd-0331-456b-b8d1-3f1d6d96cc57" "22b97cc8-37e6-4f0e-aa9d-2209843fee3b" "5eea4c5f-02ab-48fb-837c-56f3ee954821" "073a8fd8-144e-4e6a-8b7c-98a69b439311" "565b72a0-ebae-422d-8092-46720fc18e84" "79a90e09-2b7b-41b6-ac3e-6c2e949e0901" "3b7d2fbe-09af-4696-989f-2acef48af72d" "61bef566-da66-41ad-822b-7cbeec37e1c2" "fd717478-1ba4-4f1b-b473-ff4733eec652")

# for name in "${list[@]}"; do
# #   diff "max_fp_${name}.txt" "max_approximation_${name}.txt" > "diff_max_${name}"
#     scp yufenggu@mbit7.eecs.umich.edu:/z/scratch7/yufenggu/methylation_example/diff_max_${name} ../punnet/long-reads/nanopolish/
# done