import sys

data_type = sys.argv[1]

# read = sys.argv[2]
# begin = int(sys.argv[3])
# length = int(sys.argv[4])

# f_w = open("{}_{}.txt".format(read, data_type), "w")
# f = open("trace_{}.txt".format(data_type), "r")
# f_lines = f.readlines()
# for line in f_lines[begin-1:begin+1+length]:
#     f_w.write(line)
# f.close()
# f_w.close()

    

# f_w = open("trace_{}_read.txt".format(data_type), "w")

# f_read_num = 0
# f = open("trace_{}.txt".format(data_type), "r")
# f_lines = f.readlines()
# f_read = []
# for line in f_lines:
#     if len(line) > 35: 
#         f_w.write(line)
#         f_read_num += 1
#         f_read.append(line[:-1])
# f.close()
# f_w.close()

# print("{} reads: {}".format(data_type, f_read_num))

# fp_w = open("trace_fp_read.txt", "w")

# fp_read_num = 0
# f_fp = open("trace_fp.txt", "r")
# fp_lines = f_fp.readlines()
# fp_read = []
# for line in fp_lines:
#     if len(line) > 35:
#         fp_w.write(line)
#         fp_read_num += 1
#         if line[:-1] in f_read:
#             fp_read.append(line[:-1])
#         # else:
#             # print(line[:-1])
# fp_w.close()
# f_fp.close()
# print("fp reads: {}".format(fp_read_num))

# print(len(f_read), len(fp_read)) 

# for read in f_read:
#     if read not in fp_read:
#         f_read.remove(read)
#     else:
#         print(read)

# for read in fp_read:
#     if read not in f_read:
#         fp_read.remove(read)
#     else:
#         print(read)

# print(len(f_read), len(fp_read)) 

def check_error_for_all_reads(data_type):

    f = open("methylation_calls_fix_0_64_{}.tsv".format(data_type), "r")
    # f = open("methylation_calls.tsv", "r")
    fix_lines = f.readlines()
    fix_lines = fix_lines[1:]
    fixed = []
    for line in fix_lines:
        sep = list(line.split())
        fixed.append((sep[1] + sep[2] + sep[3] + sep[4], sep[5]))

    f = open("methylation_calls_fp_0_64.tsv", "r")
    # f = open("methylation_calls_fix_0_64_32_32.tsv", "r")
    fp_lines = f.readlines()
    fp_lines = fp_lines[1:]
    fp = []
    for line in fp_lines:
        sep = list(line.split())
        fp.append((sep[1] + sep[2] + sep[3] + sep[4], sep[5]))

    dic = {}
    dic_read = {}
    # for i in range(len(fixed)):
    for i in range(10000):
        if fixed[i][0] in dic.keys(): dic[fixed[i][0]] += [("fix", fixed[i][1])]
        else: dic[fixed[i][0]] = [("fix", fixed[i][1])]
        if (i%10000 == 0): print("fix", i)

    # for i in range(len(fp)):
    for i in range(10000):
        if fp[i][0] in dic.keys(): dic[fp[i][0]] += [("fp", fp[i][1])]
        if (i%10000 == 0): print("fp", i)

    # print(len(name_fix), len(name_fp))
    # print(len(fixed_lst), len(fp_lst))

    diff_ratio = 0
    diff_ratio_0 = 0
    diff_ratio_0001 = 0
    diff_ratio_0002 = 0
    diff_ratio_0004 = 0
    diff_ratio_001 = 0
    diff_ratio_002 = 0
    diff_ratio_004 = 0
    diff_ratio_01 = 0
    diff_ratio_02 = 0
    diff_ratio_04 = 0
    diff_ratio_1 = 0
    diff_else = 0

    for key in dic.keys():
        if (len(dic[key]) == 2 and dic[key][0][0] == "fix" and dic[key][1][0] == "fp"):
            diff = abs(float(dic[key][0][1])-float(dic[key][1][1]))
            print(key[-36:])
            if key[-36:] in dic_read.keys(): dic_read[key[-36:]] += [diff]
            else: dic_read[key[-36:]] = [diff]
            diff_ratio += 1
            if diff == 0: diff_ratio_0 += 1
            elif diff < 0.001: diff_ratio_0001 += 1
            elif diff < 0.002: diff_ratio_0002 += 1
            elif diff < 0.004: diff_ratio_0004 += 1
            elif diff < 0.01: diff_ratio_001 += 1
            elif diff < 0.02: diff_ratio_002 += 1
            elif diff < 0.04: diff_ratio_004 += 1
            elif diff < 0.1: diff_ratio_01 += 1
            elif diff < 0.2: diff_ratio_02 += 1
            elif diff < 0.4: diff_ratio_04 += 1
            elif diff < 1: diff_ratio_1 += 1
                # print(float(dic[key][0][1]), float(dic[key][1][1]))
            else: diff_else += 1

    for key in dic_read.keys:
        print(key, len(dic_read))
    print("Ratio", diff_ratio, diff_ratio_0, diff_ratio_0001, diff_ratio_0002, diff_ratio_0004, diff_ratio_001, diff_ratio_002, diff_ratio_004, diff_ratio_01, diff_ratio_02, diff_ratio_04, diff_ratio_1, diff_else)


def check_error_for_all_read_parts(data_type):

    f = open("methylation_calls_fix_0_64_{}.tsv".format(data_type), "r")
    # f = open("methylation_calls.tsv", "r")
    fix_lines = f.readlines()
    fix_lines = fix_lines[1:]
    fixed = []
    for line in fix_lines:
        sep = list(line.split())
        fixed.append((sep[1] + sep[2] + sep[3] + sep[4], sep[5], sep[6], sep[7]))

    f = open("methylation_calls_fp_0_64.tsv", "r")
    # f = open("methylation_calls_fix_0_64_32_32.tsv", "r")
    fp_lines = f.readlines()
    fp_lines = fp_lines[1:]
    fp = []
    for line in fp_lines:
        sep = list(line.split())
        fp.append((sep[1] + sep[2] + sep[3] + sep[4], sep[5], sep[6], sep[7]))

    dic = {}
    for i in range(len(fixed)):
    # for i in range(20000):
        if fixed[i][0] in dic.keys(): dic[fixed[i][0]] += [("fix", fixed[i][1], fixed[i][2], fixed[i][3])]
        else: dic[fixed[i][0]] = [("fix", fixed[i][1], fixed[i][2], fixed[i][3])]
        if (i%10000 == 0): print("fix", i)

    for i in range(len(fp)):
    # for i in range(20000):
        if fp[i][0] in dic.keys(): dic[fp[i][0]] += [("fp", fp[i][1], fp[i][2], fp[i][3])]
        if (i%10000 == 0): print("fp", i)

    # print(len(name_fix), len(name_fp))
    # print(len(fixed_lst), len(fp_lst))

    diff_ratio = []
    diff_ratio_0 = 0
    diff_ratio_0001 = 0
    diff_ratio_0002 = 0
    diff_ratio_0004 = 0
    diff_ratio_001 = 0
    diff_ratio_002 = 0
    diff_ratio_004 = 0
    diff_ratio_01 = 0
    diff_ratio_02 = 0
    diff_ratio_04 = 0
    diff_ratio_1 = 0
    diff_else = 0

    for key in dic.keys():
        if len(dic[key]) == 2:
            if (dic[key][0][0] == "fix" and dic[key][1][0] == "fp"):
                diff = abs(float(dic[key][0][1])-float(dic[key][1][1]))
                diff_ratio.append(diff)
                if diff == 0: diff_ratio_0 += 1
                elif diff < 0.001: diff_ratio_0001 += 1
                elif diff < 0.002: diff_ratio_0002 += 1
                elif diff < 0.004: diff_ratio_0004 += 1
                elif diff < 0.01: diff_ratio_001 += 1
                elif diff < 0.02: diff_ratio_002 += 1
                elif diff < 0.04: diff_ratio_004 += 1
                elif diff < 0.1: diff_ratio_01 += 1
                elif diff < 0.2: diff_ratio_02 += 1
                elif diff < 0.4: diff_ratio_04 += 1
                elif diff < 1: diff_ratio_1 += 1
                    # print(float(dic[key][0][1]), float(dic[key][1][1]))
                else: diff_else += 1

    print("Ratio", len(diff_ratio), diff_ratio_0, diff_ratio_0001, diff_ratio_0002, diff_ratio_0004, diff_ratio_001, diff_ratio_002, diff_ratio_004, diff_ratio_01, diff_ratio_02, diff_ratio_04, diff_ratio_1, diff_else)

    diff_ratio = []
    diff_ratio_0 = 0
    diff_ratio_0001 = 0
    diff_ratio_0002 = 0
    diff_ratio_0004 = 0
    diff_ratio_001 = 0
    diff_ratio_002 = 0
    diff_ratio_004 = 0
    diff_ratio_01 = 0
    diff_ratio_02 = 0
    diff_ratio_04 = 0
    diff_ratio_1 = 0
    diff_else = 0

    for key in dic.keys():
        if len(dic[key]) == 2:
            if (dic[key][0][0] == "fix" and dic[key][1][0] == "fp"):
                diff = abs(float(dic[key][0][2])-float(dic[key][1][2]))
                diff_ratio.append(diff)
                if diff == 0: diff_ratio_0 += 1
                elif diff < 0.001: diff_ratio_0001 += 1
                elif diff < 0.002: diff_ratio_0002 += 1
                elif diff < 0.004: diff_ratio_0004 += 1
                elif diff < 0.01: diff_ratio_001 += 1
                elif diff < 0.02: diff_ratio_002 += 1
                elif diff < 0.04: diff_ratio_004 += 1
                elif diff < 0.1: diff_ratio_01 += 1
                elif diff < 0.2: diff_ratio_02 += 1
                elif diff < 0.4: diff_ratio_04 += 1
                elif diff < 1: diff_ratio_1 += 1
                    # print(float(dic[key][0][2]), float(dic[key][1][2]))
                else: diff_else += 1

    print("Methylated Score", len(diff_ratio), diff_ratio_0, diff_ratio_0001, diff_ratio_0002, diff_ratio_0004, diff_ratio_001, diff_ratio_002, diff_ratio_004, diff_ratio_01, diff_ratio_02, diff_ratio_04, diff_ratio_1, diff_else)

    diff_ratio = []
    diff_ratio_0 = 0
    diff_ratio_0001 = 0
    diff_ratio_0002 = 0
    diff_ratio_0004 = 0
    diff_ratio_001 = 0
    diff_ratio_002 = 0
    diff_ratio_004 = 0
    diff_ratio_01 = 0
    diff_ratio_02 = 0
    diff_ratio_04 = 0
    diff_ratio_1 = 0
    diff_else = 0

    for key in dic.keys():
        if len(dic[key]) == 2:
            if (dic[key][0][0] == "fix" and dic[key][1][0] == "fp"):
                diff = abs(float(dic[key][0][3])-float(dic[key][1][3]))
                diff_ratio.append(diff)
                if diff == 0: diff_ratio_0 += 1
                elif diff < 0.001: diff_ratio_0001 += 1
                elif diff < 0.002: diff_ratio_0002 += 1
                elif diff < 0.004: diff_ratio_0004 += 1
                elif diff < 0.01: diff_ratio_001 += 1
                elif diff < 0.02: diff_ratio_002 += 1
                elif diff < 0.04: diff_ratio_004 += 1
                elif diff < 0.1: diff_ratio_01 += 1
                elif diff < 0.2: diff_ratio_02 += 1
                elif diff < 0.4: diff_ratio_04 += 1
                elif diff < 1: diff_ratio_1 += 1
                    # print(float(dic[key][0][3]), float(dic[key][1][3]))
                else: diff_else += 1

    print("Unmethylated Score", len(diff_ratio), diff_ratio_0, diff_ratio_0001, diff_ratio_0002, diff_ratio_0004, diff_ratio_001, diff_ratio_002, diff_ratio_004, diff_ratio_01, diff_ratio_02, diff_ratio_04, diff_ratio_1, diff_else)


    # if ( abs(float(fixed_lst[i][3])-float(fp_lst[i][3]))<1 ): continue
    # # print("{} {} {} {:.3f}  {:.3f}  {:.3f}".format(fixed_lst[i][0], fixed_lst[i][1], fixed_lst[i][2], float(fixed_lst[i][3])-float(fp_lst[i][3]), float(fixed_lst[i][4])-float(fp_lst[i][4]), float(fixed_lst[i][5])-float(fp_lst[i][5])))
    # print("{} {} {} {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}".format(fixed_lst[i][0], fixed_lst[i][1], fixed_lst[i][2], float(fixed_lst[i][3]), float(fp_lst[i][3]), float(fixed_lst[i][4]), float(fp_lst[i][4]), float(fixed_lst[i][5]), float(fp_lst[i][5])))




# file = "error_{}.txt".format(data_type)

# f = open(file, "r")
# lines = f.readlines()
# lines = lines[1:-1]

# fp = 0
# approximation = 0
# fp_passing = 0
# approximation_passing = 0
# fp_not_spanned = 0
# approximation_not_spanned = 0
# wrong = 0

# fp_flag = True
# approximation_flag = True
# wrong_both_pass = 0
# wrong_fp_fail = 0
# wrong_approximation_fail = 0

# for line in lines:
#     sep = list(line.split())
#     if len(sep) == 1:
#         wrong += 1
#         if (fp_flag and approximation_flag):
#             wrong_both_pass += 1
#         elif (not fp_flag and approximation_flag):
#             wrong_fp_fail += 1
#         elif (fp_flag and not approximation_flag):
#             wrong_approximation_fail += 1
#     else:
#         if sep[-1] == "fp":
#             fp += 1
#             fp_flag = False
#             if sep[0] == "Pass!":
#                 fp_flag = True
#                 fp_passing += 1
#             elif sep[0] == "Not":
#                 fp_not_spanned += 1
#         elif sep[-1] == "Approximation":
#             approximation += 1
#             approximation_flag = False
#             if sep[0] == "Pass!":
#                 approximation_flag = True
#                 approximation_passing += 1
#             elif sep[0] == "Not":
#                 approximation_not_spanned += 1
# print("fp: {}, approximation: {}, fp_passing: {}, approximation_passing: {}\nfp_not_spanned: {}, approximation_not_spanned: {}\nwrong: {}, wrong_both_pass: {}, wrong_fp_fail: {}, wrong_approximation_fail: {}\n".format(fp, approximation, fp_passing, approximation_passing, fp_not_spanned, approximation_not_spanned, wrong, wrong_both_pass, wrong_fp_fail, wrong_approximation_fail))



# file = "error_{}.txt".format(data_type)
# f = open(file, "r")
# lines = f.readlines()
# lines = lines[1:-1]
# dic = {}
# last_wrong_line = 0
# for line in lines:
#     sep = list(line.split())
#     if len(sep) == 0: continue
#     if len(sep[-1]) < 3:
#         last_wrong_line = int(sep[-1])
#     if last_wrong_line not in dic.keys():
#         dic[last_wrong_line] = 1
#     else:
#         dic[last_wrong_line] += 1
# for i in range(30):
#     if i not in dic.keys(): continue
#     else:
#         print("{} {}".format(i, dic[i]))



# file = "lp_emission.txt"
# f = open(file, "r")
# lines = f.readlines()
# print(lines[0][:-1])
# lines = lines[1:]
# lp_emission_min = 0
# score_min = 0
# lp_emission_sum = 0
# score_sum = 0
# for line in lines[:-1]:
#     item = list(line.split())
#     lp_emission_sum += float(item[0])
#     score_sum += float(item[1])
#     if float(item[0]) < lp_emission_min: lp_emission_min = float(item[0])
#     if float(item[1]) < score_min: score_min = float(item[1])
# print("lp_emission_min: {}, score_min: {}, lp_emission_avg: {}, score_avg: {}".format(lp_emission_min, score_min, lp_emission_min/len(lines), score_sum/len(lines)))



# read = sys.argv[1]
# f = open("score_approximation_{}.txt".format(read), "r")
# f_fp = open("score_fp_{}.txt".format(read), "r")
# f_lines = f.readlines()
# fp_lines = f_fp.readlines()
# f.close()
# f_fp.close()
# dic = {}
# for i in range(len(f_lines)):
#     f_sep = list(f_lines[i].split())
#     fp_sep = list(fp_lines[i].split())
#     if (i % 3) == 0:
#         f_score_d = float(f_sep[1][:-1])
#         score_d = float(fp_sep[1][:-1])
#         dic[i // 3] = [(f_score_d, score_d)]
#     elif (i % 3) == 1:
#         f_score_u = float(f_sep[1][:-1])
#         score_u = float(fp_sep[1][:-1])
#         dic[i // 3] += [(f_score_u, score_u)]
#     elif (i % 3) == 2:
#         f_score_l = float(f_sep[1][:-1])
#         score_l = float(fp_sep[1][:-1])
#         dic[i // 3] += [(f_score_l, score_l)]
# for key in dic.keys():
#     f_score = [dic[key][0][0], dic[key][1][0], dic[key][2][0]]
#     max_index_f = f_score.index(max(f_score))
#     fp_score = [dic[key][0][1], dic[key][1][1], dic[key][2][1]]
#     max_index_fp = fp_score.index(max(fp_score))
#     if (max_index_f == max_index_fp): continue
#     else: print("false", key, f_score, fp_score)



# file = "result_{}.txt".format(data_type)
# f = open(file, "r")
# lines = f.readlines()
# error_trace = []
# similar_trace = []
# not_detected = 0
# over_detected = 0
# for line in lines:
#     lst = list(line.split())
#     if len(lst) == 3:
#         if lst[0] == "error":
#             error_trace.append(lst[-1])
#         elif lst[0] == "similar":
#             similar_trace.append(lst[-1])
# for error in error_trace:
#     if error not in similar_trace:
#         not_detected += 1
# for similar in similar_trace:
#     if similar not in error_trace:
#         over_detected += 1
# f.close()
# print("threshold: 1e-{}".format(data_type), "error_trace: ", len(error_trace), "similar_trace: ", len(similar_trace), "not_detected: ", not_detected, "over_detected: ", over_detected)


check_error_for_all_reads(data_type)

